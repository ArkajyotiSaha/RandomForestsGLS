RFGLS_estimate_timeseries <- function(y, X, Xtest = NULL, nrnodes = NULL, nthsize = 20, mtry = 1, pinv_choice = 1, n_omp = 1,
                                    ntree = 50, h = 1, lag_params = 0.5, variance = 1, param_estimate = FALSE, verbose = FALSE){

  n <- length(y)
  nsample <- n
  if(is.null(nrnodes)){
    nrnodes <- 2 * nsample + 1
  }
  if(param_estimate){
    sp <- randomForest(X, y)
    sp_input_est <- predict(sp, X)
    rf_residual <- y - sp_input_est
    AR <- arima(rf_residual, order = c(length(lag_params),0, 0), include.mean = FALSE)
    lag_params <- AR$coef
  }

  ##Option for Multithreading if compiled with OpenMp support
  n.omp.threads <- as.integer(n_omp)
  storage.mode(n.omp.threads) <- "integer"

  ##type conversion
  storage.mode(n) <- "integer"
  storage.mode(verbose) <- "integer"

  q <- length(lag_params)
  res_BF <- list()
  res_BF$B <- c(rep(0, sum(1:(q-1))), rep(lag_params, n-q))
  res_BF$F <- rep(sigma.sq, n)
  res_BF$nnIndxLU <- rep(0, 2*n)
  res_BF$nnIndxLU[(n+1):(n+q)] <- (0:(q-1))
  res_BF$nnIndxLU[(n+q+1):(2*n)] <- q
  res_BF$nnIndxLU[2:n] <- cumsum(res_BF$nnIndxLU[(n+1):((2*n)-1)])
  res_BF$nnIndx <- unlist(sapply(1:(n-1), function(i) (i-1):max(0, (i - q))))
  res_Z <- .Call("RFGLS_invZcpp", as.integer(length(res_BF$nnIndxLU)/2), as.integer(res_BF$nnIndx), as.integer(res_BF$nnIndxLU), as.integer(rep(0, length(res_BF$nnIndxLU)/2)), as.integer(0*res_BF$nnIndx), as.integer(rep(0, length(res_BF$nnIndxLU)/2 + 1)), as.integer(rep(0, length(res_BF$nnIndxLU)/2)) )

  p <- ncol(X)
  storage.mode(p) <- "integer"
  storage.mode(nsample) <- "integer"
  #nthsize <- 20
  storage.mode(nthsize) <- "integer"

  storage.mode(nrnodes) <- "integer"
  #mtry <- 1
  storage.mode(mtry) <- "integer"
  treeSize <- 0
  storage.mode(treeSize) <- "integer"
  #pinv_choice <- 0
  storage.mode(pinv_choice) <- "integer"
  if(is.null(Xtest)){
    Xtest <- X
  }
  ntest <- nrow(Xtest)
  storage.mode(ntest) <- "integer"
  if(is.null(h)){h <- 4}

  q_lag <- length(lag_params)
  storage.mode(q_lag) <- "integer"


  #out_rtree <- .Call("RFGLStree_treecpp", t(X), y, res_BF$B, res_BF$F, as.integer(res_BF$nnIndx), as.integer(res_BF$nnIndxLU), as.integer(res_Z$invZ_val), as.integer(res_Z$invZ_loc), mtry, n, p, nsample, nthsize, nrnodes, treeSize, pinv_choice, t(Xtest), ntest, n.omp.threads, q_lag)
  #return(out_rtree)
  if(h > 1){
    cl <- makeCluster(h)
    clusterExport(cl=cl, varlist=c("X", "y", "res_BF", "res_Z", "mtry", "n", "p",
                                   "nsample", "nthsize", "nrnodes", "treeSize", "pinv_choice", "Xtest", "ntest",
                                   "n.omp.threads", "RFGLS_tree", "q"),envir=environment())
    if(verbose == TRUE){
      cat(paste(("----------------------------------------"), collapse="   "), "\n"); cat(paste(("\tRF Progress"), collapse="   "), "\n"); cat(paste(("----------------------------------------"), collapse="   "), "\n")
      pboptions(type = "txt", char = "=")
      result <- pblapply(1:ntree,RFGLS_tree, X, y, res_BF, res_Z, mtry, n, p,
                         nsample, nthsize, nrnodes, treeSize, pinv_choice, Xtest, ntest,
                         n.omp.threads, q, cl = cl)
    }
    if(verbose != TRUE){result <- parLapply(cl,1:ntree,RFGLS_tree, X, y, res_BF, res_Z, mtry, n, p,
                                            nsample, nthsize, nrnodes, treeSize, pinv_choice, Xtest, ntest,
                                            n.omp.threads, q)}
    stopCluster(cl)
  }
  if(h == 1){
    if(verbose == TRUE){
      cat(paste(("----------------------------------------"), collapse="   "), "\n"); cat(paste(("\tRF Progress"), collapse="   "), "\n"); cat(paste(("----------------------------------------"), collapse="   "), "\n")
      pboptions(type = "txt", char = "=")
      result <- pblapply(1:ntree,RFGLS_tree, X, y, res_BF, res_Z, mtry, n, p,
                         nsample, nthsize, nrnodes, treeSize, pinv_choice, Xtest, ntest,
                         n.omp.threads, q)
    }

    if(verbose != TRUE){
      result <- lapply(1:ntree,RFGLS_tree, X, y, res_BF, res_Z, mtry, n, p,
                       nsample, nthsize, nrnodes, treeSize, pinv_choice, Xtest, ntest,
                       n.omp.threads, q)
    }
  }
  #result_mat <- do.call(cbind, result)
  RFGLS_out <- list()
  RFGLS_out$P_matrix <- do.call(cbind, lapply(1:ntree, function(i) result[[i]]$P_index))
  RFGLS_out$predicted_matrix <- do.call(cbind, lapply(1:ntree, function(i) result[[i]]$ytest))
  RFGLS_out$predicted <- rowMeans(RFGLS_out$predicted_matrix)
  RFGLS_out$X <- X
  RFGLS_out$y <- y
  RFGLS_out$RFGLS_object <- list()
  RFGLS_out$RFGLS_object$ldaughter <- do.call(cbind, lapply(1:ntree, function(i) result[[i]]$lDaughter))
  RFGLS_out$RFGLS_object$rdaughter <- do.call(cbind, lapply(1:ntree, function(i) result[[i]]$rDaughter))
  RFGLS_out$RFGLS_object$nodestatus <- do.call(cbind, lapply(1:ntree, function(i) result[[i]]$nodestatus))
  RFGLS_out$RFGLS_object$upper <- do.call(cbind, lapply(1:ntree, function(i) result[[i]]$upper))
  RFGLS_out$RFGLS_object$avnode <- do.call(cbind, lapply(1:ntree, function(i) result[[i]]$avnode))
  RFGLS_out$RFGLS_object$mbest <- do.call(cbind, lapply(1:ntree, function(i) result[[i]]$mbest))
  return(RFGLS_out)
}
