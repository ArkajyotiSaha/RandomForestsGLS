RFGLS_estimate_spatial <- function(coords, y, X, Xtest = NULL, nrnodes = NULL, nthsize = 20, mtry = 1, pinv_choice = 1, n_omp = 1, ntree = 50, h = 1,
                                   sigma.sq = 1, tau.sq = 0.1, phi = 5, nu = 0.5, n.neighbors = 15, cov.model = "exponential", search.type = "tree",
                                   param_estimate = FALSE, verbose = FALSE){

  n <- nrow(coords)
  nsample <- n
  if(is.null(nrnodes)){
    nrnodes <- 2 * nsample + 1
  }

  if(is.null(Xtest)){
    Xtest <- X
  }
  if(ncol(Xtest) != ncol(X)){ stop(paste("error: Xtest must have ",ncol(X)," columns\n"))}

  if(param_estimate){
    sp <- randomForest(X, y, nodesize = nthsize)
    sp_input_est <- predict(sp, X)
    rf_residual <- y - sp_input_est
    if(verbose){
      cat(paste(("----------------------------------------"), collapse="   "), "\n"); cat(paste(("\tParameter Estimation"), collapse="   "), "\n"); cat(paste(("----------------------------------------"), collapse="   "), "\n")
    }
    est_theta <- BRISC_estimation(coords, x = matrix(1,n,1), y = rf_residual, verbose = verbose, cov.model = cov.model)
    sigma.sq <- est_theta$Theta[1]
    tau.sq <- est_theta$Theta[2]
    phi <- est_theta$Theta[3]
    if(cov.model =="matern"){
      nu <- est_theta$Theta[4]
    }
  }

  cov.model.names <- c("exponential","spherical","matern","gaussian")
  cov.model.indx <- which(cov.model == cov.model.names) - 1
  storage.mode(cov.model.indx) <- "integer"

  ##Parameter values
  if(cov.model!="matern"){
    initiate <- c(sigma.sq, tau.sq, phi)
    names(initiate) <- c("sigma.sq", "tau.sq", "phi")
  }
  else{
    initiate <- c(sigma.sq, tau.sq, phi, nu)
    names(initiate) <- c("sigma.sq", "tau.sq", "phi", "nu")}

  alpha.sq.starting <- sqrt(tau.sq/sigma.sq)
  phi.starting <- sqrt(phi)
  nu.starting <- sqrt(nu)

  storage.mode(alpha.sq.starting) <- "double"
  storage.mode(phi.starting) <- "double"
  storage.mode(nu.starting) <- "double"

  search.type.names <- c("brute", "tree")
  if(!search.type %in% search.type.names){
    stop("error: specified search.type '",search.type,"' is not a valid option; choose from ", paste(search.type.names, collapse=", ", sep="") ,".")
  }
  search.type.indx <- which(search.type == search.type.names)-1
  storage.mode(search.type.indx) <- "integer"


  ##Option for Multithreading if compiled with OpenMp support
  n.omp.threads <- as.integer(n_omp)
  storage.mode(n.omp.threads) <- "integer"

  fix_nugget <- 1
  ##type conversion
  storage.mode(n) <- "integer"
  storage.mode(coords) <- "double"
  storage.mode(n.neighbors) <- "integer"
  storage.mode(verbose) <- "integer"

  if(verbose){
    cat(paste(("----------------------------------------"), collapse="   "), "\n"); cat(paste(("\tRFGLS Model Fitting"), collapse="   "), "\n"); cat(paste(("----------------------------------------"), collapse="   "), "\n")
  }

  res_BF <- .Call("RFGLS_BFcpp", n, n.neighbors, coords, cov.model.indx, alpha.sq.starting, phi.starting, nu.starting, search.type.indx, n.omp.threads, verbose)
  res_Z <- .Call("RFGLS_invZcpp", as.integer(length(res_BF$nnIndxLU)/2), as.integer(res_BF$nnIndx), as.integer(res_BF$nnIndxLU), as.integer(rep(0, length(res_BF$nnIndxLU)/2)), as.integer(0*res_BF$nnIndx), as.integer(rep(0, length(res_BF$nnIndxLU)/2 + 1)), as.integer(rep(0, length(res_BF$nnIndxLU)/2)) )

  p <- ncol(X)
  storage.mode(p) <- "integer"
  storage.mode(nsample) <- "integer"

  storage.mode(nthsize) <- "integer"
  if(is.null(nrnodes)){
    nrnodes <- 2 * nsample + 1
  }
  storage.mode(nrnodes) <- "integer"

  storage.mode(mtry) <- "integer"
  treeSize <- 0
  storage.mode(treeSize) <- "integer"

  storage.mode(pinv_choice) <- "integer"

  ntest <- nrow(Xtest)
  storage.mode(ntest) <- "integer"
  if(is.null(h)){h <- 1}

  q <- 0
  storage.mode(q) <- "integer"

  local_seed <- sample(.Random.seed, 1)


  if(h > 1){
    cl <- makeCluster(h)
    clusterExport(cl=cl, varlist=c("X", "y", "res_BF", "res_Z", "mtry", "n", "p",
                                   "nsample", "nthsize", "nrnodes", "treeSize", "pinv_choice", "Xtest", "ntest",
                                   "n.omp.threads", "RFGLS_tree", "q", "local_seed"),envir=environment())
    if(verbose == TRUE){
      cat(paste(("----------------------------------------"), collapse="   "), "\n"); cat(paste(("\tRF Progress"), collapse="   "), "\n"); cat(paste(("----------------------------------------"), collapse="   "), "\n")
      pboptions(type = "txt", char = "=")
      result <- pblapply(1:ntree,RFGLS_tree, X, y, res_BF, res_Z, mtry, n, p,
                         nsample, nthsize, nrnodes, treeSize, pinv_choice, Xtest, ntest,
                         n.omp.threads, q, local_seed, cl = cl)
    }
    if(verbose != TRUE){result <- parLapply(cl,1:ntree,RFGLS_tree, X, y, res_BF, res_Z, mtry, n, p,
                                            nsample, nthsize, nrnodes, treeSize, pinv_choice, Xtest, ntest,
                                            n.omp.threads, q, local_seed)}
    stopCluster(cl)
  }
  if(h == 1){
    if(verbose == TRUE){
      cat(paste(("----------------------------------------"), collapse="   "), "\n"); cat(paste(("\tRF Progress"), collapse="   "), "\n"); cat(paste(("----------------------------------------"), collapse="   "), "\n")
      pboptions(type = "txt", char = "=")
      result <- pblapply(1:ntree,RFGLS_tree, X, y, res_BF, res_Z, mtry, n, p,
                         nsample, nthsize, nrnodes, treeSize, pinv_choice, Xtest, ntest,
                         n.omp.threads, q, local_seed)
    }

    if(verbose != TRUE){
      result <- lapply(1:ntree,RFGLS_tree, X, y, res_BF, res_Z, mtry, n, p,
                       nsample, nthsize, nrnodes, treeSize, pinv_choice, Xtest, ntest,
                       n.omp.threads, q, local_seed)
    }
  }

  RFGLS_out <- list()
  RFGLS_out$P_matrix <- do.call(cbind, lapply(1:ntree, function(i) result[[i]]$P_index))
  RFGLS_out$predicted_matrix <- do.call(cbind, lapply(1:ntree, function(i) result[[i]]$ytest))
  RFGLS_out$predicted <- rowMeans(RFGLS_out$predicted_matrix)
  RFGLS_out$X <- X
  RFGLS_out$y <- y
  RFGLS_out$coords <- coords
  RFGLS_out$RFGLS_object <- list()
  RFGLS_out$RFGLS_object$ldaughter <- do.call(cbind, lapply(1:ntree, function(i) result[[i]]$lDaughter))
  RFGLS_out$RFGLS_object$rdaughter <- do.call(cbind, lapply(1:ntree, function(i) result[[i]]$rDaughter))
  RFGLS_out$RFGLS_object$nodestatus <- do.call(cbind, lapply(1:ntree, function(i) result[[i]]$nodestatus))
  RFGLS_out$RFGLS_object$upper <- do.call(cbind, lapply(1:ntree, function(i) result[[i]]$upper))
  RFGLS_out$RFGLS_object$avnode <- do.call(cbind, lapply(1:ntree, function(i) result[[i]]$avnode))
  RFGLS_out$RFGLS_object$mbest <- do.call(cbind, lapply(1:ntree, function(i) result[[i]]$mbest))
  return(RFGLS_out)
}
