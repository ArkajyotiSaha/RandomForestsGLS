RFGLS_predict <- function(RFGLS_out, Xtest, h = 1, verbose = FALSE){
  ntest <- nrow(Xtest)
  p <- ncol(Xtest)
  if(missing(Xtest)){stop("error: Xtest must be specified\n")}
  if(!any(is.data.frame(Xtest), is.matrix(Xtest))){stop("error: Xtest must be a data.frame or matrix\n")}
  if(ncol(Xtest) != ncol(RFGLS_out$X)){ stop(paste("error: Xtest must have ",ncol(RFGLS_out$X)," columns\n"))}

  lDaughter <- RFGLS_out$RFGLS_object$ldaughter
  rDaughter <- RFGLS_out$RFGLS_object$rdaughter
  nodestatus <- RFGLS_out$RFGLS_object$nodestatus
  upper <- RFGLS_out$RFGLS_object$upper
  avnode <- RFGLS_out$RFGLS_object$avnode
  mbest <- RFGLS_out$RFGLS_object$mbest
  ntree <- ncol(RFGLS_out$RFGLS_object$ldaughter)

  if(is.null(h)){h <- 4}


  if(h > 1){
    cl <- makeCluster(h)
    clusterExport(cl=cl, varlist=c("Xtest", "ntest", "p", "lDaughter", "rDaughter", "nodestatus", "upper", "avnode", "mbest", "rfglspredict_tree"),envir=environment())
    if(verbose == TRUE){
      cat(paste(("----------------------------------------"), collapse="   "), "\n"); cat(paste(("\tRF Prediction Progress"), collapse="   "), "\n"); cat(paste(("----------------------------------------"), collapse="   "), "\n")
      pboptions(type = "txt", char = "=")
      result <- pblapply(1:ntree,rfglspredict_tree, Xtest, ntest, p, lDaughter, rDaughter, nodestatus,
                         upper, avnode, mbest, cl = cl)
    }
    if(verbose != TRUE){result <- parLapply(cl,1:ntree,rfglspredict_tree, Xtest, ntest, p, lDaughter,
                                            rDaughter, nodestatus, upper, avnode, mbest)}
    stopCluster(cl)
  }

  if(h == 1){
    if(verbose == TRUE){
      cat(paste(("----------------------------------------"), collapse="   "), "\n"); cat(paste(("\tRF Prediction Progress"), collapse="   "), "\n"); cat(paste(("----------------------------------------"), collapse="   "), "\n")
      pboptions(type = "txt", char = "=")
      result <- pblapply(1:ntree,RFGLS_predict_tree, Xtest, ntest, p, lDaughter, rDaughter, nodestatus,
                         upper, avnode, mbest)
    }
    if(verbose != TRUE){result <- lapply(1:ntree,RFGLS_predict_tree, Xtest, ntest, p, lDaughter,
                                            rDaughter, nodestatus, upper, avnode, mbest)}
  }
  result_mat <- do.call(cbind, result)
  RFGLS_prediction <- list()
  RFGLS_prediction$predicted_matrix <- result_mat
  RFGLS_prediction$predicted <- rowMeans(result_mat)
  return(RFGLS_prediction)
}
