RFGLS_tree <- function(i, X, y, res_BF, res_Z, mtry, n, p,
                       nsample, nthsize, nrnodes, treeSize, pinv_choice, Xtest, ntest,
                       n.omp.threads, q){
  set.seed(i)
  rfgls_stat <- .Call("RFGLStree_cpp", t(X), y, res_BF$B, res_BF$F, as.integer(res_BF$nnIndx), as.integer(res_BF$nnIndxLU), as.integer(res_Z$invZ_val), as.integer(res_Z$invZ_loc), mtry, n, p, nsample, nthsize, nrnodes, treeSize, pinv_choice, t(Xtest), ntest, n.omp.threads, q)
  result <- rfgls_stat
  result
}

RFGLS_predict_tree <- function(i, Xtest, ntest, p, lDaughter, rDaughter, nodestatus, upper, avnode, mbest){
  rfgls_predict <- .Call("RFGLSpredicttree_cpp", t(Xtest), as.integer(ntest), as.integer(p),
                         as.integer(lDaughter[,i]), as.integer(rDaughter[,i]), as.integer(nodestatus[,i]), upper[,i], avnode[,i], as.integer(mbest[,i]))
  return(rfgls_predict$ytest)
}
