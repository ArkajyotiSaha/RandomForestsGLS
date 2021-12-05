test_that("predict_spatial works", {
  rmvn <- function(n, mu = 0, V = matrix(1)){
    p <- length(mu)
    if(any(is.na(match(dim(V),p))))
      stop("Dimension not right!")
    D <- chol(V)
    t(matrix(rnorm(n*p), ncol=p)%*%D + rep(mu,rep(n,p)))
  }

  set.seed(5)
  n <- 200
  coords <- cbind(runif(n,0,1), runif(n,0,1))
  set.seed(2)
  x <- as.matrix(runif(n),n,1)
  sigma.sq = 10
  phi = 1
  tau.sq = 0.1
  D <- as.matrix(dist(coords))
  R <- exp(-phi*D)
  w <- rmvn(1, rep(0,n), sigma.sq*R)
  y <- rnorm(n, 10*sin(pi * x) + w, sqrt(tau.sq))
  est_known_short <- RFGLS_estimate_spatial(coords[1:160,], y[1:160],
                                            matrix(x[1:160,],160,1), ntree = 10,  cov.model = "exponential",
                                            nthsize = 20, param_estimate = TRUE)
  RFGLS_predict_spatial <- RFGLS_predict_spatial(est_known_short, coords[161:200,],
                                                 matrix(x[161:200,],40,1))
  expect_length(RFGLS_predict_spatial$prediction, 40)
  expect_equal(round(RFGLS_predict_spatial$prediction[c(1, 40)], 3), c(5.795, 7.346))
})
