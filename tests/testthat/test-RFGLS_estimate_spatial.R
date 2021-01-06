test_that("estimate_spatial works", {
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
  set.seed(1)
  est_known <- RFGLS_estimate_spatial(coords, y, x, ntree = 10, cov.model = "exponential",
                                      nthsize = 20, sigma.sq = sigma.sq, tau.sq = tau.sq,
                                      phi = phi)
  expect_true(is.matrix(est_known$predicted_matrix))
  expect_equal(dim(est_known$predicted_matrix), c(200,10))
  expect_length(est_known$predicted, 200)
  set.seed(1)
  est_known_parallel <- RFGLS_estimate_spatial(coords, y, x, ntree = 10, cov.model = "exponential",
                                      nthsize = 20, sigma.sq = sigma.sq, tau.sq = tau.sq,
                                      phi = phi, h = 2)
  expect_equal(est_known_parallel$P_matrix, est_known$P_matrix)
})
