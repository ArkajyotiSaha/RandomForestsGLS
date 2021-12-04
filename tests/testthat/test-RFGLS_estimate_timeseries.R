test_that("estimate_timeseries works", {
  n <- 200
  set.seed(2)
  x <- as.matrix(runif(n),n,1)
  sigma.sq = 10
  rho <- 0.9
  set.seed(1)
  b <- rho
  s <- sqrt(sigma.sq)
  eps = arima.sim(list(order = c(1,0,0), ar = b), n = n, rand.gen = rnorm, sd = s)
  y <- c(eps + 10*sin(pi * x))
  set.seed(1)
  est_temp_known <- RFGLS_estimate_timeseries(y, x, ntree = 10, lag_params = rho, nthsize = 20)
  expect_true(is.matrix(est_temp_known$predicted_matrix))
  expect_equal(dim(est_temp_known$predicted_matrix), c(200,10))
  expect_length(est_temp_known$predicted, 200)
  round(est_temp_known$predicted[c(1,200)],3) == c(8.317, 6.183)
  set.seed(1)
  est_temp_known_parallel <- RFGLS_estimate_timeseries(y, x, ntree = 10, lag_params = rho, nthsize = 20, h = 2)
  expect_equal( est_temp_known_parallel$P_matrix, est_temp_known$P_matrix)
})
