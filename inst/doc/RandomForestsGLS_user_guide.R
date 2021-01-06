## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(RandomForestsGLS)

## ----simulation, message=FALSE, warning=FALSE---------------------------------
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

## ----model_fit, message=FALSE, warning=FALSE----------------------------------
set.seed(1)
est_known <- RFGLS_estimate_spatial(coords, y, x, ntree = 50, cov.model = "exponential",
                                    nthsize = 20, sigma.sq = sigma.sq, tau.sq = tau.sq, 
                                    phi = phi)

## ----model_fit_unknown, message=FALSE, warning=FALSE--------------------------
set.seed(1)
est_unknown <- RFGLS_estimate_spatial(coords, y, x, ntree = 50, cov.model = "exponential", 
                                      nthsize = 20, param_estimate = TRUE)

## ----model_estimate, message=FALSE, warning=FALSE-----------------------------
Xtest <- matrix(seq(0,1, by = 1/10000), 10001, 1)
RFGLS_predict_known <- RFGLS_predict(est_known, Xtest)

## ----comparison, message=FALSE, warning=FALSE---------------------------------
library(randomForest)
set.seed(1)
RF_est <- randomForest(x, y, nodesize = 20)

RF_predict <- predict(RF_est, Xtest)
#RF MISE
mean((RF_predict - 10*sin(pi * Xtest))^2)

#RF-GLS MISE
mean((RFGLS_predict_known$predicted - 10*sin(pi * Xtest))^2)

RFGLS_predict_unknown <- RFGLS_predict(est_unknown, Xtest)
#RF-GLS unknown MISE
mean((RFGLS_predict_unknown$predicted - 10*sin(pi * Xtest))^2)

## ----plot_comparison, message=FALSE, warning=FALSE----------------------------
rfgls_loess_10 <- loess(RFGLS_predict_known$predicted ~ c(1:length(Xtest)), span=0.1)
rfgls_smoothed10 <- predict(rfgls_loess_10)

rf_loess_10 <- loess(RF_predict ~ c(1:length(RF_predict)), span=0.1)
rf_smoothed10 <- predict(rf_loess_10)

xval <- c(10*sin(pi * Xtest), rf_smoothed10, rfgls_smoothed10)
xval_tag <- c(rep("Truth", length(10*sin(pi * Xtest))), rep("RF", length(rf_smoothed10)), 
              rep("RF-GLS",length(rfgls_smoothed10)))
plot_data <- as.data.frame(xval)
plot_data$Methods <- xval_tag
coval <- c(rep(seq(0,1, by = 1/10000), 3))
plot_data$Covariate <- coval

library(ggplot2)
ggplot(plot_data, aes(x=Covariate, y=xval, color=Methods)) +
geom_point() + labs( x = "x") + labs( y = "f(x)")

## ----prediction_spatial, message=FALSE, warning=FALSE-------------------------
est_known_short <- RFGLS_estimate_spatial(coords[1:160,], y[1:160], 
                   matrix(x[1:160,],160,1), ntree = 50,  cov.model = "exponential", 
                   nthsize = 20, param_estimate = TRUE)
RFGLS_predict_spatial <- RFGLS_predict_spatial(est_known_short, coords[161:200,], 
                                               matrix(x[161:200,],40,1))
pred_mat <- as.data.frame(cbind(RFGLS_predict_spatial$prediction, y[161:200]))
colnames(pred_mat) <- c("Predicted", "Observed")
ggplot(pred_mat, aes(x=Observed, y=Predicted)) + geom_point() + 
  geom_abline(intercept = 0, slope = 1, color = "blue") +
  ylim(0, 16) + xlim(0, 16)

## ----misspec_spatial, message=FALSE, warning=FALSE----------------------------
#Data simulation from matern with nu = 1.5
nu = 3/2
R1 <- (D*phi)^nu/(2^(nu-1)*gamma(nu))*besselK(x=D*phi, nu=nu)
diag(R1) <- 1
set.seed(2)
w <- rmvn(1, rep(0,n), sigma.sq*R1)
y <- rnorm(n, 10*sin(pi * x) + w, sqrt(tau.sq))

#RF-GLS with exponential covariance
set.seed(3)
est_misspec <- RFGLS_estimate_spatial(coords, y, x, ntree = 50, cov.model = "exponential", 
                                      nthsize = 20, param_estimate = TRUE)
RFGLS_predict_misspec <- RFGLS_predict(est_misspec, Xtest)

#RF
set.seed(4)
RF_est <- randomForest(x, y,  nodesize = 20)
RF_predict <- predict(RF_est, Xtest)

#RF-GLS MISE
mean((RFGLS_predict_misspec$predicted - 10*sin(pi * Xtest))^2)
#RF MISE
mean((RF_predict - 10*sin(pi * Xtest))^2)

## ----simulation_temporal------------------------------------------------------
rho <- 0.9
set.seed(1)
b <- rho
s <- sqrt(sigma.sq)
eps = arima.sim(list(order = c(1,0,0), ar = b), n = n, rand.gen = rnorm, sd = s)
y <- c(eps + 10*sin(pi * x))

## ----model_fit_temporal, message=FALSE, warning=FALSE-------------------------
set.seed(1)
est_temp_known <- RFGLS_estimate_timeseries(y, x, ntree = 50, lag_params = rho, nthsize = 20)

## ----model_fit_temporal_unknown, message=FALSE, warning=FALSE-----------------
set.seed(1)
est_temp_unknown <- RFGLS_estimate_timeseries(y, x, ntree = 50, lag_params = rho, 
                                              nthsize = 20, param_estimate = TRUE)

## ----model_estimate_temporal, message=FALSE, warning=FALSE--------------------
Xtest <- matrix(seq(0,1, by = 1/10000), 10001, 1)
RFGLS_predict_temp_known <- RFGLS_predict(est_temp_known, Xtest)

## ----comparison_temp, message=FALSE, warning=FALSE----------------------------
library(randomForest)
set.seed(1)
RF_est_temp <- randomForest(x, y, nodesize = 20)

RF_predict_temp <- predict(RF_est_temp, Xtest)
#RF MISE
mean((RF_predict_temp - 10*sin(pi * Xtest))^2)

#RF-GLS MISE
mean((RFGLS_predict_temp_known$predicted - 10*sin(pi * Xtest))^2)

RFGLS_predict_temp_unknown <- RFGLS_predict(est_temp_unknown, Xtest)
#RF-GLS unknown MISE
mean((RFGLS_predict_temp_unknown$predicted - 10*sin(pi * Xtest))^2)

## ----misspec_temporal, message=FALSE, warning=FALSE---------------------------
#Simulation from AR(2) process
rho1 <- 0.7
rho2 <- 0.2
set.seed(2)
b <- c(rho1, rho2)
s <- sqrt(sigma.sq)
eps = arima.sim(list(order = c(2,0,0), ar = b), n = n, rand.gen = rnorm, sd = s)
y <- c(eps + 10*sin(pi * x))

#RF-GLS with AR(1)
set.seed(3)
est_misspec_temp <- RFGLS_estimate_timeseries(y, x, ntree = 50, lag_params = 0, 
                                              nthsize = 20, param_estimate = TRUE)
RFGLS_predict_misspec_temp <- RFGLS_predict(est_misspec_temp, Xtest)

#RF
set.seed(4)
RF_est_temp <- randomForest(x, y,  nodesize = 20)
RF_predict_temp <- predict(RF_est_temp, Xtest)

#RF-GLS MISE
mean((RFGLS_predict_misspec_temp$predicted - 10*sin(pi * Xtest))^2)
#RF MISE
mean((RF_predict_temp - 10*sin(pi * Xtest))^2)

## ----parallel_spatial, message=FALSE, warning=FALSE---------------------------
#simulation from exponential distribution
set.seed(5)
n <- 200
coords <- cbind(runif(n,0,1), runif(n,0,1))
set.seed(2)
x <- as.matrix(runif(n),n,1)
sigma.sq = 10
phi = 1
tau.sq = 0.1
nu = 0.5
D <- as.matrix(dist(coords))
R <- exp(-phi*D)
w <- rmvn(1, rep(0,n), sigma.sq*R)
y <- rnorm(n, 10*sin(pi * x) + w, sqrt(tau.sq))

#RF-GLS model fitting and prediction with parallel computation
set.seed(1)
est_known_pl <- RFGLS_estimate_spatial(coords, y, x, ntree = 50, cov.model = "exponential",
                                    nthsize = 20, sigma.sq = sigma.sq, tau.sq = tau.sq, 
                                    phi = phi, h = 2)
RFGLS_predict_known_pl <- RFGLS_predict(est_known_pl, Xtest, h = 2)

#MISE from single core
mean((RFGLS_predict_known$predicted - 10*sin(pi * Xtest))^2)
#MISE from parallel computation
mean((RFGLS_predict_known_pl$predicted - 10*sin(pi * Xtest))^2)

