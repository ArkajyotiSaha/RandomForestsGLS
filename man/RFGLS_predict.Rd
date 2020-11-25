\name{RFGLS_predict}
\alias{RFGLS_predict}
\title{Prediction of function value with RF-GLS}

\description{
  The function \code{RFGLS_predict} predicts the functional value at a given set of covariates.
  It uses a fitted RF-GLS model to obtain the predictions.

  Some code blocks are borrowed from the R package: randomForest: Breiman and Cutler's Random
  Forests for Classification and Regression https://CRAN.R-project.org/package=randomForest .

}

\usage{
RFGLS_predict  <- function(RFGLS_out, Xtest, h = 1, verbose = FALSE)
}

\arguments{

  \item{RFGLS_out}{An object obtained from \code{RFGLS_estimate_spatial} or \code{RFGLS_estimate_timeseries}. }

  \item{Xtest}{an \eqn{ntest \times p}{ntest x p} matrix of covariates for prediction. Its Structure should be
               identical (including intercept) with that of covariates provided for estimation purpose in \code{X}
               in RFGLS_out. }

  \item{h}{Number of core to be used in parallel computing setup for bootstrap samples. If h = 1, there is no parallelization.
           Default value is 1. }

  \item{verbose}{if TRUE, model specifications along with information regarding OpenMP support and progress of the algorithm
                 is printed to the screen. Otherwise, nothing is printed to the screen. Default value is FALSE. }
}

\value{
  A list comprising:

  \item{predicted_matrix}{an \eqn{ntest \times ntree}{n x ntree} matrix of predictions. \code{t}-th column denote the
                          predictions at \code{ntest} datapoints obtained from the \code{t}-th tree. }

  \item{predicted}{preducted values at the \code{ntest} prediction points. Average (\code{rowMeans}) of the treewise predctions
                   in \code{predicted_matrix} }

}

\references{
  Andy Liaw, and Matthew Wiener (2015). randomForest: Breiman and Cutler's Random
  Forests for Classification and Regression. R package version 4.6-14.
  https://CRAN.R-project.org/package=randomForest

}

\author{
  Arkajyoti Saha \email{arkajyotisaha93@gmail.com}, \cr
  Abhirup Datta \email{abhidatta@jhu.edu}
}

\examples{

rmvn <- function(n, mu = 0, V = matrix(1)){
  p <- length(mu)
  if(any(is.na(match(dim(V),p))))
    stop("Dimension not right!")
  D <- chol(V)
  t(matrix(rnorm(n*p), ncol=p)\%*\%D + rep(mu,rep(n,p)))
}

set.seed(2)
n <- 200
x <- as.matrix(rnorm(n),n,1)

sigma.sq <- 1
rho <- 0.5

set.seed(3)
b <- rho
s <- sqrt(sigma.sq)
eps = arima.sim(list(order = c(1,0,0), ar = b), n = n, rand.gen = rnorm, sd = s)

y <- eps + 10*sin(pi * x[,1])

estimation_result <- RFGLS_estimate_timeseries(y, x, ntree = 10)
Xtest <- matrix(seq(0,1, by = 1/1000), 1001, 1)
RFGLS_predict <- RFGLS_predict(estimation_result, Xtest)

}

\keyword{model}

