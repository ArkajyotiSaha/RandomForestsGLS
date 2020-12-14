\name{RFGLS_predict_spatial}
\alias{RFGLS_predict_spatial}
\title{Prediction at new location with RF-GLS}

\description{
  The function \code{RFGLS_predict_spatial} performs fast prediction on a set of new locations with
  RF-GLS models using Nearest Neighbor Gaussian Processes (NNGP) (Datta et al., 2016).

  Some code blocks are borrowed from the R packages: spNNGP:
  Spatial Regression Models for Large Datasets using Nearest Neighbor Gaussian Processes
  https://CRAN.R-project.org/package=spNNGP and randomForest: Breiman and Cutler's Random
  Forests for Classification and Regression https://CRAN.R-project.org/package=randomForest .

}

\usage{
RFGLS_predict_spatial <- function(RFGLS_out, coords.0, Xtest, h = 1, verbose = FALSE)
}

\arguments{

  \item{RFGLS_out}{An object obtained from \code{RFGLS_estimate_spatial} or \code{RFGLS_estimate_timeseries}. }

  \item{coords.0}{the spatial coordinates corresponding to prediction locations. Its structure should be same as that of coords
                  in \code{BRISC_estimation}. Default value is a column of \eqn{1} to adjust for the mean(intercept). }

  \item{Xtest}{an \eqn{ntest \times p}{ntest x p} matrix of covariates for prediction. Its Structure should be
               identical (including intercept) with that of covariates provided for estimation purpose in \code{X}
               in RFGLS_out. }

  \item{h}{Number of core to be used in parallel computing setup for bootstrap samples. If h = 1, there is no parallelization.
           Default value is 1. }

  \item{verbose}{if TRUE, model specifications along with information regarding OpenMP support and progress of the algorithm
                 is printed to the screen. Otherwise, nothing is printed to the screen. Default value is \code{FALSE}. }
}

\value{
  A list comprising:

  \item{prediction}{predicted response corresponding to X.0 and coords.0. }

}

\references{
  Saha, A., & Datta, A. (2018). BRISC: bootstrap for rapid inference on spatial
  covariances. Stat, e184, DOI: 10.1002/sta4.184.

  Datta, A., S. Banerjee, A.O. Finley, and A.E. Gelfand. (2016)
  Hierarchical Nearest-Neighbor Gaussian process models for large
  geostatistical datasets. Journal of the American Statistical
  Association, 111:800-812.

  Andrew Finley, Abhirup Datta and Sudipto Banerjee (2017). spNNGP: Spatial Regression Models for Large
  Datasets using Nearest Neighbor Gaussian Processes. R package version 0.1.1.
  https://CRAN.R-project.org/package=spNNGP

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

set.seed(1)
n <- 250
coords <- cbind(runif(n,0,1), runif(n,0,1))

set.seed(2)
x <- as.matrix(rnorm(n),n,1)

sigma.sq = 1
phi = 5
tau.sq = 0.1

D <- as.matrix(dist(coords))
R <- exp(-phi*D)
w <- rmvn(1, rep(0,n), sigma.sq*R)

y <- rnorm(n, 10*sin(pi * x) + w, sqrt(tau.sq))

estimation_result <- RFGLS_estimate_spatial(coords[1:200,], y[1:200], matrix(x[1:200,],200,1), ntree = 10)
prediction_result <- RFGLS_predict_spatial(estimation_result, coords[201:250,], matrix(x[201:250,],50,1))

}

\keyword{model}

