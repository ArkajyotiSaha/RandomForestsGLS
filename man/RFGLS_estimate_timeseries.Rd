\name{RFGLS_estimate_timeseries}
\alias{RFGLS_estimate_timeseries}
\title{Function for estimation in time-series data with RF-GLS}

\description{
  The function \code{RFGLS_estimate_spatial} fits univariate non-linear regression models for
  time-series data using a RF-GLS in Saha et al. 2020. \code{RFGLS_estimate_spatial} uses the sparse Cholesky representation
  corresponsinding to \code{AR(q)} process. The fitted Random Forest (RF) model is used later for
  prediction via the \code{RFGLS-predict}.

  Some code blocks are borrowed from the R packages: spNNGP:
  Spatial Regression Models for Large Datasets using Nearest Neighbor
  Gaussian Processes \cr https://CRAN.R-project.org/package=spNNGP and
  randomForest: Breiman and Cutler's Random Forests for Classification
  and Regression \cr https://CRAN.R-project.org/package=randomForest .

}

\usage{
RFGLS_estimate_timeseries(y, X, Xtest = NULL, nrnodes = NULL,
                          nthsize = 20, mtry = 1,
                          pinv_choice = 1, n_omp = 1,
                          ntree = 50, h = 1, lag_params = 0.5,
                          variance = 1,
                          param_estimate = FALSE,
                          verbose = FALSE)
}

\arguments{

  \item{y}{an \eqn{n}{n} length vector of response at the observed time points. }

  \item{X}{an \eqn{n \times p}{n x p} matrix of the covariates in the observation time points. }

  \item{Xtest}{an \eqn{ntest \times p}{ntest x p} matrix of covariates for prediction. Its Structure should be
               identical (including intercept) with that of covariates provided for estimation purpose in \code{X}.
               If \code{NULL}, will use \code{X} as \code{Xtest}. Default value is \code{NULL}. }

  \item{nrnodes}{the maximum number of nodes a tree can have. Default choice leads to the deepest tree contigent on \code{nthsize}. For significantly large \eqn{n},
                 one needs to bound it for growing shallow trees which trades off efficiency for computation time. }

  \item{nthsize}{minimum size of leaf nodes. We recommend not setting this value too small, as that will lead to very deep trees
                 that takes a lot of time to be built and can produce unstable estimaes. Default value is 20. }

  \item{mtry}{number of variables randomly sampled at each partition as a candidate split direction. We recommend using
              the value p/3 where p is the number of variables in \code{X}. Default value is 1. }

  \item{pinv_choice}{dictates the choice of method for obtaining the pseudoinverse involved in the cost function and node
                     representative evaluation. if pinv_choice = 0, SVD is used (slower but more stable), if pinv_choice = 1,
                     orthogonal decomposition (faster, may produce unstable results if \code{nthsize} is too low) is used.
                     Default value is 1. }

  \item{n_omp}{number of threads to be used, value can be more than 1 if source code is compiled with OpenMP support.
               Default is 1. }

  \item{ntree}{number of trees to be grown. This value should not be too small. Default value is 50. }

  \item{h}{number of core to be used in parallel computing setup for bootstrap samples. If h = 1, there is no parallelization.
           Default value is 1. }

  \item{lag_params}{\eqn{q}{q} length vector of AR coefficients. If the parameters need to be estimated from AR(q) process, should be
                    any numeric vector of length q. For notations please see \code{arima}. Default value is 0.5. }

  \item{variance}{variance of the white noise in temporal error. The function estimate is not affected by this. Default value is 1. }

  \item{param_estimate}{if \code{TRUE}, using the residuals obtained from fitting a classical RF default options and \code{nodesize = nthsize},
                        will estimate the coefficeints corresponding to \eqn{AR(q)} from \code{arima} with the option, \code{include.mean = FALSE}.
                        Default value is \code{FALSE}. }

  \item{verbose}{if \code{TRUE}, model specifications along with information regarding OpenMP support and progress of the algorithm
                 is printed to the screen. Otherwise, nothing is printed to the screen. Default value is \code{FALSE}. }
}

\value{
  A list comprising:

  \item{P_matrix}{an \eqn{n \times ntree}{n x ntree} matrix of zero indexed resamples. t-th column denote the
                  \eqn{n} resamples used in the t-th tree. }

  \item{predicted_matrix}{an \eqn{ntest \times ntree}{ntest x ntree} matrix of predictions. t-th column denote the
                          predictions at \eqn{ntest} datapoints obtained from the t-th tree. }

  \item{predicted}{preducted values at the \eqn{ntest} prediction points. Average (\code{rowMeans}) of the treewise predctions
                   in \code{predicted_matrix}, }

  \item{X}{the matrix \code{X}. }

  \item{y}{the vector \code{y}. }

  \item{RFGLS_Object}{object required for prediction. }

}

\references{
  Saha, A., Basu, S., & Datta, A. (2020). Random Forests for dependent data. arXiv preprint arXiv:2007.15421.

  Saha, A., & Datta, A. (2018). BRISC: bootstrap for rapid inference on spatial
  covariances. Stat, e184, DOI: 10.1002/sta4.184.

  Andy Liaw, and Matthew Wiener (2015). randomForest: Breiman and Cutler's Random
  Forests for Classification and Regression. R package version 4.6-14. \cr
  https://CRAN.R-project.org/package=randomForest

  Andrew Finley, Abhirup Datta and Sudipto Banerjee (2017). spNNGP: Spatial Regression Models for Large
  Datasets using Nearest Neighbor Gaussian Processes. R package version 0.1.1.
  https://CRAN.R-project.org/package=spNNGP

}

\author{
  Arkajyoti Saha \email{arkajyotisaha93@gmail.com}, \cr
  Sumanta Basu \email{sumbose@cornell.edu}, \cr
  Abhirup Datta \email{abhidatta@jhu.edu}
}

\examples{

set.seed(2)
n <- 200
x <- as.matrix(rnorm(n),n,1)

sigma.sq <- 1
rho <- 0.5

set.seed(3)
b <- rho
s <- sqrt(sigma.sq)
eps = arima.sim(list(order = c(1,0,0), ar = b),
                n = n, rand.gen = rnorm, sd = s)

y <- eps + 10*sin(pi * x)

estimation_result <- RFGLS_estimate_timeseries(y, x, ntree = 10)

}

\keyword{model}

