\name{RFGLS_estimate_spatial}
\alias{RFGLS_estimate_spatial}
\title{Function for estimation in spatial data with RF-GLS}

\description{
  The function \code{RFGLS_estimate_spatial} fits univariate non-linear spatial regression models for
  spatial data using RF-GLS in Saha et al. 2020. \code{RFGLS_estimate_spatial} uses the sparse Cholesky representation
  of Vecchia’s likelihood (Vecchia, 1988) developed in Datta et al., 2016 and Saha & Datta, 2018.
  The fitted Random Forest (RF) model is used later for prediction via the \code{RFGLS_predict} and \code{RFGLS_predict_spatial}.

  Some code blocks are borrowed from the R packages: spNNGP:
  Spatial Regression Models for Large Datasets using Nearest
  Neighbor Gaussian Process \cr
  https://CRAN.R-project.org/package=spNNGP and randomForest:
  Breiman and Cutler's Random
  Forests for Classification and Regression \cr
  https://CRAN.R-project.org/package=randomForest .

}

\usage{
RFGLS_estimate_spatial(coords, y, X, Xtest = NULL,
                       nrnodes = NULL, nthsize = 20,
                       mtry = 1, pinv_choice = 1,
                       n_omp = 1, ntree = 50, h = 1,
                       sigma.sq = 1, tau.sq = 0.1,
                       phi = 5, nu = 0.5,
                       n.neighbors = 15,
                       cov.model = "exponential",
                       search.type = "tree",
                       param_estimate = FALSE,
                       verbose = FALSE)
}

\arguments{

  \item{coords}{an \eqn{n \times 2}{n x 2} matrix of the observation
  coordinates in \eqn{R^2}{R^2} (e.g., easting and northing). }

  \item{y}{an \eqn{n}{n} length vector of response at the observed coordinates. }

  \item{X}{an \eqn{n \times p}{n x p} matrix of the covariates in the observation coordinates. }

  \item{Xtest}{an \eqn{ntest \times p}{ntest x p} matrix of covariates for prediction locations. Its Structure should be
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

  \item{sigma.sq}{value of sigma square. Default value is 1. }

  \item{tau.sq}{value of tau square. Default value is 0.1. }

  \item{phi}{value of phi. Default value is 5. }

  \item{nu}{value of nu, only required for matern covariance model. Default value is 0.5. }

  \item{n.neighbors}{number of neighbors used in the NNGP. Default value is 15. }

  \item{cov.model}{keyword that specifies the covariance function to be used in modelling the spatial dependence structure
                   among the observations. Supported keywords are: "exponential", "matern", "spherical", and "gaussian"
                   for exponential, matern, spherical and gaussian covariance function respectively. Default value is "exponential". }

  \item{search.type}{keyword that specifies type of nearest neighbor search algorithm to be used. Supported keywords are:
                     "tree" and "brute". Both of them provide the same result, though "tree" should be faster.
                     Default value is "tree". }

  \item{param_estimate}{if \code{TRUE}, using the residuals obtained from fitting a classical RF with default options and \code{nodesize = nthsize},
                        will estimate the coefficeints corresponding to \code{cov.model} from \code{BRISC_estimate} with the deafult options.
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

  Datta, A., S. Banerjee, A.O. Finley, and A.E. Gelfand. (2016)
  Hierarchical Nearest-Neighbor Gaussian process models for large
  geostatistical datasets. Journal of the American Statistical
  Association, 111:800-812.

  Andrew Finley, Abhirup Datta and Sudipto Banerjee (2017). spNNGP: Spatial Regression Models for Large
  Datasets using Nearest Neighbor Gaussian Processes. R package version 0.1.1.
  https://CRAN.R-project.org/package=spNNGP

  Andy Liaw, and Matthew Wiener (2015). randomForest: Breiman and Cutler's Random
  Forests for Classification and Regression.
  R package version 4.6-14. \cr https://CRAN.R-project.org/package=randomForest

}

\author{
  Arkajyoti Saha \email{arkajyotisaha93@gmail.com}, \cr
  Sumanta Basu \email{sumbose@cornell.edu}, \cr
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
n <- 200
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

estimation_result <- RFGLS_estimate_spatial(coords, y, x, ntree = 10)

}

\keyword{model}

