---
title: 'RandomForestsGLS: An R package for Random Forests for dependent data'
tags:
  - R
  - spatial statistics
  - Gaussian Processes
  - Random forests
  - generalized least-squares
authors:
  - name: Arkajyoti Saha
    affiliation: 1
  - name: Sumanta Basu
    affiliation: 2
  - name: Abhirup Datta
    affiliation: 3
affiliations:
 - name: Departments of Statistics, University of Washington
   index: 1
 - name: Department of Statistics and Data Science, Cornell University
   index: 2
 - name: Department of Biostatistics, Johns Hopkins Bloomberg School of Public Health
   index: 3
date: 23 February 2022
bibliography: paper.bib
---

# Summary
With the modern advances in geographical information systems, remote sensing technologies, and low-cost sensors, we are increasingly encountering datasets where we need to account for spatial or serial dependence. Dependent observations $(y_1, y_2, \cdots, y_n)$ with covariates $(\mathbf x_1,\ldots,\mathbf x_n)$ can be modeled non-parametrically as $y_i = m(\mathbf x_i) + \epsilon_i$, where $m(\mathbf x_i)$ is mean component and $\epsilon_i$ accounts for the dependency in data. We assume that dependence is captured through a covariance function of the correlated stochastic process $\epsilon_i$ (second order dependence). The correlation is typically a function of "spatial distance" or "time-lag" between two observations. 

Unlike linear regression, non-linear Machine Learning (ML) methods for estimating the regression function $m$ can capture complex interactions among the variables. However, they often fail to account for the dependence structure, resulting in sub-optimal estimation. On the other hand, specialized software for spatial/temporal data properly models data correlation but lacks flexibility in modeling the mean function $m$ by only focusing on linear models. `RandomForestsGLS` bridges the gap through a novel rendition of Random Forests (RF) -- namely, RF-GLS -- by explicitly modeling the spatial/serial data correlation in the RF fitting procedure to substantially improve the estimation of the mean function. Additionally, `RandomForestsGLS` leverages kriging to perform predictions at new locations for geo-spatial data.

# Statement of need
`RandomForestsGLS` is a statistical machine learning oriented [R](https://cran.r-project.org) package for fitting RF-GLS on dependent data. The RF-GLS algorithm described in [@saha2021random] involves computationally intensive linear algebra operations in nested loops which are especially slow in an interpreted language like `R`. `RandomForestsGLS` efficiently implements RF-GLS algorithm in a low-level language with a user-friendly interface in `R`, a popular computational software (free under GNU General Public License) in the statistics community. `RandomForestsGLS` focuses on fast, parallelizable implementations of RF-GLS for spatial and time series data, which includes popular choices for covariance functions for both spatial (Matérn GP) and time series (autoregressive) data. The package is primarily designed to be used by researchers associated with the fields of statistical machine learning, spatial statistics, time series analysis, and their scientific applications. A significant part of the code has already been used in @saha2021random. With the combination of speed and ease-of-use that `RandomForestsGLS` brings to the table regarding non-linear regression analysis for dependent data, we hope to see this package being used in a plethora of future scientific and methodological explorations.

# State of the field

Several `R` packages implement classical RF. Most notable of them being [randomForest](https://CRAN.R-project.org/package=randomForest), which implements Breiman's Random Forests [@breiman2001random] for Classification and
Regression using the [Fortran](https://fortran-lang.org/). Some of the other packages are [xgboost](https://CRAN.R-project.org/package=xgboost), [randomForestSRC](https://CRAN.R-project.org/package=randomForestSRC), [ranger](https://CRAN.R-project.org/package=ranger), [Rborist](https://CRAN.R-project.org/package=Rborist). For a detailed overview of it, we refer the reader to [CRAN Task View: Machine Learning & Statistical Learning](https://cran.r-project.org/web/views/MachineLearning.html). To the best of our knowledge, none of these packages explicitly account for spatial and/or temporal correlation.

Classical RF has been used in geo-spatial and temporal applications (see @saha2021random for references) without making methodological adjustments to account for spatial dependencies. ([CRAN Task View: Analysis of Spatial Data](https://cran.r-project.org/web/views/Spatial.html), [CRAN Task View: Time Series Analysis](https://cran.r-project.org/web/views/TimeSeries.html)). Two recent works that attempt to explicitly use spatial information in RF for prediction purposes, are @hengl2018random ( [GeoMLA](https://github.com/thengl/GeoMLA)) and @georganos2019geographical ( [SpatialML](https://CRAN.R-project.org/package=SpatialML)) (see @saha2021random for details). Both approaches try to account for the dependence structure by incorporating additional spatial covariates, which adversely affect the prediction performance in the presence of a dominant covariate effect (@saha2021random). Additionally, unlike RF-GLS, they cannot estimate covariate effects separately from the spatial effect, which can be of independent interest. The RF for temporal data, proposed in @BASAK2019552, suffers from similar shortcomings.

# The RandomForestsGLS package

We provide a brief overview of the functionality of the package. The package [vignette](https://cran.r-project.org/web/packages/RandomForestsGLS/vignettes/RandomForestsGLS_user_guide.pdf) demonstrates with example how to use the package for non-linear regression analysis of dependent data. Specific functions are discussed in detail in the documentation of the package.

## RF to RF-GLS: Accounting for correlation structure

In classical RF, which is an average of many regression trees, each node in a regression tree is split by optimizing the CART split criterion in @breiman1984classification. It can be rewritten in the following way: 
$$
v_{n}^{CART}((d,c)) =  \frac{1}{n} \left( \|\mathbf{Y} - \mathbf Z^{(0)}\boldsymbol{\hat{\beta}}(\mathbf Z^{(0)})\|_2^2 - \|\mathbf Y - \mathbf Z \boldsymbol{\hat{\beta}}(\mathbf Z)\|_2^2 \right).
$$

where $\mathbf Z^{(0)}$ and $\mathbf Z$ are the binary membership matrices for the leaf nodes of the tree before and after the potential node split. $(d,c)$ denotes a potential cut (location of the split), with $d$ and $c$ being the cut direction (choice of the covariate) and cutoff point (value of the covariate) respectively, $\boldsymbol{\hat{\beta}} (\mathbf Z)$ are the leaf node representatives given by OLS estimates corresponding to a design matrix $\mathbf Z$ and can be written as: 

$$\boldsymbol{\hat{\beta}} (\mathbf Z) = \left(\mathbf Z ^\top \mathbf Z \right)^{-1} \mathbf Z ^\top \mathbf y$$

We observe that the split criterion is the difference of OLS loss functions before and after the cut with the respective design matrices of membership of the leaf nodes. We can incorporate the correlation structure of the data in the split criterion by replacing the OLS loss with GLS loss as is traditionally done in linear models. The modified split criterion can be rewritten as:

$$
\begin{aligned}
v_{n,\mathbf Q}^{DART}((d,c)) = 
&\frac{1}{n} \Bigg[\left(\mathbf{Y} - \mathbf{Z}^{(0)}\boldsymbol{\hat{\beta}}_{GLS}(\mathbf Z^{(0)}) \right)^\top \mathbf Q\left(\mathbf{Y} - \mathbf{Z}^{(0)}\boldsymbol{\hat{\beta}}_{GLS}(\mathbf Z^{(0)}) \right)\\ &-\left(\mathbf{Y} - \mathbf{Z}\boldsymbol{\hat{\beta}}_{GLS}(\mathbf Z) \right)^\top \mathbf Q\left(\mathbf{Y} - \mathbf{Z}\boldsymbol{\hat{\beta}}_{GLS}(\mathbf Z) \right) \Bigg].
\end{aligned}
$$

where $\mathbf Q$ is the inverse of the working covariance matrix that models the spatial/serial dependence and $\boldsymbol{\hat{\beta}}_{GLS} (\mathbf Z)$ are the leaf node representatives given by the GLS estimates corresponding to a design matrix $\mathbf Z$ and can be written as follows:

$$\boldsymbol{\hat{\beta}}_{GLS} (\mathbf Z) = \left(\mathbf Z ^\top \mathbf Q \mathbf Z \right)^{-1} \mathbf Z ^\top \mathbf Q \mathbf y.$$

## Spatial Data
### Model
We consider spatial point-referenced data with the following mixed model:

$$
y_i = m(\mathbf{x}_i) + w(\mathbf{s}_i) + \epsilon_i;
$$

where $y_i, \mathbf{x}_i$ respectively denotes the observed response and the covariate corresponding to the $i^{th}$ observed location $\mathbf{s}_i$. $m(\mathbf{x}_i)$ denotes the covariate effect; spatial random effect, $w (\mathbf{s})$ accounts for spatial dependence beyond covariates modeled using a GP, and $\mathbf{\epsilon}$ accounts for independent and identically distributed random Gaussian noise. 



### Fitting & Prediction
Spatial random effects are modeled using a Gaussian process (GP) as is the practice. We use the computationally convenient Nearest Neighbor GP (NNGP) [@nngp]. Model parameters, if unknown, are estimated from the data (@saha2021random). Alongside predicting covariate effects for a new covariate value, we also offer spatial prediction at new locations with non-linear kriging by combining the non-linear mean estimate and spatial kriging estimate from the [BRISC](https://CRAN.R-project.org/package=BRISC) [@brisc] package.


## Autoregressive (AR) Time Series Data
### Model

RF-GLS can also be used for function estimation in a time series setting under autoregressive (AR) errors. We consider time series data with errors from an AR($q$) (autoregressive process of lag $q$) process as follows:

$$
y_t = m(\mathbf{x}_t) + e_t;\: e_t = \sum_{i = 1}^q\rho_i e_{t-i} + \eta_t
$$

where $y_i, \mathbf{x}_i$ denote the response and the covariate corresponding to the $t^{th}$ time point, respectively, $e_t$ is an AR(q) process, $\eta_t$ denotes i.i.d. white noise and $(\rho_1, \cdots, \rho_q)$ are the model parameters that capture the dependence of $e_t$ on $(e_{t - 1}, \cdots, e_{t-q})$.

### Fitting & Prediction
RF-GLS exploits the sparsity of the closed form precision matrix of the AR process for model fitting and prediction of the mean function $m(.)$. If the AR model parameters (coefficients and order of autoregressive process) are unknown, the code automatically estimates them from AR models with specified lags. The prediction of covariate effect for time series data is like that of spatial data. 

# Discussion

This package provides an efficient, parallel implementation of the RF-GLS method proposed in @saha2021random which accounts for correlated data with modified node splitting criteria and node representative update rule. The package accounts for spatial correlation via Matérn GP or serial autocorrelation. It provides parameter estimation for both covariance structures. Efficient implementation through C/C++ takes advantage of the NNGP approximation and scalable node splitting update rules which reduces the execution time. More details on package features, parameter estimation, and parallelization can be found in the package [README](https://github.com/ArkajyotiSaha/RandomForestsGLS#readme). Since the computational complexity of evaluating potential splits is cubic in the number of leaf nodes for RF-GLS, but constant for standard RF, improving the computational complexity of the algorithm is of independent research interest. We also plan to implement and validate additional forms of time series dependency.

# Acknowledgements

AS and AD were supported by NSF award DMS-1915803. AD was supported by NIEHS award R01ES033739. SB was supported by an NSF award DMS-1812128, and an NIH award R01GM135926.

# References

