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
    affiliation: 1
affiliations:
 - name: Department of Biostatistics, Johns Hopkins Bloomberg School of Public Health
   index: 1
 - name: Department of Statistics and Data Science, Cornell University
   index: 2
date: 29 December 2020
bibliography: paper.bib
---

# Summary

Linear mixed-models, consisting of a linear covariate effect and a Gaussian Process (GP) distributed random effect, are widely used for analyses of dependant data. We consider the setting where the covariate effect is non-linear. Random forests (RF) are widely popular for estimating such non-linear regression functions, however, the classical RF ignores the dependence structure in a mixed-model setting. A novel and well-principled rendition of RF under dependent data settings was developed in `@saha2020random`. The `RandomForestsGLS` fits non-linear regression models on dependent data with the Generalised Least Square (GLS) based Random Forest (RF-GLS) that was proposed in `@saha2020random`. Our package offers fitting a parallelizable and scalable RF that accounts for spatial and temporal correlation in the data. We predict covariate effects and also offer spatial response prediction at new locations.

# Statement of need

The Ordinary Least-Square (OLS) loss used in splitting of RF is suboptimal in dependant settings. Additionally, "bagging" (bootstrap aggregation) `[@breiman1996bagging]` in RF resamples data, hence under correlated data, violates the assumption of independent and identically distributed (i.i.d.) data unit fundamental to bootstrapping. RF-GLS implemented in `RandomForestsGLS` solves the problem by using a dependency adjusted split criterion which incorporates the covariance structure in determining the optimal partition as well as the leaf node representatives. Additionally, RF-GLS felicitates resampling decorrelated contrasts, which is suitable for bootstrap sampling `[@saha2020random]`. 

`RandomForestsGLS` is a statistical machine learning oriented [R](https://cran.r-project.org) package for fitting RF-GLS on dependent data. `R` being a free software available under GNU General Public License with tailor-made statistical operators and vast library of statistical softwares, is very popular in the statistics community. The RF-GLS algorithm described in `[@saha2020random]` involves computationally intensive linear algebra operations in nested loops which are especially slow in an interpretaed language like `R`. In order to optimize the implementation of the algorithm, it becomes a necessity to use low-level languages which have a significantly steeper learning curve than that of `R` and lacks a majority of the inbuilt statistical operators in `R`. `RandomForestsGLS` brings the best of both worlds togeher by optimizing the implementation the computationally challenging RF-GLS algorithm in the background with an easy to use `R` user-interface.

`RandomForestsGLS` focuses on fast, parallelizable implementations of RF-GLS for spatial and time series data which includes popular choices for covariance functions for both spatial (Matern) and time series (autoregressive) model. The package is primarily designed to be used by researchers associated with the fields of statistical machine learning, spatial statistics, time series analysis and scientific applications of them. A significant part of the code has already being used in `@saha2020random`. With the combination of speed and ease-of-use that `RandomForestsGLS` brings to the table regarding non-linear regression analysis in dependant data, we hope to see this package used in plethora of future scientific and methodological explorations.


 
# State of the field

There are a number of `R` packages that provides the option for implementing classical RF. Most popular of them being [randomForest](https://CRAN.R-project.org/package=randomForest), which implements Breiman Random Forests `[@breiman2001random]` for Classification and
Regression using the [Fortran](https://fortran-lang.org/) original code by Leo Breiman and Adele Cutler. Some of the other notable members are [xgboost](https://CRAN.R-project.org/package=xgboost), [randomForestSRC](https://CRAN.R-project.org/package=randomForestSRC), [ranger](https://CRAN.R-project.org/package=ranger), [Rborist](https://CRAN.R-project.org/package=Rborist). For a detailed overview of the we refer the reader to [CRAN Task View: Machine Learning & Statistical Learning](https://cran.r-project.org/web/views/MachineLearning.html). To the best of our knowledge, none of these packages explicitly account for correlated data.

While RF has also been used for a number of geo-spatial applications (see `@saha2020random` for details),
most of the aforementioned work do not make procedural considerations for the RF algorithm to address the
spatial correlation in the data. Two recent works attempt to explicitly use spatial information in RF. `@hengl2018random` adds all pairwise spatial-distances as additional covariates. This doesn't have a dedicated package as it uses the classical RF. A detailed tutorial to implement this in `R` is avialble in [GeoMLA](https://github.com/thengl/GeoMLA). The other approach, `@georganos2019geographical` proposes geographically local estimation of RF. Thoug this approach also employs the classical RF implemented in `randomForest`, it is implemented in a standalone package [SpatialML](https://CRAN.R-project.org/package=SpatialML). Both approaches abandon the spatial mixed model framework. A downside to this is unnecessary escalation of the problem to high-dimensional settings, being unable to parsimoniously encode structured spatial dependence via the GP covariance. Additionally these only focus on prediction and are unable to estimate the covariate effect (mean function) which can be of independent interest to the researcher. As far as time eries analysis with RF is concerned, to the best of our knowledge no dedicated `R` package is avialble at the moment [CRAN Task View: Time Series Analysis](https://cran.r-project.org/web/views/TimeSeries.html). The standard practice involves incorporating prior responses at desired number of lags as covariates `[@BASAK2019552]` and using block bootstrap for "bagging" to retain the local covariance structure. This data is then used with classical RF for the purpose of forecasting. These approachs also suffer from similar problems as discussed before that are explicitly addressed in RF-GLS.

# The RandomForestsGLS package

We provide a brief overview of the user functionality of the package. The package vignette delves deeper into this and demonstrates with example how the functions available in `RandomForestsGLS` can be used for non-linear regression analysis of dependent data. Specific functions are discussed in much detail in the code documentation of the package.

## Spatial Data
### Model
We consider spatial point referenced data with the follwing model:

$$
y_i = m(\mathbf{x}_i) + w(\mathbf{s}_i) + \epsilon_i;
$$
where, $y_i, \mathbf{x}_i$ respectively denotes the observed response and the covariate corresponding to the $i^{th}$ observed location $\mathbf{s}_i$. $m(\mathbf{x}_i)$ denotes the covariate effect, spatial random effect, $w (\mathbf{s})$ accounts for spatial dependence beyond covariates, and $\mathbf{\epsilon}$ accounts for the independent and identically distributed random Gaussian noise. 

### Fitting
In the spatial mixture model setting, the package `RandomForestsGLS` allows for fiting $m(.)$ using RF-GLS using `RFGLS_estimate_spatial`. Spatial random effects are modeled using Gaussian Process as is the practice. For model fitting, we use the computationally convenient Nearest Neighbor Gaussian Process (NNGP) `@nngp`. If the covariance parameters are unknown the code automatically estimates them from the specified covariance model and uses them in model fitting.

### Prediction
Given a fitted model and the covariates `RFGLS_predict` predicts the covariate effects. We also offer spatial prediction at new locations with `RFGLS_predict_spatial` which performs a non-linear krigging by combining the non-linear mean estimate from the covariate effects in `RFGLS_predict` and spatial krigging estimate from `BRISC`.

## Autoregressive (AR) Time Series Data
### Model

RF-GLS can also be used for function estimation in a time series setting under autoregressive errors. We consider time series data with errors from a AR(q) process as follows:

$$
y_t = m(\mathbf{x}_t) + e_t;\: e_t = \sum_{i = 1}^q\rho_i e_{t-i} + \eta_t
$$

where, $y_i, \mathbf{x}_i$ denotes the response and the covariate corresponding to the $t^{th}$ time point, $e_t$ is an AR(q) pprocess, $\eta_t$ denotes the i.i.d. white noise and $(\rho_1, \cdots, \rho_q)$ are the model parameters that captures the dependence of $e_t$ on $(e_{t - 1}, \cdots, e_{t-q})$.

### Fitting
In the AR time series scenario, the package `RandomForestsGLS` allows for fiting $m(.)$ using RF-GLS with `RFGLS_estimate_timeseries`. RF-GLS exploits the sparsity of the closed form precision matrix of the AR process for model fitting and prediction of mean function $m(.)$. If the AR model parameters are unknown the code automatically estimates them from AR models with specificed lags.

### Prediction of covariate effects
Prediction of covariate effects in AR process models are similar to that of spatial data and are performed with `RFGLS_predict`.

### Parallelization

For `RFGLS_estimate_spatial`, `RFGLS_estimate_timeseries`, `RFGLS_predict` and `RFGLS_predict_spatial` one can also take the advantage of parallelization, contingent upon the availability of multiple cores. With very small dataset and small number of trees, communication overhead between the nodes for parallelization outweighs the benefits of the parallel computing hence it is recommended to parallelize only for moderately large dataset and/or number of trees.



# Package Features
The source code of the package are written in [C](https://en.cppreference.com/w/c/language)/[C++](https://isocpp.org/) for sake of optimization. The functions available to the user are wrappers around the source code, built with `R`'s foreign language interface. For the basic structure of the code we make use of the source code of the regression trees in `R` based implementation of the classical RF in `randomForest`. As the split criterion in RF-GLS involves computationally intensive linear algeba operation in nested loops, we use `Fortran`'s Basic Linear Algebra Subprograms ([BLAS](http://www.netlib.org/blas/)) and Linear Algebra Package ([LAPACK](http://www.netlib.org/lapack/)). This is achieved by storing all matrices in contiguous memory column-major format. We also offer multicore computation by building each regression tree independly.

As the split in RF-GLS requires optimizing cost function involving the cholsky of the precision matrix, use of the full dense precision matrix in spatial processes becomes taxing on typical personal computing resources both in terms of computational cost ($O(n^3)$) and storage cost ($O(n^2)$). Inorder to circumvent this problem, we use NNGP to replace the dense graph among spatial locations with a nearest neighbor graphical model. This was shown to directly yield a sparse Cholesky factor that offers an excellent approximation to its original dense counterpart. 
We implement a convenient nearest neighbor search following [spNNGP](https://CRAN.R-project.org/package=spNNGP) `[@spnngppaper]` and efficient sparse matrix multiplication of the associated form in [BRISC](https://CRAN.R-project.org/package=BRISC)`[@brisc]`. The structure of the loops used in the process fecilitaes parallelization using `openMP` for this stage of calculation. In time series analysis, the sparsity in the precision matrix is inherently induced by AR covariance structure. 

# References
