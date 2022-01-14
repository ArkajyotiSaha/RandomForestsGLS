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
date: 14 January 2022
bibliography: paper.bib
---

# Summary
With the modern advancements in geographical information systems, remote sensing technologies, low-cost sensors, we are increasingly encountering massive datasets indexed by geo-locations and time-stamps. Modeling such high throughput spatio-temporal data needs to carefully account for spatial/serial dependence in the data, i.e., model the structures (patterns) of the data along space or time. A general model for describing dependent observations $(y_1, y_2, \cdots, y_n)$ and covariates $(\mathbf x_1,\ldots,\mathbf x_n)$ can be formulated as $y_i = f(\mathbf x_i) + \epsilon_i$, where $f(\mathbf x_i)$ denotes the mean component and $\epsilon_i$ is the residual error, which accounts for the dependency in the data, beyond what is explained by the covariates. In this article, we work with second order dependency, i.e. we assume that the dependence is captured through the covariance function of the correlated stochastic process $\epsilon_i$. Since the data are indexed by locations (spatial data) or time-points (temporal data), the correlation is a function of the "spatial distance" or "time-lag" between two observations.

The primary purpose of the package `RandomForestsGLS` is to fit nonlinear regression for spatial and temporal data. This package performs estimation through a novel rendition of Random Forests (RF); namely, RF-GLS, which makes use of the dependence structure in the data. With increasing computing capacity of personal computers, non-linear Machine Learning (ML) methods have become increasingly commonplace for regression and classification tasks due to their ability of capture complex interactions among the variables, that cannot be modeled by linear regression. However, many modern ML software lack the capability to efficiently account for the dependence structure in the data which leads to sub-optimal estimation. On the other hand, specialized software for spatial/temporal data are capable of properly modeling data correlation, but usually assume the relationships between the response and the covariates.  `RandomForestsGLS` brings the best of both worlds together by bridging the gap between these two approaches by explicitly modeling the spatial/serial data correlation in the RF fitting procedure to substantially improve estimation of the regression function. In ML workflow, `RandomForestsGLS` can substitute existing ML methods in model training to take care of the dependence structure. Similarly, for traditional spatial modeling  frameworks, `RandomForestsGLS` can be used instead of linear model based methods to account for nonlinearity. Additionally, `RandomForestsGLS` seamlessly leverages kriging to perform predictions at new locations for geo-spatial data, a primary objective in many spatial analysis.

# Statement of need
RF is an ensemble of regression trees built on subsamples of the data. Regression trees iteratively partition the covariate space by greedily optimizing a split criterion, that minimizes the sum of intra-node variances. This split criterion can also be expressed as an Ordinary Least-Square (OLS) loss with design matrix corresponding to the leaf nodes. OLS loss used in split criterion and node representative update rule in RF does not account for the dependence in the data, hence is not optimal. Additionally, responses are resampled in the "bagging" (bootstrap aggregation) [@breiman1996bagging] step in RF, which, under dependence setting violates the primary requirement of bootstrapping that the data are independent and identically distributed. RF-GLS, implemented in `RandomForestsGLS`, mitigates the issue by using a dependency adjusted split criterion which incorporates the covariance structure in determining the optimal partition as well as the leaf node representatives. Additionally, RF-GLS facilitates resampling decorrelated contrasts, which is suitable for bootstrap sampling [@saha2021random]. The dependency is modeled using GP in a mixed model framework, allowing both estimation and prediction.

`RandomForestsGLS` is a statistical machine learning oriented [R](https://cran.r-project.org) package for fitting RF-GLS on dependent data. `R` being a free software available under GNU General Public License with tailor-made statistical operators and vast library of statistical software, is very popular in the statistics community. The RF-GLS algorithm described in [@saha2021random] involves computationally intensive linear algebra operations in nested loops which are especially slow in an interpreted language like `R`. In order to optimize the implementation of the algorithm, it becomes a necessity to use low-level languages which have a significantly steeper learning curve than that of `R` and lacks a majority of the inbuilt statistical operators in `R`.\
`RandomForestsGLS` brings the best of both worlds together by optimizing the implementation of the computationally challenging RF-GLS algorithm in the background with an easy to use `R` user interface.

`RandomForestsGLS` focuses on fast, parallelizable implementations of RF-GLS for spatial and time series data, which includes popular choices for covariance functions for both spatial (Matérn GP) and time series (autoregressive) model. The package is primarily designed to be used by researchers associated with the fields of statistical machine learning, spatial statistics, time series analysis, and scientific applications of them. A significant part of the code has already been used in @saha2021random. With the combination of speed and ease-of-use that `RandomForestsGLS` brings to the table regarding non-linear regression analysis in dependent data, we hope to see this package being used in a plethora of future scientific and methodological explorations.


 
# State of the field

There are a number of `R` packages that provides the option for implementing classical RF. Most notable of them being [randomForest](https://CRAN.R-project.org/package=randomForest), which implements Breiman's Random Forests [@breiman2001random] for Classification and
Regression using the [Fortran](https://fortran-lang.org/) original code by Leo Breiman and Adele Cutler. Some of the other packages are [xgboost](https://CRAN.R-project.org/package=xgboost), [randomForestSRC](https://CRAN.R-project.org/package=randomForestSRC), [ranger](https://CRAN.R-project.org/package=ranger), [Rborist](https://CRAN.R-project.org/package=Rborist). For a detailed overview of the we refer the reader to [CRAN Task View: Machine Learning & Statistical Learning](https://cran.r-project.org/web/views/MachineLearning.html). To the best of our knowledge, none of these packages explicitly account for model correlation, common in spatial and time-series data.

While the classical RF has also been used for a number of geo-spatial applications (see @saha2021random for references),
most of them do not make methodological adjustments in RF to account for the
spatial correlation in the data ([CRAN Task View: Analysis of Spatial Data](https://cran.r-project.org/web/views/Spatial.html)). Two recent works attempt to explicitly use spatial information in RF. @hengl2018random adds all pairwise spatial-distances as additional covariates. This does not need a dedicated package as it uses the classical RF. A detailed tutorial to implement this in `R` is available in [GeoMLA](https://github.com/thengl/GeoMLA). The other approach, @georganos2019geographical proposes geographically local estimation of RF. Though this approach also employs the classical RF implemented in `randomForest`, it is implemented in a standalone package [SpatialML](https://CRAN.R-project.org/package=SpatialML). Both approaches abandon the spatial mixed model framework and try to account for the dependence structure in the data by incorporating additional spatial covariates. A downside to this is unnecessary escalation of the problem to high-dimensional settings, being unable to parsimoniously encode structured spatial dependence via the GP covariance. This adversely affects the prediction performance of these methods compared to that of RF-GLS when there is a dominant covariate effect.

Additionally, these only focus on prediction and are unable to estimate the covariate effect (mean function) separately from the spatial effect, which can be of independent interest to the researcher. As far as time series analysis with RF is concerned, to the best of our knowledge, no dedicated `R` package is available at the moment ([CRAN Task View: Time Series Analysis](https://cran.r-project.org/web/views/TimeSeries.html)). The standard practice involves incorporating prior responses for a desired number of lags as covariates [@BASAK2019552] and using block bootstrap for "bagging" to retain the local covariance structure. This data is then used with the classical RF for the purpose of forecasting. These approaches also suffer from similar problems as discussed before that are explicitly addressed in RF-GLS.

# The RandomForestsGLS package

We provide a brief overview of the user functionality of the package. The package [vignette](https://cran.r-project.org/web/packages/RandomForestsGLS/vignettes/RandomForestsGLS_user_guide.pdf) delves deeper into this and demonstrates with example how the functions available in the package can be used for non-linear regression analysis of dependent data. Specific functions are discussed in detail in the code documentation of the package.

## RF to RF-GLS: Accounting for correlation structure

In classical RF, which is an average of many regression trees, each node in a regression tree is split by optimizing the CART split criterion in @breiman1984classification. It can be rewritten in the following way. 
$$
v_{n}^{CART}((d,c)) =  \frac{1}{n} \left( \|\mathbf{Y} - \mathbf Z^{(0)}\boldsymbol{\hat{\beta}}(\mathbf Z^{(0)})\|_2^2 - \|\mathbf Y - \mathbf Z \boldsymbol{\hat{\beta}}(\mathbf Z)\|_2^2 \right).
$$
where, $\mathbf Z^{(0)}$ and $\mathbf Z$ are the membersip matrices for the leaf nodes of the tree before and after the potential node split. $(d,c)$ denotes a potential cut (location of the split), with $d$ and $c$ being the cut direction (choice of the covariate) and cutoff point (value of the covariate) respectively, $\boldsymbol{\hat{\beta}} (\mathbf Z)$ are the leaf node representatives given by OLS estimates corresponding to design matrix $\mathbf Z$ and can be written as: 

$$\boldsymbol{\hat{\beta}} (\mathbf Z) = \left(\mathbf Z ^\top \mathbf Z \right)^{-1} \mathbf Z ^\top \mathbf y$$

We observe that the split criterion is the difference of OLS loss functions before and after the cut with the design matrix of membership of the leaf nodes. We can incorporate the correlation structure of the data in the split criterion by replacing the OLS loss with GLS loss. The modified split criterion can be rewritten as:

$$
\begin{aligned}
v_{n,\mathbf Q}^{DART}((d,c)) = 
&\frac{1}{n} \Bigg[\left(\mathbf{Y} - \mathbf{Z}^{(0)}\boldsymbol{\hat{\beta}}_{GLS}(\mathbf Z^{(0)}) \right)^\top \mathbf Q\left(\mathbf{Y} - \mathbf{Z}^{(0)}\boldsymbol{\hat{\beta}}_{GLS}(\mathbf Z^{(0)}) \right)\\ &-\left(\mathbf{Y} - \mathbf{Z}\boldsymbol{\hat{\beta}}_{GLS}(\mathbf Z) \right)^\top \mathbf Q\left(\mathbf{Y} - \mathbf{Z}\boldsymbol{\hat{\beta}}_{GLS}(\mathbf Z) \right) \Bigg].
\end{aligned}
$$

where, $\mathbf Q$ is the inverse of the working covariance matrix that models the spatial/serial dependence and $\boldsymbol{\hat{\beta}}_{GLS} (\mathbf Z)$ are the leaf node representatives given by the GLS estimates corresponding to design matrix $\mathbf Z$ and can be written as follows:

$$\boldsymbol{\hat{\beta}}_{GLS} (\mathbf Z) = \left(\mathbf Z ^\top \mathbf Q \mathbf Z \right)^{-1} \mathbf Z ^\top \mathbf Q \mathbf y.$$

## Spatial Data
### Model
We consider spatial point referenced data with the following mixed model:

$$
y_i = m(\mathbf{x}_i) + w(\mathbf{s}_i) + \epsilon_i;
$$
where, $y_i, \mathbf{x}_i$ respectively denotes the observed response and the covariate corresponding to the $i^{th}$ observed location $\mathbf{s}_i$. $m(\mathbf{x}_i)$ denotes the covariate effect; spatial random effect, $w (\mathbf{s})$ accounts for spatial dependence beyond covariates modeled using a GP, and $\mathbf{\epsilon}$ accounts for the independent and identically distributed random Gaussian noise. 



### Fitting
In the spatial mixture model setting, the package `RandomForestsGLS` allows for fitting $m(.)$ using RF-GLS using `RFGLS_estimate_spatial`. Spatial random effects are modeled using GP as is the practice. For model fitting, we use the computationally convenient Nearest Neighbor Gaussian Process (NNGP) [@nngp]. If the covariance parameters are unknown, they are automatically estimated from the specified covariance model and used in model fitting.

With coordinates `coords`, response `y` and the covariates `x`, we can fit RF-GLS for spatial data as:

```
est <- RFGLS_estimate_spatial(coords, y, x)
```


### Prediction
Given a fitted model and the covariates `RFGLS_predict` predicts the covariate effects given a new covariate value. We also offer spatial prediction at new locations with `RFGLS_predict_spatial` which performs a non-linear kriging by combining the non-linear mean estimate from the covariate effects in `RFGLS_predict` and spatial kriging estimate from the [BRISC](https://CRAN.R-project.org/package=BRISC) package.

With new covariates `Xtest` and fitted RF-GLS model `est`, the covariate effect can be predicted using

```
RFGLS_predict <- RFGLS_predict(est, Xtest)
```

The spatial predictions at new locations `Coordstest`, can be obtained through 

```
RFGLS_predict_spatial <- RFGLS_predict_spatial(est, Xtest, Coordstest)
```

## Autoregressive (AR) Time Series Data
### Model

RF-GLS can also be used for function estimation in a time series setting under autoregressive (AR) errors. We consider time series data with errors from an AR(q) (autoregressive process of lag q) process as follows:

$$
y_t = m(\mathbf{x}_t) + e_t;\: e_t = \sum_{i = 1}^q\rho_i e_{t-i} + \eta_t
$$

where, $y_i, \mathbf{x}_i$ denotes the response and the covariate corresponding to the $t^{th}$ time point, $e_t$ is an AR(q) process, $\eta_t$ denotes the i.i.d. white noise and $(\rho_1, \cdots, \rho_q)$ are the model parameters that captures the dependence of $e_t$ on $(e_{t - 1}, \cdots, e_{t-q})$.

### Fitting
In the AR time series scenario, the package `RandomForestsGLS` allows for fitting $m(.)$ using RF-GLS with `RFGLS_estimate_timeseries`. RF-GLS exploits the sparsity of the closed form precision matrix of the AR process for model fitting and prediction of mean function $m(.)$. If the AR model parameters (coefficients and order of autoregressive process) are unknown the code automatically estimates them from AR models with specified lags.

With response `y` and the covariates `x`, we can fit RF-GLS for temporal data as:

```
est <- RFGLS_estimate_timeseries(y, x)
```

### Prediction of covariate effects
Prediction of covariate effects in AR process models are similar to that of spatial data and are performed with `RFGLS_predict` as in the case of spatial data.


# Discussion

In this package, we have developed an efficient, parallel implementation of the RF-GLS method proposed in @saha2021random. This accounts for the correlation in the data by incorporating it in the node splitting criteria and the node representative update rule. The package accounts for spatial correlation, modeled with Matérn GP and serial autocorrelation. More often than not, the model parameters are unknown to the users, hence the package has inbuilt parameter estimation corresponding to both the covariance structures. Efficient implementation through C/C++ takes advantage of the NNGP approximation and scalable node splitting update rules which reduces the execution time. For details regarding package features, unknown parameter estimation and parallelization, we refer the readers to the package [README](https://github.com/ArkajyotiSaha/RandomForestsGLS#readme). Present implementation of RF-GLS is slower than that of RF due to the additional computational complexity associated with the split criteria evaluation. In RF-GLS, evaluation of cost function, corresponding to a potential split has $O(t^3)$ computational complexity, where $t$ is the number of leaf nodes at the present split. For very deep trees, this would lead to significant added computational burden over RF, where this step has $O(1)$ complexity. As a future extension of the package, implementation of covariate specific parallel optimization corresponding to each node splitting and improving overall computational complexity of the algorithm can be of independent research interest. We also plan to implement additional forms of time series dependency and perform thorough empirical validation, prior to making it available in the software.

# Acknowledgements

AS and AD were supported by NSF award DMS-1915803. SB was supported by an NSF award DMS-1812128, and an NIH award R01GM135926.

# References

