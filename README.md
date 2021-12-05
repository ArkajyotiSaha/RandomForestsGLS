<!-- badges: start -->
[![R-CMD-check](https://github.com/ArkajyotiSaha/RandomForestsGLS/workflows/R-CMD-check/badge.svg)](https://github.com/ArkajyotiSaha/RandomForestsGLS/actions)
<!-- badges: end -->

# RandomForestsGLS
====

## Overview
The R package `RandomForestsGLS: Random Forests for dependent data` fits non-linear regression models on dependent data with Generalized Least Square (GLS) based Random Forest (RF-GLS). Classical Random forests ignore the correlation structure in the data for purpose of greedy partition, mean estimation and resampling for each regression tree. The package implements a RF-GLS proposed in Saha, Basu and Datta (2021) which circumvents the aforementioned problems by incorporating a working correlation structure of the data in for partition, mean estimation and resampling. In this article, it is shown that the greedy split criterion of classical regression trees can be written as an Ordinary Least Square (OLS) optimization with membership in current leaf nodes forming the design matrix. The article extends RF to RF-GLS in a similar fashion to how OLS is extended to GLS by incorporating the covariance structure of the data in the cost function. This ensures that the node splitting and node representatives involve contribution from points belonging to other nodes, weighed by their respective spatial correlations. In classical Random Forest (RF), data points are resampled/subsampled for each of the regression trees, without accounting for their inherent correlation structure. RF-GLS circumvents this problem by resampling/subsampling uncorrelated contrasts instead of original data points. `RandomForestsGLS` implements a fast version of RF-GLS which approximates the working correlation structure using nearest neighbors which makes it suitable for larger datasets.

## Installation
In order to install the development version of the package, please run the following command in R:

```{r }
if (!require("devtools")) install.packages("devtools")
devtools::install_github("ArkajyotiSaha/RandomForestsGLS", ref ="HEAD")
```
For installation of the CRAN version of the package, please use the following command in R:

```{r}
install.packages("RandomForestsGLS")
```

## Example usage: Vignette
The package vignette, available at https://cran.rstudio.com/web/packages/RandomForestsGLS/vignettes/RandomForestsGLS_user_guide.pdf demonstrates with example how the functions available in `RandomForestsGLS` can be used for non-linear regression analysis of dependent data. Specific functions are discussed in much detail in the code documentation of the package. 

## Function documentation

For detailed help on the functions in `RandomForestsGLS` please use the following:
```{r }
?RFGLS_estimate_spatial #(for estimation in spatial data)
?RFGLS_estimate_timeseries #(for estimation in timeseries data)
?RFGLS_predict #(for prediction of mean function)
?RFGLS_predict_spatial #(for prediction of Spatial Response)
```
The function input and outputs are described in detail in the reference manual documentation, available in https://cran.rstudio.com/web/packages/RandomForestsGLS/RandomForestsGLS.pdf .

## Community guidelines

Please report issues, bugs or problem with the software at https://github.com/ArkajyotiSaha/RandomForestsGLS/issues . For contribution to the software and support please get in touch with the maintainer Arkajyoti Saha (arkajyotisaha93@gmail.com).

## Note
Some code blocks are borrowed from the R packages: `spNNGP: Spatial Regression Models for Large Datasets using Nearest Neighbor Gaussian Processes` https://CRAN.R-project.org/package=spNNGP and `randomForest: Breiman and Cutler's Random Forests for Classification and Regression` https://CRAN.R-project.org/package=randomForest .
RF-GLS uses nearest neighbor Gaussian process (NNGP) to approximate the covariance structure in the data, tools necessary to implement the NNGP are borrowed from the spNNGP package, which include `util.cpp` and parts of `updateBF_org` and `RFGLS_BFcpp` in `RFGLS.cpp`. The basic building blocks for Random Forest is borrowed from `randomForest` which include parts of `RFGLStree_cpp`, `findBestSplit` `RFGLSpredicttree_cpp` in `RFGLS.cpp`.


## Citation
Please cite the following paper when you use RF-GLS

Saha, Arkajyoti, Sumanta Basu, and Abhirup Datta. "Random forests for spatially dependent data." Journal of the American Statistical Association (2021): 1-19. https://www.tandfonline.com/doi/full/10.1080/01621459.2021.1950003 .

