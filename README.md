
# RandomForestsGLS
====

## Overview
The R package `RandomForestsGLS: Random Forests for dependent data` fits non-linear regression models on dependent data with Generalised Least Square (GLS) based Random Forest (RF-GLS). Classical Random forests ignore the correlation structure in the data for purpose of greedy partition, mean estimation and resampling for each regression tree. The package implements a RF-GLS which circumvents the aforementioned problems by incorporating a working correlation structure of the data in for partition, mean estimation and resampling.


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

## Function description & documentation

For detailed help on the functions in `RandomForestsGLS` please use the following:
```{r }
?RFGLS_estimate_spatial #(for estimation in spatial data)
?RFGLS_estimate_timeseries #(for estimation in timeseries data)
?RFGLS_predict #(for prediction of mean functtion)
?RFGLS_predict_spatial #(for prediction of Spatial Response)
```
The function input and outputs are described in detail in the reference manual documentation, available in https://cran.rstudio.com/web/packages/RandomForestsGLS/RandomForestsGLS.pdf

## Vignette
The package vignette, available at https://cran.rstudio.com/web/packages/RandomForestsGLS/vignettes/RandomForestsGLS_user_guide.pdf demonstrates with example how the functions available in `RandomForestsGLS` can be used for non-linear regression analysis of dependent data. Specific functions are discussed in much detail in the code documentation of the package. 

## Note
Some code blocks are borrowed from the R packages: spNNGP: Spatial Regression Models for Large Datasets using Nearest Neighbor Gaussian Processes https://CRAN.R-project.org/package=spNNGP and randomForest: Breiman and Cutler's Random Forests for Classification and Regression https://CRAN.R-project.org/package=randomForest 


## Citation
Please cite the following paper when you use RF-GLS

Saha, Arkajyoti, Sumanta Basu, and Abhirup Datta. "Random forests for spatially dependent data." Journal of the American Statistical Association (2021): 1-19.

