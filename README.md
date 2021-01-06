
# RandomForestsGLS
====

## Overview
RandomForestsGLS: Random Forests for dependent data: Fits non-linear regression models on dependent data with Generalised Least Square (GLS) based Random Forest (RF-GLS) detailed in Saha, Basu and Datta (2020) <https://arxiv.org/abs/2007.15421>.


## Installation
In order to install the development version of the package, please run the following command in R:

```{r }
if (!require("devtools")) install.packages("devtools")
devtools::install_github("ArkajyotiSaha/RandomForestsGLS", ref ="HEAD")
```

## Function description
For help on the functions in Brisc please use the following:
```{r }
?RFGLS_estimate_spatial #(for estimation in spatial data)
?RFGLS_estimate_timeseries #(for estimation in timeseries data)
?RFGLS_predict #(for prediction of mean functtion)
?RFGLS_predict_spatial #(for prediction of Spatial Response)
```

## Vignette
The package vignette, available at https://github.com/ArkajyotiSaha/RandomForestsGLS/blob/main/inst/doc/RandomForestsGLS_user_guide.pdf demonstrates with example how the functions available in `RandomForestsGLS` can be used for non-linear regression analysis of dependent data. Specific functions are discussed in much detail in the code documentation of the package. 

## Note
Some code blocks are borrowed from the R packages: spNNGP: Spatial Regression Models for Large Datasets using Nearest Neighbor Gaussian Processes https://CRAN.R-project.org/package=spNNGP and randomForest: Breiman and Cutler's Random Forests for Classification and Regression https://CRAN.R-project.org/package=randomForest 


## Citation
Please cite the following paper when you use RF-GLS

Saha, A., Basu, S., & Datta, A. (2020). Random Forests for dependent data. arXiv preprint arXiv:2007.15421.

