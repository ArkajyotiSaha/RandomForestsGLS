
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
A vignette for the package that describes how to use RF-GLS for estimating non-linear effects under dependent errors is avilable at https://github.com/ArkajyotiSaha/Article/blob/main/docs/RandomForestsGLS_user_guide.pdf

## Note
Some code blocks are borrowed from the R packages: spNNGP: Spatial Regression Models for Large Datasets using Nearest Neighbor Gaussian Processes https://CRAN.R-project.org/package=spNNGP and randomForest: Breiman and Cutler's Random Forests for Classification and Regression https://CRAN.R-project.org/package=randomForest 


## Citation
Please cite the following paper when you use RF-GLS

Saha, A., Basu, S., & Datta, A. (2020). Random Forests for dependent data. arXiv preprint arXiv:2007.15421.

