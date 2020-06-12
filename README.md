# cfcausal
An R package for conformal inference of counterfactuals and individual treatment effects

## Overview
This R package implements weighted conformal inference-based procedures for counterfactuals and individual treatment effects proposed in our paper: [Conformal Inference of Counterfactuals and Individual Treatment Effects](https://arxiv.org/abs/2006.06138). It includes both the split conformal inference and cross-validation+. For each type of conformal inference, both conformalized quantile regression (CQR) and standard conformal inference are supported. It provides a pool of convenient learners and allows flexible user-defined learners for conditional mean and quantiles. 

- `conformalCf()` produces intervals for counterfactuals or outcomes with missing values in general.
- `conformalIte()` produces intervals for individual treatment effects with a binary treatment under the potential outcome framework. 
- `conformal()` provides a generic framework of weighted conformal inference for continuous outcomes.
- `conformalInt()` provides a generic framework of weighted conformal inference for interval outcomes.

## Installation         

```
if (!require("devtools")){
    install.packages("devtools")
}
devtools::install_github("lihualei71/cfcausal")
```
To access the vignette, run the following code to build it. 
```
devtools::install_github("lihualei71/cfcausal", build_vignettes = TRUE)
```

We suggest installing [grf](https://cran.r-project.org/web/packages/grf/grf.pdf), [randomForest](https://cran.r-project.org/web/packages/randomForest/randomForest.pdf), [gbm](https://cran.r-project.org/web/packages/gbm/gbm.pdf) and [bartMachine](https://cran.r-project.org/web/packages/bartMachine/bartMachine.pdf) to take advantage of the built-in learners. 

## Usage Examples
We illustrate the usage of cfcausal package using simple synthetic datasets. For details please read the vignette (`vignette("cfcausal_demo", package = "cfcausal")`) and the manual.

```
#### Conformal inference of counterfactuals
library("cfcausal")

# Generate data
set.seed(1)
n <- 1000
d <- 5
X <- matrix(rnorm(n * d), nrow = n)
beta <- rep(1, 5)
Y <- X %*% beta + rnorm(n)

# Generate missing indicators
missing_prob <- pnorm(X[, 1])
if_missing <- missing_prob < runif(n)
Y[if_missing] <- NA

# Generate testing data
ntest <- 5
Xtest <- matrix(rnorm(ntest * d), nrow = ntest)

# Run weighted split CQR
obj <- conformalCf(X, Y, type = "CQR", 
                   quantiles = c(0.05, 0.95),
                   outfun = "quantRF", useCV = FALSE)
predict(obj, Xtest, alpha = 0.1)

# Run weighted CQR-CV+
obj <- conformalCf(X, Y, type = "CQR", 
                   quantiles = c(0.05, 0.95),
                   outfun = "quantRF", useCV = TRUE)
predict(obj, Xtest, alpha = 0.1)
```

```
#### Conformal inference of individual treatment effects
library("cfcausal")

# Generate potential outcomes from two linear models
set.seed(1)
n <- 1000
d <- 5
X <- matrix(rnorm(n * d), nrow = n)
beta <- rep(1, 5)
Y1 <- X %*% beta + rnorm(n)
Y0 <- rnorm(n)

# Generate treatment indicators
ps <- pnorm(X[, 1])
T <- as.numeric(ps < runif(n))
Y <- ifelse(T == 1, Y1, Y0)

# Generate testing data
ntest <- 5
Xtest <- matrix(rnorm(ntest * d), nrow = ntest)

# Inexact nested method
CIfun <- conformalIte(X, Y, T, alpha = 0.1, 
                      algo = "nest", exact = FALSE, type = "CQR",
                      quantiles = c(0.05, 0.95), 
                      outfun = "quantRF", useCV = FALSE)
CIfun(Xtest)

# Exact nested method
CIfun <- conformalIte(X, Y, T, alpha = 0.1, 
                      algo = "nest", exact = TRUE, type = "CQR",
                      quantiles = c(0.05, 0.95), 
                      outfun = "quantRF",  useCV = FALSE)
CIfun(Xtest)

# Naive method
CIfun <- conformalIte(X, Y, T, alpha = 0.1, 
                      algo = "naive", type = "CQR",
                      quantiles = c(0.05, 0.95), 
                      outfun = "quantRF",  useCV = FALSE)
CIfun(Xtest)

# Counterfactual method, Y and T needs to be observed
pstest <- pnorm(Xtest[, 1])
Ttest <- as.numeric(pstest < runif(ntest))
Y1test <- Xtest %*% beta + rnorm(ntest)
Y0test <- rnorm(ntest)
Ytest <- ifelse(Ttest == 1, Y1test, Y0test)
CIfun <- conformalIte(X, Y, T, alpha = 0.1, 
                      algo = "counterfactual", type = "CQR",
                      quantiles = c(0.05, 0.95), 
                      outfun = "quantRF",  useCV = FALSE)
CIfun(Xtest, Ytest, Ttest)
```
