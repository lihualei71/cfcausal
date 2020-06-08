## -----------------------------------------------------------------------------
# Install the "cfcausal" package from github.
# if (!require("devtools")){
#     install.packages("devtools")
# }
# devtools::install_github("lihualei71/cfcausal")
library("cfcausal")

## -----------------------------------------------------------------------------
# Generate data
set.seed(2020)
genY <- function(X){
    2 / (1 + exp(-12 * (X[, 1] - 0.5))) * 2 / (1 + exp(-12 * (X[, 2] - 0.5))) + rnorm(n)
}
n <- 1000
d <- 10
X <- matrix(runif(n * d), nrow = n, ncol = d)
Y <- genY(X)
ps <- (1 + pbeta(X[, 1], 2, 4)) / 4
T <- as.numeric(ps < runif(n))
Y[!T] <- NA
summary(Y)

## -----------------------------------------------------------------------------
obj <- conformalCf(X, Y, type = "CQR",
                   quantiles = c(0.05, 0.95),
                   outfun = "quantRF", useCV = FALSE)
class(obj)

## -----------------------------------------------------------------------------
ntest <- 5
Xtest <- matrix(runif(ntest * d), nrow = ntest, ncol = d)

## -----------------------------------------------------------------------------
CI <- predict(obj, Xtest, alpha = 0.1)
CI

## -----------------------------------------------------------------------------
ntest_large <- 10000
Xtest_large <- matrix(runif(ntest_large * d), nrow = ntest_large, ncol = d)
Ytest_large <- genY(Xtest_large)
CI_large <- predict(obj, Xtest_large, alpha = 0.1)
mean(CI_large[, 1] <= Ytest_large & CI_large[, 2] >= Ytest_large)

## -----------------------------------------------------------------------------
obj <- conformalCf(X, Y, type = "CQR",
                   quantiles = c(0.05, 0.95),
                   outfun = "quantRF", useCV = TRUE,
                   nfolds = 10)
predict(obj, Xtest, alpha = 0.1)

## -----------------------------------------------------------------------------
obj <- conformalCf(X, Y, type = "mean",
                   quantiles = c(0.05, 0.95),
                   outfun = "RF", useCV = FALSE)
predict(obj, Xtest, alpha = 0.1)

## -----------------------------------------------------------------------------
obj <- conformalCf(X, Y, type = "CQR", side = "above",
                   quantiles = 0.95,
                   outfun = "quantRF", useCV = FALSE)
predict(obj, Xtest, alpha = 0.1)
obj <- conformalCf(X, Y, type = "CQR", side = "below",
                   quantiles = 0.05,
                   outfun = "quantRF", useCV = FALSE)
predict(obj, Xtest, alpha = 0.1)

## -----------------------------------------------------------------------------
# Install grf package
if (!require("grf")){
    install.packages("grf")
}
# User-defined quantile random forest
quantRF <- function(Y, X, Xtest, quantiles, ...){
    fit <- grf::quantile_forest(X, Y, quantiles = quantiles, ...)
    res <- predict(fit, Xtest, quantiles = quantiles)
    if (length(quantiles) == 1){
        res <- as.numeric(res)
    } else {
        res <- as.matrix(res)
    }
    return(res)
}
# conformalCf with user-defined quantRF
obj <- conformalCf(X, Y, type = "CQR",
                   quantiles = c(0.05, 0.95),
                   outfun = quantRF, useCV = FALSE)
predict(obj, Xtest, alpha = 0.1)

## -----------------------------------------------------------------------------
linearReg <- function(Y, X, Xtest){
    X <- as.data.frame(X)
    Xtest <- as.data.frame(Xtest)
    data <- data.frame(Y = Y, X)
    fit <- lm(Y ~ ., data = data)
    as.numeric(predict(fit, Xtest))
}
obj <- conformalCf(X, Y, type = "mean", 
                   outfun = linearReg, useCV = FALSE)
predict(obj, Xtest, alpha = 0.1)

## -----------------------------------------------------------------------------
logitReg <- function(Y, X, Xtest, ...){
    X <- as.data.frame(X)
    Xtest <- as.data.frame(Xtest)
    data <- data.frame(Y = Y, X)
    fit <- glm(Y ~ ., data = data, family = "binomial", ...)
    as.numeric(predict(fit, Xtest, type = "response"))
}
obj <- conformalCf(X, Y, type = "CQR", quantiles = c(0.05, 0.95),
                   outfun = "quantRF", psfun = logitReg, useCV = FALSE)
predict(obj, Xtest, alpha = 0.1)

## -----------------------------------------------------------------------------
# Generate training data
set.seed(2020)
genY <- function(X){
    2 / (1 + exp(-12 * (X[, 1] - 0.5))) * 2 / (1 + exp(-12 * (X[, 2] - 0.5))) + rnorm(n)
}
n <- 1000
d <- 10
X <- matrix(runif(n * d), nrow = n, ncol = d)
Y1 <- genY(X)
Y0 <- rnorm(n)
ps <- (1 + pbeta(X[, 1], 2, 4)) / 4
T <- as.numeric(ps < runif(n))
Y <- ifelse(T == 1, Y1, Y0)

# Generate testing data
ntest <- 5
Xtest <- matrix(runif(ntest * d), nrow = ntest, ncol = d)
pstest <- (1 + pbeta(Xtest[, 1], 2, 4)) / 4
Ttest <- as.numeric(pstest < runif(ntest))
Y1test <- genY(Xtest)
Y0test <- rnorm(ntest)
Ytest <- ifelse(Ttest == 1, Y1test, Y0test)

## -----------------------------------------------------------------------------
# Inexact nested method
CIfun <- conformalIte(X, Y, T, alpha = 0.1, algo = "nest", exact = FALSE, type = "CQR",
                      quantiles = c(0.05, 0.95), outfun = "quantRF", useCV = FALSE)
CIfun(Xtest)

## -----------------------------------------------------------------------------
# Exact nested method
CIfun <- conformalIte(X, Y, T, alpha = 0.1, algo = "nest", exact = TRUE, type = "CQR",
                      quantiles = c(0.05, 0.95), outfun = "quantRF",  useCV = FALSE)
CIfun(Xtest)

