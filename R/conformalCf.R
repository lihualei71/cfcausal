#' Conformal counterfactual inference
#'
#' \code{conformalCf} computes intervals for counterfactuals or outcome with missing values in general.
#' It supports both split conformal inference and CV+,
#' including weighted Jackknife+ as a special case. For each type, it supports both conformalized
#' quantile regression (CQR) and standard conformal inference based on mean regression.
#'
#' @details The outcome \code{Y} must comprise both observed values and missing values encoded as NA.
#' The missing values are used to estimate the propensity score \eqn{P(missing | X)}.
#'
#' \code{estimand} controls the type of coverage to be guaranteed:
#' \itemize{
#' \item (default) when \code{estimand = "unconditional"}, the interval has
#' \eqn{P(Y \in \hat{C}(X))\ge 1 - \alpha};
#' \item when \code{estimand = "nonmissing"}, the interval has
#' \eqn{P(Y \in \hat{C}(X) | nonmissing) \ge 1 - \alpha};
#' \item when \code{estimand = "missing"}, the interval has
#' \eqn{P(Y \in \hat{C}(X) | missing) \ge 1 - \alpha}.
#' }
#'
#' @param X covariates.
#' @param Y outcome with missing values encoded as NA; see Details
#' @param estimand a string that is either "unconditional" or "nonmissing" or "missing"; see Details.
#' @param type a string that is either "CQR" or "mean".
#' @param side a string that is either "two" or "above" or "below".
#' @param quantiles for covariates in the training data; only necessary when \code{type = "CQR"}.
#' @param outfun a function that models the conditional mean or quantiles or a valid string; see Details.
#'               Default to be random forest when \code{type = "mean"} and quantile random forest when
#'               \code{type = "CQR"}.
#' @param outparams a list of other parameters to be passed into \code{outfun}.
#' @param psfun a function that models the missing mechanism (probability of missing given X); see Details.
#'              Default to be gradient boosting.
#' @param psparams a list of other parameters to be passed into \code{psfun}.
#' @param useCV FALSE for split conformal inference and TRUE for CV+.
#' @param trainprop proportion of units for training \code{outfun}.
#' @param nfolds number of folds; 10 by default.
#'
#' @return a \code{conformalSplit} object when \code{useCV = FALSE} or a \code{conformalCV} object
#' when \code{useCV = TRUE}. See \code{\link{conformal}} for details.
#'
#' @seealso
#' \code{\link{conformal}}, \code{\link{conformalIte}}
#'
#' @examples
#' \donttest{# Generate data from a linear model
#' set.seed(1)
#' n <- 1000
#' d <- 5
#' X <- matrix(rnorm(n * d), nrow = n)
#' beta <- rep(1, 5)
#' Y <- X %*% beta + rnorm(n)
#'
#' # Generate missing indicators
#' missing_prob <- pnorm(X[, 1])
#' if_missing <- missing_prob < runif(n)
#' Y[if_missing] <- NA
#'
#' # Generate testing data
#' ntest <- 5
#' Xtest <- matrix(rnorm(ntest * d), nrow = ntest)
#'
#' # Run weighted split CQR
#' obj <- conformalCf(X, Y, type = "CQR", quantiles = c(0.05, 0.95),
#'                    outfun = "quantRF", useCV = FALSE)
#' predict(obj, Xtest, alpha = 0.1)
#'
#' # Run weighted standard conformal inference
#' obj <- conformalCf(X, Y, type = "mean", quantiles = c(0.05, 0.95),
#'                    outfun = "RF", useCV = FALSE)
#' predict(obj, Xtest, alpha = 0.1)
#'
#' # Run split CQR with a self-defined quantile random forest
#' # Y, X, Xtest, quantiles should be included in the inputs
#' quantRF <- function(Y, X, Xtest, quantiles, ...){
#'     fit <- grf::quantile_forest(X, Y, quantiles = quantiles, ...)
#'     res <- predict(fit, Xtest, quantiles = quantiles)
#'     if (length(quantiles) == 1){
#'         res <- as.numeric(res)
#'     } else {
#'         res <- as.matrix(res)
#'     }
#'     return(res)
#' }
#' obj <- conformalCf(X, Y, type = "CQR", quantiles = c(0.05, 0.95),
#'                    outfun = quantRF, useCV = FALSE)
#' predict(obj, Xtest, alpha = 0.1)
#'
#' # Run standard split conformal inference with a self-defined linear regression
#' # Y, X, Xtest should be included in the inputs
#' linearReg <- function(Y, X, Xtest){
#'     data <- data.frame(Y = Y, X)
#'     fit <- lm(Y ~ ., data = data)
#'     as.numeric(predict(fit, Xtest))
#' }
#' X <- as.data.frame(X)
#' Xtest <- as.data.frame(Xtest)
#' obj <- conformalCf(X, Y, type = "mean", quantiles = c(0.05, 0.95),
#'                    outfun = linearReg, useCV = FALSE)
#' predict(obj, Xtest, alpha = 0.1)
#' }
#' @export
conformalCf <- function(X, Y,
                        estimand = c("unconditional",
                                     "nonmissing",
                                     "missing"),
                        type = c("CQR", "mean"),
                        side = c("two", "above", "below"),
                        quantiles = NULL,
                        outfun = NULL,
                        outparams = list(),
                        psfun = NULL,
                        psparams = list(),
                        useCV = FALSE,
                        trainprop = 0.75,
                        nfolds = 10){
    if (!useCV){
        conformalCf_split(X, Y,
                          estimand,
                          type, side,
                          quantiles,
                          outfun, outparams,
                          psfun, psparams,
                          trainprop)
    } else {
        conformalCf_CV(X, Y,
                       estimand,
                       type, side,
                       quantiles,
                       outfun, outparams,
                       psfun, psparams,
                       nfolds)
    }
}
