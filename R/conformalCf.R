#' Conformal inference for counterfactuals
#'
#' \code{conformalCf} computes intervals for counterfactuals or outcomes with ignorable missing values in general.
#' It supports both split conformal inference and CV+,
#' including weighted Jackknife+ as a special case. For each type, it supports both conformalized
#' quantile regression (CQR) and standard conformal inference based on conditional mean regression.
#'
#' @details The outcome \code{Y} must comprise both observed values and missing values encoded as NA.
#' The missing values are used to estimate the propensity score \eqn{P(missing | X)}.
#'
#' \code{estimand} controls the type of coverage to be guaranteed:
#' \itemize{
#' \item (Default) when \code{estimand = "unconditional"}, the interval has
#' \eqn{P(Y \in \hat{C}(X))\ge 1 - \alpha}.
#' \item When \code{estimand = "nonmissing"}, the interval has
#' \eqn{P(Y \in \hat{C}(X) | nonmissing) \ge 1 - \alpha}.
#' \item When \code{estimand = "missing"}, the interval has
#' \eqn{P(Y \in \hat{C}(X) | missing) \ge 1 - \alpha}.
#' }
#'
#' When \code{side = "above"},
#' intervals are of form [-Inf, a(x)] and when \code{side = "below"} the intervals are of form [a(x), Inf].
#' 
#' \code{outfun} can be a valid string, including
#' \itemize{
#' \item "RF" for random forest that predicts the conditional mean, a wrapper built on \code{randomForest} package.
#'   Used when \code{type = "mean"}.
#' \item "quantRF" for quantile random forest that predicts the conditional quantiles, a wrapper built on
#'   \code{grf} package. Used when \code{type = "CQR"}.
#' \item "Boosting" for gradient boosting that predicts the conditional mean, a wrapper built on \code{gbm}
#'    package. Used when \code{type = "mean"}.
#' \item "quantBoosting" for quantile gradient boosting that predicts the conditional quantiles, a wrapper built on
#'   \code{gbm} package. Used when \code{type = "CQR"}.
#' \item "BART" for gradient boosting that predicts the conditional mean, a wrapper built on \code{bartMachine}
#'    package. Used when \code{type = "mean"}.
#' \item "quantBART" for quantile gradient boosting that predicts the conditional quantiles, a wrapper built on
#'   \code{bartMachine} package. Used when \code{type = "CQR"}.
#' }
#' or a function object whose input must include, but not limited to
#' \itemize{
#' \item \code{Y} for outcome in the training data.
#' \item \code{X} for covariates in the training data.
#' \item \code{Xtest} for covariates in the testing data.
#' }
#' When \code{type = "CQR"}, \code{outfun} should also include an argument \code{quantiles} that is either
#' a vector of length 2 or a scalar, depending on the argument \code{side}. The output of \code{outfun}
#' must be a matrix with two columns giving the conditional quantile estimates when \code{quantiles} is
#' a vector of length 2; otherwise, it must be a vector giving the conditional quantile estimate or
#' conditional mean estimate. Other optional arguments can be passed into \code{outfun} through \code{outparams}.
#'
#' \code{psfun} can be a valid string, including
#' \itemize{
#' \item "RF" for random forest that predicts the propensity score, a wrapper built on \code{randomForest} package.
#'   Used when \code{type = "mean"}.
#' \item "Boosting" for gradient boosting that predicts the propensity score, a wrapper built on \code{gbm}
#'    package. Used when \code{type = "mean"}.
#' }
#' or a function object whose input must include, but not limited to
#' \itemize{
#' \item \code{Y} for treatment assignment, a binary vector, in the training data.
#' \item \code{X} for covariates in the training data.
#' \item \code{Xtest} for covariates in the testing data.
#' }
#' The output of \code{psfun} must be a vector of predicted probabilities. Other optional arguments
#' can be passed into \code{psfun} through \code{psparams}.
#' 
#' @param X covariates.
#' @param Y outcome vector with missing values encoded as NA. See Details.
#' @param estimand a string that takes values in \{"unconditional", "nonmissing", "missing"\}. See Details.
#' @param type a string that takes values in \{"CQR", "mean"\}.
#' @param side a string that takes values in \{"two", "above", "below"\}. See Details.
#' @param quantiles a scalar or a vector of length 2 depending on \code{side}. Used only when \code{type = "CQR"}. See Details. 
#' @param outfun a function that models the conditional mean or quantiles, or a valid string. 
#'               The default is random forest when \code{type = "mean"} and quantile random forest when
#'               \code{type = "CQR"}. See Details.
#' @param outparams a list of other parameters to be passed into \code{outfun}.
#' @param psfun a function that models the missing mechanism (probability of missing given X), or a valid string.
#'              The default is "Boosting". See Details. 
#' @param psparams a list of other parameters to be passed into \code{psfun}.
#' @param useCV FALSE for split conformal inference and TRUE for CV+.
#' @param trainprop proportion of units for training \code{outfun}. The default if 75\%. Used only when \code{useCV = FALSE}.
#' @param nfolds number of folds. The default is 10. Used only when \code{useCV = TRUE}. 
#'
#' @return a \code{conformalSplit} object when \code{useCV = FALSE} or a \code{conformalCV} object
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
#' obj <- conformalCf(X, Y, type = "mean",
#'                    outfun = "RF", useCV = FALSE)
#' predict(obj, Xtest, alpha = 0.1)
#'
#' # Run one-sided weighted split CQR
#' obj1 <- conformalCf(X, Y, type = "CQR", side = "above",
#'                     quantiles = 0.95, outfun = "quantRF", useCV = FALSE)
#' predict(obj1, Xtest, alpha = 0.1)
#' obj2 <- conformalCf(X, Y, type = "CQR", side = "below",
#'                     quantiles = 0.05, outfun = "quantRF", useCV = FALSE)
#' predict(obj2, Xtest, alpha = 0.1)
#'
#' # Run split CQR with a self-defined quantile random forest
#' # Y, X, Xtest, quantiles should be included in the inputs
#' quantRF <- function(Y, X, Xtest, quantiles, ...){
#'     fit <- grf::quantile_forest(X, Y, quantiles = quantiles, ...)
#'     res <- predict(fit, Xtest, quantiles = quantiles)
#'     if (is.list(res) && !is.data.frame(res)){
#'     # for the recent update of \code{grf} package that
#'     # changes the output format
#'         res <- res$predictions 
#'     }
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
#'     X <- as.data.frame(X)
#'     Xtest <- as.data.frame(Xtest)
#'     data <- data.frame(Y = Y, X)
#'     fit <- lm(Y ~ ., data = data)
#'     as.numeric(predict(fit, Xtest))
#' }
#' obj <- conformalCf(X, Y, type = "mean",
#'                    outfun = linearReg, useCV = FALSE)
#' predict(obj, Xtest, alpha = 0.1)
#'
#' # Run split CQR with a built-in psfun
#' # Y, X, Xtest, should be included in the inputs
#' obj <- conformalCf(X, Y, type = "CQR", quantiles = c(0.05, 0.95),
#'                    outfun = "quantRF", psfun = "RF", useCV = FALSE)
#' predict(obj, Xtest, alpha = 0.1)
#' 
#' # Run split CQR with a self-defined function to estimate propensity scores
#' # Y, X, Xtest, should be included in the inputs
#' logitReg <- function(Y, X, Xtest, ...){
#'     X <- as.data.frame(X)
#'     Xtest <- as.data.frame(Xtest)
#'     data <- data.frame(Y = Y, X)
#'     fit <- glm(Y ~ ., data = data, family = "binomial", ...)
#'     as.numeric(predict(fit, Xtest, type = "response"))
#' }
#' obj <- conformalCf(X, Y, type = "CQR", quantiles = c(0.05, 0.95),
#'                    outfun = "quantRF", psfun = logitReg, useCV = FALSE)
#' predict(obj, Xtest, alpha = 0.1)
#'
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
    type <- type[1]
    stopifnot(type %in% c("CQR", "mean"))
    side <- side[1]
    stopifnot(side %in% c("two", "above", "below"))
    estimand <- estimand[1]
    stopifnot(estimand %in% c("unconditional",
                              "nonmissing",
                              "missing"))

    if (is.null(outfun)){
        outfun <- switch(type,
                         CQR = quantRF,
                         mean = RF)
    } else if (is.character(outfun)){
        outfun <- str_outfun(outfun[1])
    } else if (is.function(outfun)){
        check_outfun(outfun, type)
    } else {
        stop("outfun must be NULL or a string or a function")
    }
    
    if (is.null(psfun)){
        psfun <- Boosting
    } else if (is.character(psfun)){
        psfun <- str_psfun(psfun[1])
    } else if (is.function(psfun)){
        check_psfun(psfun)
    } else {
        stop("psfun must be NULL or a string or a function")
    }
    
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
