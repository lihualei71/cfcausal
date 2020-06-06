#' Generic Conformal inference
#'
#' \code{conformal} is a framework for weighted and unweighted conformal inference with continuous
#' outcomes. It supports both weighted split conformal inference and weighted CV+,
#' including weighted Jackknife+ as a special case. For each type, it supports both conformalized
#' quantile regression (CQR) and standard conformal inference based on mean regression.
#'
#'  @details When \code{side = "two"}, CQR produces intervals in the form of
#' \deqn{[q_{\alpha_{lo}}(x) - \eta, q_{\alpha_{hi}}(x) + \eta]}
#' where \eqn{q_{\alpha_{lo}}(x)} and \eqn{q_{\alpha_{hi}}(x)} are estimates of conditional
#' quantiles of Y given X and the standard conformal inference produces (two-sided) intervals in the form of
#' \deqn{[m(x) - \eta, m(x) + \eta]}
#' where \eqn{m(x)} is an estimate of conditional mean/median of Y given X. When \code{side = "above"},
#' intervals are of form [-Inf, a(x)] and when \code{side = "below"} the intervals are of form [a(x), Inf].
#'
#' \code{quantiles} should be given when \code{type = "CQR"}. When \code{side = "two"}, \code{quantiles}
#' should be a vector of length 2, giving \eqn{\alpha_{lo}} and \eqn{\alpha_{hi}}. When \code{side = "above"}
#' or \code{side = "below"}, only one quantile should be given.
#'
#' \code{outfun} can be a valid string, including
#' \itemize{
#' \item "RF" for random forest that predicts the conditional mean, a wrapper from \code{randomForest} package.
#'   Used when \code{type = "mean"};
#' \item "quantRF" for quantile random forest that predicts the conditional quantiles, a wrapper from
#'   \code{grf} package. Used when \code{type = "CQR"};
#' \item "Boosting" for gradient boosting that predicts the conditional mean, a wrapper from \code{gbm}
#'    package. Used when \code{type = "mean"};
#' \item "quantBoosting" for quantile gradient boosting that predicts the conditional quantiles, a wrapper from
#'   \code{gbm} package. Used when \code{type = "CQR"};
#' \item "BART" for gradient boosting that predicts the conditional mean, a wrapper from \code{bartMachine}
#'    package. Used when \code{type = "mean"};
#' \item "quantBART" for quantile gradient boosting that predicts the conditional quantiles, a wrapper from
#'   \code{bartMachine} package. Used when \code{type = "CQR"};
#' }
#'
#' or a function object whose input must include, but not limited to
#' \itemize{
#' \item \code{Y} for outcome in the training data;
#' \item \code{X} for covariates in the training data;
#' \item \code{Xtest} for covariates in the testing data.
#' }
#' When \code{type = "CQR"}, \code{outfun} should also include an argument \code{quantiles} that is either
#' a vector of length 2 or a scalar, depending on the argument \code{side}. Other optional arguments can be
#' passed into \code{outfun} through \code{outparams}.
#'
#' \code{wtfun} is NULL for unweighted conformal inference. For weighted split conformal inference, it is a
#' function with required input \code{X} that produces a vector of numeric values in [0, 1] of length \code{nrow(X)}.
#' For weighted CV+, it can be a function as in the case \code{useCV = FALSE} so that the same function will
#' apply to each fold or a list of functions of length \code{nfolds} so that \code{wtfun[[k]]} is applied to fold \code{k}.
#'
#' @param X covariates.
#' @param Y outcome.
#' @param type a string that is either "CQR" or "mean".
#' @param side a string that is either "two" or "above" or "below".
#' @param quantiles for covariates in the training data; only necessary when \code{type = "CQR"}.
#' @param outfun a function that models the conditional mean or quantiles or a valid string; see Details.
#'               Default to be random forest when \code{type = "mean"} and quantile random forest when
#'               \code{type = "CQR"}.
#' @param outparams a list of other parameters to be passed into \code{outfun}.
#' @param wtfun NULL for unweighted conformal inference or a function for weighted conformal inference
#'              when \code{useCV = FALSE} and a list of functions for weighted conformal inference when \code{useCV = TRUE};
#'              see Details.
#' @param useCV FALSE for split conformal inference and TRUE for CV+.
#' @param trainprop proportion of units for training \code{outfun}.
#' @param trainid indices of training units; NULL by default and indices will be randomly sampled.
#' @param nfolds number of folds; 10 by default.
#' @param idlist list of indices of length \code{nfolds}; NULL by default and indices will be randomly sampled.
#'
#' @return a \code{conformalSplit} object when \code{useCV = FALSE} with the following attributes:
#' \itemize{
#' \item{Yscore:}{ a vector of non-conformity score on the calibration fold}
#' \item{wt:}{ a vector of weights on the calibration fold}
#' \item{Ymodel:}{ a function with required argument \code{X} that produces the estimates the conditional
#'                mean or quantiles of \code{X}}
#' \item{wtfun, type, side, quantiles, trainprop, trainid:}{ the same as inputs}
#' }
#'
#' or a \code{conformalCV} object when \code{useCV = TRUE} with the following attributes:
#' \itemize{
#' \item{info: }{ a list of length \code{nfolds} with each element being a list with attributes
#'               \code{Yscore}, \code{wt} and \code{Ymodel} described above for each fold}
#' \item{wtfun, type, side, quantiles, nfolds, idlist:}{ the same as inputs}
#' }
#'
#' @seealso
#' \code{\link{predict.conformalSplit}}, \code{\link{predict.conformalCV}}.
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
#' # Generate testing data
#' ntest <- 5
#' Xtest <- matrix(rnorm(ntest * d), nrow = ntest)
#'
#' # Run unweighted split CQR with the built-in quantile random forest learner
#' # grf package needs to be installed
#' obj <- conformal(X, Y, type = "CQR", quantiles = c(0.05, 0.95),
#'                  outfun = "quantRF", wtfun = NULL, useCV = FALSE)
#' predict(obj, Xtest, alpha = 0.1)
#'
#' # Run unweighted standard split conformal inference with the built-in random forest learner
#' # randomForest package needs to be installed
#' obj <- conformal(X, Y, type = "mean",
#'                  outfun = "RF", wtfun = NULL, useCV = FALSE)
#' predict(obj, Xtest, alpha = 0.1)
#'
#' # Run unweighted CQR-CV+ with the built-in quantile random forest learner
#' # grf package needs to be installed
#' obj <- conformal(X, Y, type = "CQR", quantiles = c(0.05, 0.95),
#'                  outfun = "quantRF", wtfun = NULL, useCV = TRUE)
#' predict(obj, Xtest, alpha = 0.1)
#'
#' # Run unweighted standard CV+ with the built-in random forest learner
#' # randomForest package needs to be installed
#' obj <- conformal(X, Y, type = "mean",
#'                  outfun = "RF", wtfun = NULL, useCV = TRUE)
#' predict(obj, Xtest, alpha = 0.1)
#'
#' # Run weighted split CQR with w(x) = pnorm(x1)
#' wtfun <- function(X){pnorm(X[, 1])}
#' obj <- conformal(X, Y, type = "CQR", quantiles = c(0.05, 0.95),
#'                  outfun = "quantRF", wtfun = wtfun, useCV = FALSE)
#' predict(obj, Xtest, alpha = 0.1)
#'
#' # Run unweighted split CQR with a self-defined quantile random forest
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
#' obj <- conformal(X, Y, type = "CQR", quantiles = c(0.05, 0.95),
#'                  outfun = quantRF, wtfun = NULL, useCV = FALSE)
#' predict(obj, Xtest, alpha = 0.1)
#'
#' # Run unweighted standard split conformal inference with a self-defined linear regression
#' # Y, X, Xtest should be included in the inputs
#' linearReg <- function(Y, X, Xtest){
#'     data <- data.frame(Y = Y, X)
#'     fit <- lm(Y ~ ., data = data)
#'     as.numeric(predict(fit, Xtest))
#' }
#' X <- as.data.frame(X)
#' Xtest <- as.data.frame(Xtest)
#' obj <- conformal(X, Y, type = "mean",
#'                  outfun = linearReg, wtfun = NULL, useCV = FALSE)
#' predict(obj, Xtest, alpha = 0.1)
#' }
#' @export
conformal <- function(X, Y,
                      type = c("CQR", "mean"),
                      side = c("two", "above", "below"),
                      quantiles = NULL,
                      outfun = NULL,
                      outparams = list(),
                      wtfun = NULL,
                      useCV = FALSE,
                      trainprop = 0.75, trainid = NULL,
                      nfolds = 10, idlist = NULL){
    type <- type[1]
    stopifnot(type %in% c("CQR", "mean"))
    side <- side[1]
    stopifnot(side %in% c("two", "above", "below"))

    if (is.null(outfun)){
        outfun <- switch(type,
                         CQR = quantRF,
                         mean = RF)
    } else if (is.character(outfun)){
        outfun <- str_fun(outfun[1])
    } else if (is.function(outfun)){
        check_outfun(outfun, type)
    } else {
        stop("outfun must be NULL or a string or a function")
    }

    if (!useCV){
        conformalSplit(X, Y,
                       type, side,
                       quantiles,
                       outfun, outparams,
                       wtfun,
                       trainprop, trainid)
    } else {
        conformalCV(X, Y,
                    type, side,
                    quantiles,
                    outfun, outparams,
                    wtfun,
                    nfolds, idlist)
    }
}
