#' Conformal inference for interval outcomes
#'
#' \code{conformalInt} is a framework for weighted and unweighted conformal inference with interval
#' outcomes. It supports both weighted split conformal inference and weighted CV+,
#' including weighted Jackknife+ as a special case. For each type, it supports both conformalized
#' quantile regression (CQR) and standard conformal inference based on mean regression.
#'
#' @details The conformal interval for a testing point x is in the form of
#' \eqn{[\hat{m}^{L}(x) - \eta, \hat{m}^{R}(x) + \eta]} where \eqn{\hat{m}^{L}(x)} is fit by \code{lofun}
#' and \eqn{\hat{m}^{R}(x)} is fit by \code{upfun}.
#'
#' \code{lofun} or \code{upfun} can be a valid string, including
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
#' When \code{type = "CQR"}, \code{lofun} or \code{upfun} should also include an argument \code{quantiles} that is either
#' a vector of length 2 or a scalar, depending on the argument \code{side}. Other optional arguments can be
#' passed into \code{lofun} and \code{upfun} through \code{loparams} and \code{upparams}.
#'
#' @param X covariates.
#' @param Y interval outcomes, a matrix with two columns.
#' @param type a string that is either "CQR" or "mean".
#' @param lofun the function to fit the lower bound or a valid string; see Details.
#' @param loquantile the quantile to fit for \code{lofun}.
#' @param loparams a list of other parameters to be passed into \code{lofun}.
#' @param upfun the function to fit the upper bound or a valid string; see Details.
#' @param upquantile the quantile to fit for \code{upfun}.
#' @param upparams a list of other parameters to be passed into \code{upfun}.
#' @param wtfun NULL for unweighted conformal inference or a function for weighted conformal inference
#'              when \code{useCV = FALSE} and a list of functions for weighted conformal inference when \code{useCV = TRUE};
#'              see Details.
#' @param useCV FALSE for split conformal inference and TRUE for CV+.
#' @param trainprop proportion of units for training \code{outfun}.
#' @param trainid indices of training units; NULL by default and indices will be randomly sampled.
#' @param nfolds number of folds; 10 by default.
#' @param idlist list of indices of length \code{nfolds}; NULL by default and indices will be randomly sampled.
#'
#' @return a \code{conformalIntSplit} object when \code{useCV = FALSE} with the following attributes:
#' \itemize{
#' \item{Yscore:}{ a vector of non-conformity score on the calibration fold}
#' \item{wt:}{ a vector of weights on the calibration fold}
#' \item{Ymodel:}{ a function with required argument \code{X} that produces the estimates the conditional
#'                mean or quantiles of \code{X}}
#' \item{wtfun, type, loquantile, upquantile, trainprop, trainid:}{ the same as inputs}
#' }
#'
#' or a \code{conformalIntCV} object when \code{useCV = TRUE} with the following attributes:
#' \itemize{
#' \item{info: }{ a list of length \code{nfolds} with each element being a list with attributes
#'               \code{Yscore}, \code{wt} and \code{Ymodel} described above for each fold}
#' \item{wtfun, type, loquantile, upquantile, nfolds, idlist:}{ the same as inputs}
#' }
#'
#' @seealso
#' \code{\link{predict.conformalIntSplit}}, \code{\link{predict.conformalIntCV}}.
#'
#' @examples
#' \donttest{# Generate data from a linear model
#' set.seed(1)
#' n <- 1000
#' d <- 5
#' X <- matrix(rnorm(n * d), nrow = n)
#' beta <- rep(1, 5)
#' Ylo <- X %*% beta + rnorm(n)
#' Yup <- Ylo + pmax(1, 2 * rnorm(n))
#' Y <- cbind(Ylo, Yup)
#'
#' # Generate testing data
#' ntest <- 5
#' Xtest <- matrix(rnorm(ntest * d), nrow = ntest)
#'
#' # Run unweighted split CQR with the built-in quantile random forest learner
#' # grf package needs to be installed
#' obj <- conformalInt(X, Y, type = "CQR",
#'                     lofun = "quantRF", upfun = "quantRF",
#'                     wtfun = NULL, useCV = FALSE)
#' predict(obj, Xtest, alpha = 0.1)
#'
#' # Run unweighted standard split conformal inference with the built-in random forest learner
#' # randomForest package needs to be installed
#' obj <- conformalInt(X, Y, type = "mean",
#'                     lofun = "RF", upfun = "RF",
#'                     wtfun = NULL, useCV = FALSE)
#' predict(obj, Xtest, alpha = 0.1)
#'
#' # Run unweighted CQR-CV+ with the built-in quantile random forest learner
#' # grf package needs to be installed
#' obj <- conformalInt(X, Y, type = "CQR",
#'                     lofun = "quantRF", upfun = "quantRF",
#'                     wtfun = NULL, useCV = TRUE)
#' predict(obj, Xtest, alpha = 0.1)
#'
#' # Run unweighted standard CV+ with the built-in random forest learner
#' # randomForest package needs to be installed
#' obj <- conformalInt(X, Y, type = "mean",
#'                     lofun = "RF", upfun = "RF",
#'                     wtfun = NULL, useCV = TRUE)
#' predict(obj, Xtest, alpha = 0.1)
#'
#' # Run weighted split CQR with w(x) = pnorm(x1)
#' wtfun <- function(X){pnorm(X[, 1])}
#' obj <- conformalInt(X, Y, type = "CQR",
#'                    lofun = "quantRF", upfun = "quantRF",
#'                    wtfun = wtfun, useCV = FALSE)
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
#' obj <- conformalInt(X, Y, type = "CQR",
#'                     lofun = quantRF, upfun = quantRF,
#'                     wtfun = NULL, useCV = FALSE)
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
#' obj <- conformalInt(X, Y, type = "mean",
#'                     lofun = linearReg, upfun = linearReg,
#'                     wtfun = NULL, useCV = FALSE)
#' predict(obj, Xtest, alpha = 0.1)
#' }
#' @export
conformalInt <- function(X, Y,
                         type = c("CQR", "mean"),
                         lofun = NULL,
                         loquantile = 0.5,
                         loparams = list(),
                         upfun = NULL,
                         upquantile = 0.5,
                         upparams = list(),
                         wtfun = NULL,
                         useCV = FALSE,
                         trainprop = 0.75, trainid = NULL,
                         nfolds = 10, idlist = NULL){
    if (!is.matrix(Y) || ncol(Y) != 2){
        stop("Y must a matrix with 2 columns")
    }
    type <- type[1]
    stopifnot(type %in% c("CQR", "mean"))

    if (is.null(lofun)){
        lofun <- switch(type,
                         CQR = quantRF,
                         mean = RF)
    } else if (is.character(lofun)){
        lofun <- str_fun(lofun[1])
    } else if (is.function(lofun)){
        check_outfun(lofun, type)
    } else {
        stop("lofun must be NULL or a string or a function")
    }

    if (is.null(upfun)){
        upfun <- switch(type,
                        CQR = quantRF,
                        mean = RF)
    } else if (is.character(upfun)){
        upfun <- str_fun(upfun[1])
    } else if (is.function(upfun)){
        check_outfun(upfun, type)
    } else {
        stop("upfun must be NULL or a string or a function")
    }

    if (!useCV){
        conformalIntSplit(X, Y,
                          type,
                          lofun, loquantile, loparams,
                          upfun, upquantile, upparams,
                          wtfun,
                          trainprop, trainid)
    } else {
        conformalIntCV(X, Y,
                       type,
                       lofun, loquantile, loparams,
                       upfun, upquantile, upparams,
                       wtfun,
                       nfolds, idlist)
    }
}
