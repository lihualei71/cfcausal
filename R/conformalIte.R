#' Conformal inference for individual treatment effects
#'
#' \code{conformalIte} supports four algorithms: the nested approach, with exact and inexact
#' calibration, and the naive approach for cases with both potential outcomes missing and the counterfactual
#' inference for cases with only one potential outcome missing. For each algorithm, it supports both
#' split conformal inference and CV+, including weighted Jackknife+ as a special case. For each type, it
#' supports both conformalized quantile regression (CQR) and standard conformal inference based on mean regression.
#'
#' @details The algorithm to be used is controlled by \code{algo} and \code{exact}:
#' \itemize{
#' \item (Default) when \code{algo = "nest"} and \code{exact = FALSE}, the inexact nested approach is used. It
#' first split the data into two folds, with the second fold including \code{cfprop} fraction of units. Then it applies
#' \code{conformalCf} on the first fold to compute counterfactual intervals on the second fold, which further yields
#' interval estimates of ITE \eqn{\hat{C}(X_i)}. It then fits \eqn{\hat{C}^{L}(X_i)} and \eqn{\hat{C}^{R}(X_i)} on \eqn{X_i}'s
#' \item When \code{algo = "nest"} and \code{exact = TRUE}, the exact nested approach is used. It has the same step to produce
#' ITE intervals \eqn{\hat{C}(X_i)}'s on the second fold but then applies \code{\link{conformalInt}} to calibrate them
#' \item When \code{algo = "naive"}, the naive approach is used. It applies \code{\link{conformalCf}} on the data and
#' produce counterfactual intervals for both Y(1) and Y(0). The ITE intervals are computed by contrasting two counterfactual
#' intervals.
#' \item When \code{algo = "counterfactual"}, it handles the case where the treatment assignments and the observed outcome are
#' both available for each testing point. As with the naive approach, it applies \code{\link{conformalCf}} on the data and
#' produce counterfactual intervals for both Y(1) and Y(0). The ITE intervals are computed by contrasting the observed outcome
#' and the interval for the missing potential outcome.
#' }
#'
#' @param X covariates.
#' @param Y observed outcome.
#' @param T treatment indicator, a binary vector.
#' @param alpha confidence level.
#' @param algo a string that is either "nest" or "naive" or "counterfactual"; see Details.
#' @param exact a logical indicating whether the exact calibration is used for nested approach;
#' used only when \code{algo = "nest"}.
#' @param type a string that is either "CQR" or "mean".
#' @param side a string that is either "two" or "above" or "below".
#' @param quantiles for covariates in the training data; only necessary when \code{type = "CQR"}.
#' @param outfun a function that models the conditional mean or quantiles or a valid string; see \code{\link{conformal}}.
#'               Default to be random forest when \code{type = "mean"} and quantile random forest when
#'               \code{type = "CQR"}.
#' @param outparams a list of other parameters to be passed into \code{outfun}.
#' @param psfun a function that models the missing mechanism (probability of missing given X); see Details.
#'               Default to be gradient boosting.
#' @param psparams a list of other parameters to be passed into \code{psfun}.
#' @param cfprop the proportion of units to be used to compute ITE intervals in nested approach; used only when
#' \code{algo = "nest"}.
#' @param citype the type of interval conformal inference used in the nested approach with exact calibration;
#' used only when \code{algo = "nest"} and \code{exact = TRUE}.
#' @param lofun the function to fit the lower bound or a valid string; see \code{\link{conformalInt}}; used only when
#' \code{algo = "nest"}.
#' @param loquantile the quantile to fit for \code{lofun}; see \code{\link{conformalInt}}; used only when
#' \code{algo = "nest"}, \code{exact = TRUE} and \code{citype = "CQR"}.
#' @param loparams a list of other parameters to be passed into \code{lofun}.
#' @param upfun the function to fit the upper bound or a valid string; see \code{\link{conformalInt}}; used only when
#' \code{algo = "nest"}.
#' @param upquantile the quantile to fit for \code{upfun}; see \code{\link{conformalInt}}; used only when
#' \code{algo = "nest"}, \code{exact = TRUE} and \code{citype = "CQR"}.
#' @param upparams a list of other parameters to be passed into \code{upfun}.
#' @param useCV FALSE for split conformal inference and TRUE for CV+.
#' @param trainprop proportion of units for training \code{outfun}.
#' @param nfolds number of folds; 10 by default.
#' @param wthigh upper truncation level of weights; see \code{\link{predict.conformalSplit}} or \code{\link{predict.conformalCV}}.
#' @param wtlow lower truncation level of weights; see \code{\link{predict.conformalSplit}} or \code{\link{predict.conformalCV}}.
#'
#' @return a function that outputs the interval estimates on a given dataset. When \code{algo = "nest"} or \code{"naive"}, it takes
#' a single input \code{X}; when \code{algo = "counterfactual"}, it takes three inputs \code{X}, \code{Y} and \code{T}.
#'
#' #' @seealso
#' \code{\link{conformal}}, \code{\link{conformalInt}}, \code{\link{conformalCf}}
#'
#' @examples
#' \donttest{# Generate potential outcomes from two linear models
#' set.seed(1)
#' n <- 1000
#' d <- 5
#' X <- matrix(rnorm(n * d), nrow = n)
#' beta <- rep(1, 5)
#' Y1 <- X %*% beta + rnorm(n)
#' Y0 <- rnorm(n)
#'
#' # Generate treatment indicators
#' ps <- pnorm(X[, 1])
#' T <- as.numeric(ps < runif(n))
#' Y <- ifelse(T == 1, Y1, Y0)
#'
#' # Generate testing data
#' ntest <- 5
#' Xtest <- matrix(rnorm(ntest * d), nrow = ntest)
#'
#' # Inexact nested method
#' CIfun <- conformalIte(X, Y, T, alpha = 0.1, algo = "nest", exact = FALSE, type = "CQR",
#'                       quantiles = c(0.05, 0.95), outfun = "quantRF", useCV = FALSE)
#' CIfun(Xtest)
#'
#' # Exact nested method
#' CIfun <- conformalIte(X, Y, T, alpha = 0.1, algo = "nest", exact = TRUE, type = "CQR",
#'                       quantiles = c(0.05, 0.95), outfun = "quantRF",  useCV = FALSE)
#' CIfun(Xtest)
#'
#' # naive method
#' CIfun <- conformalIte(X, Y, T, alpha = 0.1, algo = "naive", type = "CQR",
#'                       quantiles = c(0.05, 0.95), outfun = "quantRF",  useCV = FALSE)
#' CIfun(Xtest)
#'
#' # counterfactual method, Y and T needs to be observed
#' pstest <- pnorm(Xtest[, 1])
#' Ttest <- as.numeric(pstest < runif(ntest))
#' Y1test <- Xtest %*% beta + rnorm(ntest)
#' Y0test <- rnorm(ntest)
#' Ytest <- ifelse(Ttest == 1, Y1test, Y0test)
#' CIfun <- conformalIte(X, Y, T, alpha = 0.1, algo = "counterfactual", type = "CQR",
#'                       quantiles = c(0.05, 0.95), outfun = "quantRF",  useCV = FALSE)
#' CIfun(Xtest, Ytest, Ttest)
#' }
#'
#' @export
conformalIte <- function(X, Y, T,
                         alpha = 0.1,
                         algo = c("nest",
                                  "naive",
                                  "counterfactual"),
                         exact = FALSE,
                         type = c("CQR", "mean"),
                         side = c("two", "above", "below"),
                         quantiles = NULL,
                         outfun = NULL,
                         outparams = list(),
                         psfun = NULL,
                         psparams = list(),
                         cfprop = 0.5,
                         citype = c("CQR", "mean"),
                         lofun = NULL,
                         loquantile = 0.5,
                         loparams = list(),
                         upfun = NULL,
                         upquantile = 0.5,
                         upparams = list(),
                         useCV = FALSE,
                         trainprop = 0.75,
                         nfolds = 10,
                         wthigh = 20, wtlow = 0.05){
    algo <- algo[1]
    stopifnot(algo %in% c("nest", "naive", "counterfactual"))
    if (algo == "nest"){
        obj <- conformalIteNest(X, Y, T,
                                alpha,
                                type, side,
                                quantiles,
                                outfun, outparams,
                                psfun, psparams,
                                exact,
                                cfprop, citype,
                                lofun, loquantile, loparams,
                                upfun, upquantile, upparams,
                                useCV,
                                trainprop, nfolds,
                                wthigh, wtlow)
        function(X){
            predict(obj, X)
        }
    } else if (algo == "naive"){
        obj <- conformalIteNaive(X, Y, T,
                                 type, side,
                                 quantiles,
                                 outfun, outparams,
                                 psfun, psparams,
                                 useCV,
                                 trainprop, nfolds)
        function(X){
            predict(obj, X, alpha, wthigh, wtlow)$Ite
        }
    } else if (algo == "counterfactual"){
        obj <- conformalIteCf(X, Y, T,
                              type, side,
                              quantiles,
                              outfun, outparams,
                              psfun, psparams,
                              useCV,
                              trainprop, nfolds)
        function(X, Y, T){
            predict(obj, X, Y, T, alpha, wthigh, wtlow)
        }
    }
}
