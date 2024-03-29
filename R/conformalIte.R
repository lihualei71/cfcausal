#' Conformal inference for individual treatment effects
#'
#' \code{conformalIte} supports four algorithms: the nested approach with exact and inexact
#' calibration for cases with both potential outcomes missing, the naive approach for cases with both potential outcomes missing and the counterfactual
#' inference for cases with only one potential outcome missing. For each algorithm, it supports both
#' split conformal inference and CV+, including weighted Jackknife+ as a special case. For each type, it
#' supports both conformalized quantile regression (CQR) and standard conformal inference based on conditional mean regression.
#'
#' @details The algorithm to be used is controlled by \code{algo} and \code{exact}:
#' \itemize{
#' \item (Default) when \code{algo = "nest"} and \code{exact = FALSE}, the inexact nested approach is used. It
#' first splits the data into two folds, with the second fold including \code{cfprop} fraction of units. Then it applies
#' \code{conformalCf} on the first fold to compute counterfactual intervals on the second fold, which further yields
#' interval estimates of ITE \eqn{\hat{C}(X_i)}. Finally it fits \eqn{\hat{C}^{L}(X_i)} and \eqn{\hat{C}^{R}(X_i)} on \eqn{X_i}'s.
#' \item When \code{algo = "nest"} and \code{exact = TRUE}, the exact nested approach is used. It has the same steps as the inexact nested approach to produce
#' ITE intervals \eqn{\hat{C}(X_i)}'s on the second fold but then applies \code{\link{conformalInt}} to calibrate them.
#' \item When \code{algo = "naive"}, the naive approach is used. It applies \code{\link{conformalCf}} on the data and
#' produce counterfactual intervals for both Y(1) and Y(0). The ITE intervals are computed by contrasting two counterfactual intervals.
#' \item When \code{algo = "counterfactual"}, it handles the case where the treatment assignments and the observed outcome are
#' both available for each testing point. As with the naive approach, it applies \code{\link{conformalCf}} on the data and
#' produce counterfactual intervals for both Y(1) and Y(0). The ITE intervals are then computed by contrasting the observed outcome
#' and the interval for the missing potential outcome.
#' }
#'
#' When \code{side = "above"},
#' intervals are of form [-Inf, a(x)] and when \code{side = "below"} the intervals are of form [a(x), Inf].
#'
#' When \code{type = "CQR"}, \code{quantiles} must be a vector of 2, regardless of \code{side}. When \code{side = "two"}, \code{quantiles} will be used in \code{outfun} for both Y(1) and Y(0); when \code{side = "above"} or \code{"below"}, \code{quantiles[1]} will be used for Y(0) and \code{quantiles[2]} will be used for Y(1).
#' 
#' \code{outfun} is applied to both Y(1) and Y(0). \code{outfun} can be a valid string, including
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
#' \code{lofun} and \code{upfun} have the same forms as \code{outfun} except that the input \code{quantiles}
#' must be scalar when \code{citype = "CQR"}, instead of a vector of 2, because only one conditional quantile
#' is fitted. The argument \code{loquantile} is used for \code{lofun} and the argument \code{hiquantile} is used
#' for \code{upfun}. Moreover, the output must be a vector giving the conditional quantile estimate or conditional mean
#' estimate. Other optional arguments can be passed into \code{lofun} through \code{loparams} and \code{upfun}
#' through \code{upparams}.
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
#' @param Y observed outcome vector.
#' @param T treatment indicator, a binary vector.
#' @param alpha confidence level.
#' @param algo a string that takes values in \{"nest", "naive", "counterfactual"\}. See Details.
#' @param exact a logical indicating whether the exact calibration is used for nested approach. Used only when \code{algo = "nest"}. See Details.
#' @param type a string that takes values in \{"CQR", "mean"\}.
#' @param side a string that takes values in \{"two", "above", "below"\}. See Details.
#' @param quantiles for covariates in the training data. Used only when \code{type = "CQR"}. See Details.
#' @param outfun a function that models the conditional mean or quantiles, or a valid string. 
#'               The default is random forest when \code{type = "mean"} and quantile random forest when
#'               \code{type = "CQR"}. See Details.
#' @param outparams a list of other parameters to be passed into \code{outfun}.
#' @param psfun a function that models the missing mechanism (probability of missing given X), or a valid string. 
#'               The default is "Boosting". See Details.
#' @param psparams a list of other parameters to be passed into \code{psfun}.
#' @param cfprop the proportion of units to be used to compute ITE intervals in nested approach. Used only when
#' \code{algo = "nest"}.
#' @param citype the type of interval conformal inference used in the nested approach with exact calibration.
#' Used only when \code{algo = "nest"} and \code{exact = TRUE}.
#' @param lofun the function to fit the lower bound, or a valid string. Used only when
#' \code{algo = "nest"}. See Details.
#' @param loquantile the quantile to fit for \code{lofun}; see Details. Used only when
#' \code{algo = "nest"} and \code{citype = "CQR"}. See Details.
#' @param loparams a list of other parameters to be passed into \code{lofun}.
#' @param upfun the function to fit the upper bound, or a valid string. Used only when
#' \code{algo = "nest"}. See Details.
#' @param upquantile the quantile to fit for \code{upfun}. Used only when
#' \code{algo = "nest"} and \code{citype = "CQR"}. See Details.
#' @param upparams a list of other parameters to be passed into \code{upfun}.
#' @param useCV FALSE for split conformal inference and TRUE for CV+.
#' @param trainprop proportion of units for training \code{outfun}. The default if 75\%. Used only when \code{useCV = FALSE}.
#' @param nfolds number of folds. The default is 10. Used only when \code{useCV = TRUE}. 
#' @param wthigh upper truncation level of weights. See \code{\link{predict.conformalSplit}} or \code{\link{predict.conformalCV}}.
#' @param wtlow lower truncation level of weights. See \code{\link{predict.conformalSplit}} or \code{\link{predict.conformalCV}}.
#' @param useInf if FALSE then replace infinity by the maximum conformity score.
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
                         loquantile = 0.4,
                         loparams = list(),
                         upfun = NULL,
                         upquantile = 0.6,
                         upparams = list(),
                         useCV = FALSE,
                         trainprop = 0.75,
                         nfolds = 10,
                         wthigh = 20, wtlow = 0.05,
                         useInf = FALSE){
    ## Check the format
    type <- type[1]
    stopifnot(type %in% c("CQR", "mean"))
    side <- side[1]
    stopifnot(side %in% c("two", "above", "below"))
    citype <- citype[1]
    stopifnot(citype %in% c("CQR", "mean"))
    algo <- algo[1]
    stopifnot(algo %in% c("nest", "naive", "counterfactual"))

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

    if (is.null(lofun)){
        lofun <- switch(citype,
                        CQR = quantRF,
                        mean = RF)
    } else if (is.character(lofun)){
        lofun <- str_outfun(lofun[1])
    } else if (is.function(lofun)){
        check_outfun(lofun, citype)
    } else {
        stop("lofun must be NULL or a string or a function")
    }

    if (is.null(upfun)){
        upfun <- switch(citype,
                        CQR = quantRF,
                        mean = RF)
    } else if (is.character(upfun)){
        upfun <- str_outfun(upfun[1])
    } else if (is.function(upfun)){
        check_outfun(upfun, citype)
    } else {
        stop("upfun must be NULL or a string or a function")
    }    
    
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
                                wthigh, wtlow,
                                useInf)
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
            predict(obj, X, alpha, wthigh, wtlow, useInf)$Ite
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
            predict(obj, X, Y, T, alpha, wthigh, wtlow, useInf)
        }
    }
}
