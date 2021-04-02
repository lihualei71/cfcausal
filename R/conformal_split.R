## Generic weighted split conformal inference. See ?conformal
conformalSplit <- function(X, Y,
                           type, side,
                           quantiles,
                           outfun, outparams,
                           wtfun,
                           trainprop, trainid){
    if (is.null(wtfun)){
        wtfun <- function(X){
            rep(1, nrow(X))
        }
    } else if (is.function(wtfun)){
        check_wtfun(wtfun)
    }

    if (type == "CQR"){
        if (is.null(quantiles)){
            stop("Quantiles should be provided if CQR is used.")
        } else if (length(quantiles) > 2){
            warning("At most two quantiles can be provided. Use the first two by default")
            quantiles <- quantiles[1:2]
        }
        if (side %in% c("above", "below") && length(quantiles) > 1){
            warning("At most one quantile can be provided when side = \"above\" or \"below\". Use the first one by default")
            quantiles <- quantiles[1]
        }
        outparams <- c(outparams, list(quantiles = quantiles))
    }

    n <- length(Y)
    if (is.null(trainid)){
        trainid <- sample(n, floor(n * trainprop))
    }
    Xtrain <- X[trainid, ,drop=FALSE]
    Ytrain <- Y[trainid]
    Xval <- X[-trainid, ,drop=FALSE]
    Yval <- Y[-trainid]

    outparams <- c(list(Y = Ytrain, X = Xtrain), outparams)

    Ymodel <- function(X){
        do.call(outfun, c(outparams, list(Xtest = X)))
    }
    Yhat <- Ymodel(Xval)
    Yscore <- conformalScore(Yval, Yhat, type, side)
    wt <- wtfun(Xval)

    obj <- list(Yscore = Yscore, wt = wt,
                Ymodel = Ymodel, wtfun = wtfun,
                type = type,
                side = side,
                quantiles = quantiles,
                trainprop = trainprop,
                trainid = trainid)
    class(obj) <- "conformalSplit"
    return(obj)
}

#' Predict Method for conformalSplit objects
#'
#' Obtains predictive intervals on a testing dataset based on a \code{conformalSplit} object
#' from \code{\link{conformal}} with \code{useCV = FALSE}.
#'
#' Given a testing set \eqn{X_1, X_2, \ldots, X_n} and a weight function \eqn{w(x)}, the
#' weight of the weighted distribution \eqn{p_j = w(X_j) / (w(X_1) + \cdots + w(X_n))}.
#' In cases where some of \eqn{p_j} are extreme, we truncate \eqn{p_j} at level \code{wthigh}
#' and \code{wtlow} to ensure stability. If \code{wthigh = Inf, wtlow = 0}, no truncation
#' is being used.
#'
#' @param object an object of class \code{conformalSplit}; see \code{\link{conformal}}.
#' @param Xtest testing covariates.
#' @param alpha confidence level.
#' @param wthigh upper truncation level of weights; see Details.
#' @param wtlow lower truncation level of weights; see Details.
#' @param ... other arguments
#'
#' @return predictive intervals. A data.frame with \code{nrow(Xtest)} rows and two columns:
#' "lower" for the lower bound and "upper" for the upper bound.
#'
#' @seealso
#' \code{\link{predict.conformalCV}}, \code{\link{conformal}}.
#'
#' @export
predict.conformalSplit <- function(object, Xtest,
                                   alpha = 0.1,
                                   wthigh = 20, wtlow = 0.05,
                                   ...){
    type <- object$type
    side <- object$side
    Yhat_test <- object$Ymodel(Xtest)
    wt_test <- object$wtfun(Xtest)

    avg_wt <- mean(c(object$wt, wt_test))
    wt <- censoring(object$wt / avg_wt, wthigh, wtlow)
    wt_test <- censoring(wt_test / avg_wt, wthigh, wtlow)

    totw <- sum(wt)
    wt <- wt / totw
    qt <- (1 + wt_test / totw) * (1 - alpha)
    qt <- pmin(qt, 1)
    Yslack <- weightedConformalCutoff(object$Yscore, wt, qt)

    if (type == "CQR" && side == "two"){
        Ylo <- Yhat_test[, 1] - Yslack
        Yup <- Yhat_test[, 2] + Yslack
    } else if (type == "mean" && side == "two"){
        Ylo <- Yhat_test - Yslack
        Yup <- Yhat_test + Yslack
    } else if (side == "above"){
        Ylo <- -Inf
        Yup <- Yhat_test + Yslack
    } else if (side == "below"){
        Ylo <- Yhat_test - Yslack
        Yup <- Inf
    }

    res <- data.frame(lower = Ylo, upper = Yup)
    return(res)
}
