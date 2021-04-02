## Generic weighted interval conformal inference. See ?conformalInt
conformalIntSplit <- function(X, Y,
                              type = c("CQR", "mean"),
                              lofun = NULL,
                              loquantile = 0.5,
                              loparams = list(),
                              upfun = NULL,
                              upquantile = 0.5,
                              upparams = list(),
                              wtfun = NULL,
                              trainprop = 0.75,
                              trainid = NULL){
    n <- nrow(Y)

    if (is.null(wtfun)){
        wtfun <- function(Xtest){
            rep(1, nrow(Xtest))
        }
    }

    if (type == "CQR"){
        if (is.null(loquantile) || is.null(upquantile)){
            stop("loquantile and upquantile should be provided if CQR is used.")
        }
    }

    if (is.null(trainid)){
        trainid <- sample(n, floor(n * trainprop))
    }
    Xtrain <- X[trainid, ,drop=FALSE]
    Ytrain <- Y[trainid, ]
    Xval <- X[-trainid, ,drop=FALSE]
    Yval <- Y[-trainid, ]

    if (type == "CQR"){
        loparams <- c(list(Y = Ytrain[, 1], X = Xtrain, quantiles = loquantile), loparams)
        upparams <- c(list(Y = Ytrain[, 2], X = Xtrain, quantiles = upquantile), upparams)
    } else if (type == "mean"){
        loparams <- c(list(Y = Ytrain[, 1], X = Xtrain), loparams)
        upparams <- c(list(Y = Ytrain[, 2], X = Xtrain), upparams)
    }

    Ylo_model <- function(X){
        do.call(lofun, c(loparams, list(Xtest = X)))
    }
    Yup_model <- function(X){
        do.call(upfun, c(upparams, list(Xtest = X)))
    }
    Ymodel <- list(Ylo_model, Yup_model)
    Yhat <- cbind(Ylo_model(Xval), Yup_model(Xval))
    Yscore <- conformalScore(Yval, Yhat, type, "two")
    wt <- wtfun(Xval)

    obj <- list(Yscore = Yscore, wt = wt,
                Ymodel = Ymodel, wtfun = wtfun,
                type = type,
                loquantile = loquantile,
                upquantile = upquantile,
                trainprop = trainprop,
                trainid = trainid)
    class(obj) <- "conformalIntSplit"
    return(obj)
}

#' Predict Method for conformalIntSplit objects
#'
#' Obtains predictive intervals on a testing dataset based on a \code{conformalIntSplit} object
#' from \code{\link{conformalInt}} with \code{useCV = FALSE}.
#'
#' Given a testing set \eqn{X_1, X_2, \ldots, X_n} and a weight function \eqn{w(x)}, the
#' weight of the weighted distribution \eqn{p_j = w(X_j) / (w(X_1) + \cdots + w(X_n))}.
#' In cases where some of \eqn{p_j} are extreme, we truncate \eqn{p_j} at level \code{wthigh}
#' and \code{wtlow} to ensure stability. If \code{wthigh = Inf, wtlow = 0}, no truncation
#' is being used.
#'
#' @param object an object of class \code{conformalIntSplit}; see \code{\link{conformalInt}}.
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
#' \code{\link{predict.conformalIntCV}}, \code{\link{conformalInt}}.
#'
#' @export
predict.conformalIntSplit <- function(object, Xtest,
                                      alpha = 0.1,
                                      wthigh = 20, wtlow = 0.05,
                                      ...){
    type <- object$type
    Yhat_test <- cbind(object$Ymodel[[1]](Xtest),
                       object$Ymodel[[2]](Xtest))
    wt_test <- object$wtfun(Xtest)

    avg_wt <- mean(c(object$wt, wt_test))
    wt <- censoring(object$wt / avg_wt, wthigh, wtlow)
    wt_test <- censoring(wt_test / avg_wt, wthigh, wtlow)

    totw <- sum(wt)
    wt <- wt / totw
    qt <- (1 + wt_test / totw) * (1 - alpha)
    qt <- pmin(qt, 1)
    Yslack <- weightedConformalCutoff(object$Yscore, wt, qt)

    Ylo <- Yhat_test[, 1] - Yslack
    Yup <- Yhat_test[, 2] + Yslack

    res <- data.frame(lower = Ylo, upper = Yup)
    return(res)
}
