#' @export
conformal <- function(X, Y,
                      Xtest = NULL,
                      type = c("CQR", "mean"),
                      side = c("two", "above", "below"),
                      alpha = 0.05,
                      fun = NULL,
                      params = list(),
                      wtfun = NULL,
                      trainprop = 0.75,
                      trainid = NULL){
    n <- length(Y)
    type <- type[1]
    stopifnot(type %in% c("CQR", "mean"))
    side <- side[1]
    stopifnot(side %in% c("two", "above", "below"))
    if (is.null(fun)){
        fun <- switch(type,
                      CQR = quantRF,
                      mean = RF)
    }
    if (is.null(wtfun)){
        wtfun <- function(Xtest){
            rep(1, nrow(Xtest))
        }
    }

    if (is.null(trainid)){
        trainid <- sample(n, floor(n * trainprop))
    }
    Xtrain <- X[trainid, ]
    Ytrain <- Y[trainid]
    Xval <- X[-trainid, ]
    Yval <- Y[-trainid]

    params <- c(list(Y = Ytrain, X = Xtrain), params)
    if (type == "CQR" && side == "two"){
        params <- c(params, list(quantiles = c(alpha / 2, 1 - alpha / 2)))
    } else if (type == "CQR" && side == "above"){
        params <- c(params, list(quantiles = 1 - alpha))
    } else if (type == "CQR" && side == "below"){
        params <- c(params, list(quantiles = alpha))
    }

    Ymodel <- function(Xtest){
        do.call(fun, c(params, list(Xtest = Xtest)))
    }
    Yhat <- Ymodel(Xval)
    Yscore <- conformalScore(Yval, Yhat, type, side)
    wt <- wtfun(Xval)

    if (!is.null(Xtest)){
        Yhat_test <- Ymodel(Xtest)
        wt_test <- wtfun(Xtest)
    } else{
        Yhat_test <- NULL
        wt_test <- NULL
    }

    obj <- list(Yscore = Yscore, wt = wt,
                Ymodel = Ymodel, wtfun = wtfun,
                Yhat_test = Yhat_test,
                wt_test = wt_test,
                type = type, side = side)
    class(obj) <- "conformal"
    return(obj)
}

predict.conformal <- function(obj, Xtest = NULL,
                              alpha = 0.05,
                              wthigh = 20, wtlow = 0.05,
                              ...){
    type <- obj$type
    side <- obj$side
    if (!is.null(Xtest)){
        Yhat_test <- obj$Ymodel(Xtest)
        wt_test <- obj$wtfun(Xtest)
    } else {
        Yhat_test <- obj$Yhat_test
        wt_test <- obj$wt_test
    }

    wt <- censoring(obj$wt, wthigh, wtlow)
    wt_test <- censoring(wt_test, wthigh, wtlow)
    Yslack <- weightedConformal(obj$Yscore, wt, wt_test, alpha)
    if (type == "CQR" && side == "two"){
        Ylow <- Yhat_test[, 1] - Yslack
        Yhigh <- Yhat_test[, 2] + Yslack
    } else if (type == "mean" && side == "two"){
        Ylow <- Yhat_test - Yslack
        Yhigh <- Yhat_test + Yslack
    } else if (side == "above"){
        Ylow <- -Inf
        Yhigh <- Yhat_test + Yslack
    } else if (side == "below"){
        Ylow <- Yhat_test - Yslack
        Yhigh <- Inf
    }

    res <- data.frame(Ylow = Ylow, Yhigh = Yhigh)
    return(res)
}
