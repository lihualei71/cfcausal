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
    if (!is.matrix(Y) || ncol(Y) != 2){
        stop("Y must a matrix with 2 columns")
    }
    n <- nrow(Y)
    type <- type[1]
    stopifnot(type %in% c("CQR", "mean"))

    if (is.null(lofun)){
        lofun <- switch(type,
                        CQR = quantRF,
                        mean = RF)
    }
    if (is.null(upfun)){
        upfun <- switch(type,
                        CQR = quantRF,
                        mean = RF)
    }
    
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
    Xtrain <- X[trainid, ]
    Ytrain <- Y[trainid, ]
    Xval <- X[-trainid, ]
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
                type = type)
    class(obj) <- "conformalIntSplit"
    return(obj)
}

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
