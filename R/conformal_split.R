conformalSplit <- function(X, Y,
                           type = c("CQR", "mean"),
                           side = c("two", "above", "below"),
                           quantiles = NULL,
                           outfun = NULL,
                           outparams = list(),
                           wtfun = NULL,
                           trainprop = 0.75,
                           trainid = NULL){
    n <- length(Y)
    type <- type[1]
    stopifnot(type %in% c("CQR", "mean"))
    side <- side[1]
    stopifnot(side %in% c("two", "above", "below"))

    if (is.null(outfun)){
        outfun <- switch(type,
                         CQR = quantRF,
                         mean = RF)
    }

    if (is.null(wtfun)){
        wtfun <- function(X){
            rep(1, nrow(X))
        }
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

    if (is.null(trainid)){
        trainid <- sample(n, floor(n * trainprop))
    }
    Xtrain <- X[trainid, ]
    Ytrain <- Y[trainid]
    Xval <- X[-trainid, ]
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
                type = type, side = side)
    class(obj) <- "conformalSplit"
    return(obj)
}

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
