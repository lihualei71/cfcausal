conformalIteCf <- function(X, Y, T,
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
    
    if (is.null(outfun)){
        outfun <- switch(type,
                         CQR = quantRF,
                         mean = RF)
    }
    if (is.null(psfun)){
        psfun <- Boosting
    }
    
    obj <- conformalIteNaive(X, Y, T,
                             type, side,
                             quantiles,
                             outfun, outparams,
                             psfun, psparams,
                             useCV,
                             trainprop,
                             nfolds)

    CIfun <- function(X, Y, T,
                      alpha, wthigh, wtlow){
        res <- predict(obj, X, 2 * alpha, wthigh, wtlow)
        CI <- matrix(NA, nrow(X), 2)
        CI[T == 1, 1] <- Y[T == 1] - res$Y0[T == 1, 2]
        CI[T == 1, 2] <- Y[T == 1] - res$Y0[T == 1, 1]
        CI[T == 0, 1] <- res$Y1[T == 0, 1] - Y[T == 0]
        CI[T == 0, 2] <- res$Y1[T == 0, 2] - Y[T == 0]
        return(CI)
    }

    res <- list(CIfun = CIfun)
    class(res) <- "conformalIteCf"
    return(res)
}

predict.conformalIteCf <- function(object,
                                   Xtest, Ytest, Ttest,
                                   alpha = 0.1,
                                   wthigh = 20, wtlow = 0.05){
    object$CIfun(Xtest, Ytest, Ttest, alpha, wthigh, wtlow)
}
