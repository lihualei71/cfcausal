## Naive methods of Conformal inference for individual treatment effects for subjects with both
## missing potential outcome. See ?conformalIte
conformalIteNaive <- function(X, Y, T,
                              type, side,
                              quantiles,
                              outfun, outparams,
                              psfun, psparams,
                              useCV,
                              trainprop,
                              nfolds){
    n <- length(Y)

    Y1 <- Y0 <- Y
    Y1[T == 0] <- NA
    Y0[T == 1] <- NA
    inds <- which(T == 1)
    Xtrain <- X

    estimand1 <- "missing"
    side1 <- switch(side,
                    two = "two",
                    above = "above",
                    below = "below")
    if (side == "two"){
        quantiles1 <- quantiles
    } else {
        quantiles1 <- quantiles[2]
    }
    obj1 <- conformalCf(Xtrain, Y1,
                        estimand1,
                        type, side1,
                        quantiles1,
                        outfun, outparams,
                        psfun, psparams,
                        useCV,
                        trainprop,
                        nfolds)
    Y1_CIfun <- function(X, alpha, wthigh, wtlow){
        predict(obj1, X, alpha = alpha / 2,
                wthigh = wthigh, wtlow = wtlow)
    }

    estimand0 <- "missing"
    side0 <- switch(side,
                    two = "two",
                    above = "below",
                    below = "above")
    if (side == "two"){
        quantiles0 <- quantiles
    } else {
        quantiles0 <- quantiles[1]
    }
    obj0 <- conformalCf(Xtrain, Y0,
                        estimand0,
                        type, side0,
                        quantiles0,
                        outfun, outparams,
                        psfun, psparams,
                        useCV,
                        trainprop,
                        nfolds)
    Y0_CIfun <- function(X, alpha, wthigh, wtlow){
        predict(obj0, X, alpha = alpha / 2,
                wthigh = wthigh, wtlow = wtlow)
    }

    Ite_CIfun <- function(X, alpha, wthigh, wtlow){
        Y1_CI <- Y1_CIfun(X, alpha, wthigh, wtlow)
        Y0_CI <- Y0_CIfun(X, alpha, wthigh, wtlow)
        CI <- data.frame(lower = Y1_CI[, 1] - Y0_CI[, 2],
                         upper = Y1_CI[, 2] - Y0_CI[, 1])
    }

    res <- list(Ite_CIfun = Ite_CIfun,
                Y1_CIfun = Y1_CIfun,
                Y0_CIfun = Y0_CIfun)
    class(res) <- "conformalIteNaive"
    return(res)
}

predict.conformalIteNaive <- function(object, Xtest,
                                      alpha = 0.1,
                                      wthigh = 20, wtlow = 0.05){
    Ite_CI <- object$Ite_CIfun(Xtest, alpha, wthigh, wtlow)
    Y1_CI <- object$Y1_CIfun(Xtest, alpha, wthigh, wtlow)
    Y0_CI <- object$Y0_CIfun(Xtest, alpha, wthigh, wtlow)
    list(Ite = Ite_CI, Y1 = Y1_CI, Y0 = Y0_CI)
}
