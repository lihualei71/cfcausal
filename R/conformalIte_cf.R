## Conformal inference for individual treatment effects for subjects with only one
## missing potential outcome. See ?conformalIte
conformalIteCf <- function(X, Y, T,
                           type, side,
                           quantiles,
                           outfun, outparams,
                           psfun, psparams,
                           useCV,
                           trainprop,
                           nfolds){
    obj <- conformalIteNaive(X, Y, T,
                             type, side,
                             quantiles,
                             outfun, outparams,
                             psfun, psparams,
                             useCV,
                             trainprop,
                             nfolds)

    CIfun <- function(X, Y, T,
                      alpha, wthigh, wtlow, useInf){
        res <- predict(obj, X, 2 * alpha, wthigh, wtlow, useInf)
        CI <- matrix(NA, nrow(X), 2)
        CI[T == 1, 1] <- Y[T == 1] - res$Y0[T == 1, 2]
        CI[T == 1, 2] <- Y[T == 1] - res$Y0[T == 1, 1]
        CI[T == 0, 1] <- res$Y1[T == 0, 1] - Y[T == 0]
        CI[T == 0, 2] <- res$Y1[T == 0, 2] - Y[T == 0]
        CI <- as.data.frame(CI)
        names(CI) <- c("lower", "upper")
        return(CI)
    }

    res <- list(CIfun = CIfun)
    class(res) <- "conformalIteCf"
    return(res)
}

predict.conformalIteCf <- function(object,
                                   Xtest, Ytest, Ttest,
                                   alpha = 0.1,
                                   wthigh = 20, wtlow = 0.05,
                                   useInf = FALSE){
    object$CIfun(Xtest, Ytest, Ttest, alpha, wthigh, wtlow, useInf)
}
