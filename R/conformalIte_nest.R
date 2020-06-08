## Nested methods of Conformal inference for individual treatment effects for subjects with both
## missing potential outcome. See ?conformalIte
conformalIteNest <- function(X, Y, T,
                             alpha = 0.1,
                             type = c("CQR", "mean"),
                             side = c("two", "above", "below"),
                             quantiles = NULL,
                             outfun = NULL,
                             outparams = list(),
                             psfun = NULL,
                             psparams = list(),
                             exact = FALSE,
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
    ## Set default values for functions
    if (is.null(outfun)){
        outfun <- switch(type,
                         CQR = quantRF,
                         mean = RF)
    }
    if (is.null(psfun)){
        psfun <- Boosting
    }
    if (is.null(lofun)){
        lofun <- switch(citype,
                        CQR = quantRF,
                        mean = RF)
    }
    if (is.null(upfun)){
        upfun <- switch(citype,
                        CQR = quantRF,
                        mean = RF)
    }

    ## Reset effective alpha based on exact
    if (exact){
        alpha <- alpha / 2
    }

    ## Data splitting
    n <- length(Y)
    tr_id <- which(T == 1)
    ct_id <- which(T == 0)
    ntr <- length(tr_id)
    nct <- length(ct_id)
    ntr_cf <- ceiling(ntr * cfprop)
    nct_cf <- ceiling(nct * cfprop)
    tr_cfid <- sample(ntr, ntr_cf)
    ct_cfid <- sample(nct, nct_cf)
    cfid <- c(tr_id[tr_cfid], ct_id[ct_cfid])

    Xtrain <- X[-cfid, ]
    Ytrain <- Y[-cfid]
    Ttrain <- T[-cfid]
    Xcf <- X[cfid, ]
    Ycf <- Y[cfid]
    Tcf <- T[cfid]

    ## Get the counterfactual intervals and transform
    ## them into ITE intervals on the second fold
    obj <- conformalIteCf(Xtrain, Ytrain, Ttrain,
                          type, side,
                          quantiles,
                          outfun, outparams,
                          psfun, psparams,
                          useCV,
                          trainprop,
                          nfolds)
    CI_cf <- predict(obj, Xcf, Ycf, Tcf,
                     alpha, wthigh, wtlow)

    ## Get ITE intervals
    if (exact && (side == "two")){
        ## Exact two-sided intervals
        CIfun <- function(X){
            res <- conformalInt(Xcf, as.matrix(CI_cf),
                                citype,
                                lofun, loquantile, loparams,
                                upfun, upquantile, upparams,
                                NULL,
                                useCV,
                                trainprop, NULL,
                                nfolds, NULL)
            predict(res, X, alpha, wthigh, wtlow)
        }
    } else if (exact && (side == "above")){
        ## Exact right-sided intervals
        CIfun <- function(X){
            res <- conformal(Xcf, as.numeric(CI_cf[, 2]),
                             citype, "above", upquantile,
                             upfun, upparams,
                             NULL,
                             useCV,
                             trainprop, NULL,
                             nfolds, NULL)
            upper <- predict(res, X, alpha, wthigh, wtlow)[, 2]
            data.frame(lower = -Inf, upper = upper)
        }
    } else if (exact && (side == "below")){
        ## Exact left-sided intervals
        CIfun <- function(X){
            res <- conformal(Xcf, as.numeric(CI_cf[, 1]),
                             citype, "below", loquantile,
                             lofun, loparams,
                             NULL,
                             useCV,
                             trainprop, NULL,
                             nfolds, NULL)
            lower <- predict(res, X, alpha, wthigh, wtlow)[, 1]
            data.frame(lower = lower, upper = Inf)
        }
    } else if (!exact && (side == "two")){
        ## Inexact two-sided intervals
        if (citype == "CQR"){
            loparams <- c(list(Y = CI_cf[, 1], X = Xcf, quantiles = loquantile), loparams)
            upparams <- c(list(Y = CI_cf[, 2], X = Xcf, quantiles = upquantile), upparams)
        } else if (citype == "mean"){
            loparams <- c(list(Y = CI_cf[, 1], X = Xcf), loparams)
            upparams <- c(list(Y = CI_cf[, 2], X = Xcf), upparams)
        }
        CIfun <- function(X){
            CI_lo <- do.call(lofun, c(loparams, list(Xtest = X)))
            CI_up <- do.call(upfun, c(upparams, list(Xtest = X)))
            CI <- data.frame(lower = CI_lo,
                             upper = CI_up)
            return(CI)
        }
    } else if (!exact && (side == "above")){
        ## Inexact right-sided intervals
        if (citype == "CQR"){
            upparams <- c(list(Y = CI_cf[, 2], X = Xcf, quantiles = upquantile), upparams)
        } else if (citype == "mean"){
            upparams <- c(list(Y = CI_cf[, 2], X = Xcf), upparams)
        }
        CIfun <- function(X){
            CI_up <- do.call(upfun, c(upparams, list(Xtest = X)))
            CI <- data.frame(lower = -Inf,
                             upper = CI_up)
            return(CI)
        }
    } else if (!exact && (side == "below")){
        ## Inexact left-sided intervals
        if (citype == "CQR"){
            loparams <- c(list(Y = CI_cf[, 1], X = Xcf, quantiles = loquantile), loparams)
        } else if (citype == "mean"){
            loparams <- c(list(Y = CI_cf[, 1], X = Xcf), loparams)
        }
        CIfun <- function(X){
            CI_lo <- do.call(lofun, c(loparams, list(Xtest = X)))
            CI <- data.frame(lower = CI_lo,
                             upper = Inf)
            return(CI)
        }
    } 

    res <- list(CIfun = CIfun,
                CI_cf = CI_cf, cfid = cfid)
    class(res) <- "conformalIteNest"
    return(res)
}

predict.conformalIteNest <- function(object, X){
    object$CIfun(X)
}
