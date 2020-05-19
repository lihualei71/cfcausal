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
    ## Check the format
    type <- type[1]
    stopifnot(type %in% c("CQR", "mean"))
    side <- side[1]
    stopifnot(side %in% c("two", "above", "below"))
    citype <- citype[1]
    stopifnot(citype %in% c("CQR", "mean"))
    ntry <- 10

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

    ## ## Sign conformal inference
    ## tausign <- rep(0, nrow(CI_cf))
    ## tausign[CI_cf[, 1] > 0] <- 1
    ## tausign[CI_cf[, 2] < 0] <- -1
    ## tausign <- factor(tausign,
    ##                   levels = c(-1, 0, 1),
    ##                   labels = c("neg", "zero", "pos"))
    ## success <- FALSE
    ## for (i in 1:ntry){
    ##     ## use ntry for potentially randomized learners
    ##     obj <- try(conformalClass(Xcf, tausign, Xtest),
    ##                silent = TRUE)
    ##     if (class(obj) != "try-error"){
    ##         success <- TRUE
    ##         break
    ##     }
    ## }
    ## if (success){
    ##     signfun <- function(X){
    ##         predclass <- predict(obj, alpha = alpha,
    ##                              wthigh = wthigh, wtlow = wtlow)
    ##         predclass <- as.data.frame(predclass)
    ##         names(predclass) <- c("neg", "zero", "pos")
    ##         predsign <- rep(0, nrow(X))
    ##         predsign[predclass$pos & !predclass$neg & !predclass$zero] <- 1
    ##         predsign[predclass$neg & !predclass$pos & !predclass$zero] <- -1
    ##     }
    ## } else {
    ##     ## Report all zeros if it fails for ntry times
    ##     ## It is likely due to small number of +1/-1
    ##     ## in tausign
    ##     signfun <- function(X){
    ##         predsign <- rep(0, nrow(X))
    ##     }
    ##     warning("Sign conformal inference fails! Check the number of positives/negatives detected in CI_cf")
    ## }

    ## Get ITE intervals
    if (exact){
        ## Exact intervals
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
    } else {
        ## Inexact intervals 
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
            cbind(CI_lo, CI_up)
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
