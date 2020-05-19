conformalClassSplit <- function(X, Y,
                                type = c("weighted", "unweighted"),
                                clsp = FALSE,
                                outfun = NULL,
                                outparams = list(),
                                wtfun = NULL,
                                trainprop = 0.75,
                                trainid = NULL){
    n <- length(Y)
    if (class(Y) != "factor"){
        stop("Y must be a factor")
    }
    if (is.null(outfun)){
        outfun <- Boosting
    }
    if (is.null(wtfun)){
        wtfun <- function(X){
            rep(1, nrow(X))
        }
    }
    type <- type[1]
    stopifnot(type %in% c("weighted", "unweighted"))
    ncl <- nlevels(Y)

    if (is.null(trainid)){
        trainid <- sample(n, floor(n * trainprop))
    }
    Xtrain <- X[trainid, ]
    Ytrain <- Y[trainid]
    Xval <- X[-trainid, ]
    Yval <- Y[-trainid]

    prmodel <- function(X){
        do.call(outfun, c(outparams, list(Y = Ytrain, X = Xtrain, Xtest = X)))
    }
    phat <- prmodel(Xval)
    if (ncl == 2 && !is.matrix(phat)){
        phat <- cbind(1 - phat, phat)
        prmodel <- function(X){
            res <- do.call(outfun, c(outparams, list(Y = Ytrain, X = Xtrain, Xtest = X)))
            cbind(1 - res, res)
        }
    }
    if (type == "weighted"){
        mpr <- freq(Ytrain)
    } else {
        mpr <- NULL
    }

    prscore <- conformalScoreClass(Yval, phat, type, mpr)
    wt <- wtfun(Xval)    
    if (clsp){
        prscore <- lapply(levels(Y), function(lv){
            prscore[Yval == lv]
        })
        wt <- lapply(levels(Y), function(lv){
            wt[Yval == lv]
        })
    }

    obj <- list(prscore = prscore,
                wt = wt,
                phat = phat,
                prmodel = prmodel,
                wtfun = wtfun,
                type = type,
                clsp = clsp,
                mpr = mpr,
                labels = levels(Y))
    class(obj) <- "conformalClassSplit"
    return(obj)
}

#' @export
predict.conformalClassSplit <- function(object, Xtest,
                                        alpha = 0.1,
                                        wthigh = 20, wtlow = 0.05,
                                        ...){
    clsp <- object$clsp
    nclass <- length(object$labels)
    if (clsp && length(alpha) > 1 && length(alpha) < nclass){
        stop("alpha must be a scalar or a vector with length = the number of classes when clsp = TRUE.")
    }
    if (clsp && length(alpha) == 1){
        alpha <- rep(alpha, nclass)
    }
    phat_test <- object$prmodel(Xtest)    
    wt_test <- object$wtfun(Xtest)
    type <- object$type
    if (type == "weighted"){
        mpr <- object$mpr
        phat_test <- row_quo(phat_test, mpr)
    }

    if (!clsp){
        avg_wt <- mean(c(object$wt, wt_test))
        wt <- censoring(object$wt / avg_wt, wthigh, wtlow)
        wt_test <- censoring(wt_test / avg_wt, wthigh, wtlow)
        totw <- sum(wt)
        wt <- wt / totw
        qt <- (1 + wt_test / totw) * (1 - alpha)
        qt <- pmin(qt, 1)
        cutoff <- -weightedConformalCutoff(-object$prscore, wt, qt)
        cutoff <- matrix(rep(cutoff, nclass), ncol = nclass)
    } else {
        wt_test0 <- wt_test
        cutoff <- sapply(1:nclass, function(j){
            avg_wt <- mean(c(object$wt[[j]], wt_test))
            wt <- censoring(object$wt[[j]] / avg_wt, wthigh, wtlow)
            wt_test <- censoring(wt_test / avg_wt, wthigh, wtlow)
            totw <- sum(wt)
            weight <- wt / totw
            qt <- (1 + wt_test / totw) * (1 - alpha)
            qt <- pmin(qt, 1)
            -weightedConformalCutoff(-object$prscore[[j]], weight, qt)
        })
    }

    res <- phat_test >= cutoff
    return(res)
}
