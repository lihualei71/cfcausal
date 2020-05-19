conformalClassCV <- function(X, Y,
                             type = c("weighted", "unweighted"),
                             clsp = FALSE,
                             outfun = NULL,
                             outparams = list(),
                             wtfun = NULL,
                             nfolds = 10,
                             idlist = NULL){
    n <- length(Y)
    if (class(Y) != "factor"){
        stop("Y must be a factor")
    }
    
    ncl <- nlevels(Y)    
    type <- type[1]
    stopifnot(type %in% c("weighted", "unweighted"))

    if (is.null(idlist)){
        ## idlist <- gen_cv_ids(n, nfolds)
        tmp <- lapply(levels(Y), function(y){
            inds <- which(Y == y)
            ny <- length(inds)
            ids <- gen_cv_ids(ny, nfolds, offset = 0)
            lapply(ids, function(id){
                inds[id]
            })
        })
        idlist <- list()
        for (k in 1:nfolds){
            ids <- lapply(tmp, function(x){
                x[[k]]
            })
            idlist[[k]] <- unlist(ids)
        }
    }
    if (!is.list(idlist) || length(idlist) != nfolds){
        stop("idlist needs to a list of length 'nfolds'")
    }

    if (is.null(outfun)){
        outfun <- Boosting
    }

    wtfun0 <- NULL
    if (is.null(wtfun)){
        wtfun <- lapply(1:nfolds, function(k){
            function(X){
                rep(1, nrow(X))
            }
        })
    } else if (is.function(wtfun)){
        wtfun0 <- wtfun
        wtfun <- rep(list(wtfun), nfolds)
    } else if (!is.list(wtfun) || length(wtfun) != nfolds){
        stop("wtfun must be a function or a list (of functions) of length nfolds")
    }
    if (is.null(wtfun0)){
        wtfun0 <- wtfun
    }

    outparams0 <- outparams
    info <- list()
    for (k in 1:nfolds){
        testid <- idlist[[k]]
        Xtrain <- X[-testid, ]
        Ytrain <- Y[-testid]
        Xval <- X[testid, ]
        Yval <- Y[testid]

        outparams <- c(list(Y = Ytrain, X = Xtrain), outparams0)
        prmodel <- function(X){
            do.call(outfun, c(outparams, list(Xtest = X)))
        }
        phat <- prmodel(Xval)
        if (ncl == 2 && !is.matrix(phat)){
            phat <- cbind(1 - phat, phat)
            prmodel <- function(X){
                res <- do.call(outfun, c(outparams, list(Xtest = X)))
                cbind(1 - res, res)
            }
        }
        if (type == "weighted"){
            mpr <- freq(Ytrain)
        } else {
            mpr <- NULL
        }
        prscore <- conformalScoreClass(Yval, phat, type, mpr)
        wt <- wtfun[[k]](Xval)        
        if (clsp){
            prscore <- lapply(levels(Y), function(lv){
                prscore[Yval == lv]
            })
            wt <- lapply(levels(Y), function(lv){
                wt[Yval == lv]
            })
        }

        obj <- list(prscore = prscore,
                    phat = phat,
                    mpr = mpr,
                    wt = wt,
                    prmodel = prmodel)
        info[[k]] <- obj
    }

    res <- list(info = info,
                wtfun = wtfun0,
                type = type,
                clsp = clsp,
                labels = levels(Y),
                nfolds = nfolds,
                idlist = idlist)
    class(res) <- "conformalClassCV"
    return(res)
}

#' @export
predict.conformalClassCV <- function(object, Xtest,
                                     alpha = 0.1,
                                     wthigh = 20, wtlow = 0.05,
                                     ...){
    type <- object$type    
    clsp <- object$clsp
    nclass <- length(object$labels)
    nfolds <- object$nfolds
    info <- object$info

    if (clsp && length(alpha) > 1 && length(alpha) < nclass){
        stop("alpha must be a scalar or a vector with length = the number of classes when clsp = TRUE.")
    }
    if (clsp && length(alpha) == 1){
        alpha <- rep(alpha, nclass)
    }
    
    for (k in 1:nfolds){
        phat_test <- info[[k]]$prmodel(Xtest)
        if (type == "weighted"){
            mpr <- info[[k]]$mpr
            phat_test <- row_quo(phat_test, mpr)
        }
        info[[k]]$phat_test <- phat_test
    }

    if (is.function(object$wtfun)){
        wt_test <- object$wtfun(Xtest)
    } else {
        wt_test <- sapply(object$wtfun, function(wtfun){
            wtfun(Xtest)
        })
        wt_test <- rowMeans(wt_test)
    }
    
    if (!clsp){
        wt <- do.call(c, lapply(info, function(x){x$wt}))
        avg_wt <- mean(c(wt, wt_test))
        wt <- censoring(wt / avg_wt, wthigh, wtlow)
        wt_test <- censoring(wt_test / avg_wt, wthigh, wtlow)
        totw <- sum(wt)
        wt <- wt / totw    
        qt <- (1 + wt_test / totw) * (1 - alpha)
        qt <- pmin(qt, 1)
        CI <- sapply(1:nclass, function(j){
            sapply(1:length(qt), function(i){
                phat <- lapply(info, function(x){
                    x$phat_test[i, j] > x$prscore
                })
                phat <- unlist(phat)
                sum(phat * wt) < qt[i]
            })
        })
    } else {
        wt_test0 <- wt_test
        CI <- sapply(1:nclass, function(j){
            wt <- do.call(c, lapply(info, function(x){x$wt[[j]]}))
            avg_wt <- mean(c(wt, wt_test0))
            wt <- censoring(wt / avg_wt, wthigh, wtlow)
            wt_test <- censoring(wt_test0 / avg_wt, wthigh, wtlow)
            totw <- sum(wt)
            wt <- wt / totw    
            qt <- (1 + wt_test / totw) * (1 - alpha)
            qt <- pmin(qt, 1)
            sapply(1:length(qt), function(i){
                phat <- lapply(info, function(x){
                    x$phat_test[i, j] > x$prscore[[j]]
                })
                phat <- unlist(phat)
                sum(phat * wt) < qt[i]
            })
        })
    }

    res <- as.data.frame(CI)
    names(res) <- object$labels
    return(res)
}
