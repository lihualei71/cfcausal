## Generic weighted interval CV+. See ?conformalInt
conformalIntCV <- function(X, Y,
                           type = c("CQR", "mean"),
                           lofun = NULL,
                           loquantile = 0.5,
                           loparams = list(),
                           upfun = NULL,
                           upquantile = 0.5,
                           upparams = list(),
                           wtfun = NULL,
                           nfolds = 10,
                           idlist = NULL){
    n <- nrow(Y)

    if (is.null(idlist)){
        idlist <- gen_cv_ids(n, nfolds)
    }
    if (!is.list(idlist) || length(idlist) != nfolds){
        stop("idlist needs to a list of length 'nfolds'")
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

    if (type == "CQR"){
        if (is.null(loquantile) || is.null(upquantile)){
            stop("loquantile and upquantile should be provided if CQR is used.")
        }
        loparams <- c(list(quantiles = loquantile), loparams)
        upparams <- c(list(quantiles = upquantile), upparams)
    }

    loparams0 <- loparams
    upparams0 <- upparams
    info <- list()
    for (k in 1:nfolds){
        testid <- idlist[[k]]
        Xtrain <- X[-testid, ,drop=FALSE]
        Ytrain <- Y[-testid, ]
        Xval <- X[testid, ,drop=FALSE]
        Yval <- Y[testid, ]

        loparams <- c(list(Y = Ytrain[, 1], X = Xtrain), loparams0)
        upparams <- c(list(Y = Ytrain[, 2], X = Xtrain), upparams0)
        Ylo_model <- function(X){
            do.call(lofun, c(loparams, list(Xtest = X)))
        }
        Yup_model <- function(X){
            do.call(upfun, c(upparams, list(Xtest = X)))
        }
        Ymodel <- list(Ylo_model, Yup_model)
        Yhat <- cbind(Ylo_model(Xval), Yup_model(Xval))
        Yscore <- conformalScore(Yval, Yhat, type, "two")
        wt <- wtfun[[k]](Xval)

        obj <- list(Yscore = Yscore,
                    wt = wt,
                    Ymodel = Ymodel)
        info[[k]] <- obj
    }

    res <- list(info = info,
                wtfun = wtfun0,
                type = type,
                loquantile = loquantile,
                upquantile = upquantile,
                nfolds = nfolds,
                idlist = idlist)
    class(res) <- "conformalIntCV"
    return(res)
}

#' Predict Method for conformalIntCV objects
#'
#' Obtains predictive intervals on a testing dataset based on a \code{conformalIntCV} object
#' from \code{\link{conformalInt}} with \code{useCV = TRUE}.
#'
#' Given a testing set \eqn{X_1, X_2, \ldots, X_n} and a weight function \eqn{w(x)}, the
#' weight of the weighted distribution \eqn{p_j = w(X_j) / (w(X_1) + \cdots + w(X_n))}.
#' In cases where some of \eqn{p_j} are extreme, we truncate \eqn{p_j} at level \code{wthigh}
#' and \code{wtlow} to ensure stability. If \code{wthigh = Inf, wtlow = 0}, no truncation
#' is being used.
#'
#' @param object an object of class \code{conformalIntCV}; see \code{\link{conformalInt}}.
#' @param Xtest testing covariates.
#' @param alpha confidence level.
#' @param wthigh upper truncation level of weights; see Details.
#' @param wtlow lower truncation level of weights; see Details.
#' @param ... other arguments
#'
#' @return predictive intervals. A data.frame with \code{nrow(Xtest)} rows and two columns:
#' "lower" for the lower bound and "upper" for the upper bound.
#'
#' @seealso
#' \code{\link{predict.conformalIntSplit}}, \code{\link{conformalInt}}.
#'
#' @export
predict.conformalIntCV <- function(object, Xtest,
                                   alpha = 0.1,
                                   wthigh = 20, wtlow = 0.05,
                                   ...){
    type <- object$type
    nfolds <- object$nfolds
    info <- object$info


    for (k in 1:nfolds){
        info[[k]]$Yhat_test <- cbind(
            info[[k]]$Ymodel[[1]](Xtest),
            info[[k]]$Ymodel[[2]](Xtest))
    }

    wt <- do.call(c, lapply(info, function(x){x$wt}))
    if (is.function(object$wtfun)){
        wt_test <- object$wtfun(Xtest)
    } else {
        wt_test <- sapply(object$wtfun, function(wtfun){
            wtfun(Xtest)
        })
        wt_test <- rowMeans(wt_test)
    }
    avg_wt <- mean(c(wt, wt_test))
    wt <- censoring(wt / avg_wt, wthigh, wtlow)
    wt_test <- censoring(wt_test / avg_wt, wthigh, wtlow)

    totw <- sum(wt)
    wt <- wt / totw
    qt <- (1 + wt_test / totw) * (1 - alpha)
    qt <- pmin(qt, 1)

    CI <- sapply(1:length(qt), function(i){
        Ylo <- lapply(info, function(x){
            x$Yhat_test[i, 1] - x$Yscore
        })
        Ylo <- do.call(c, Ylo)
        Ylo <- -weightedConformalCutoff(-Ylo, wt, qt[i])
        Yup <- lapply(info, function(x){
            x$Yhat_test[i, 2] + x$Yscore
        })
        Yup <- do.call(c, Yup)
        Yup <- weightedConformalCutoff(Yup, wt, qt[i])
        c(Ylo, Yup)
    })

    res <- data.frame(lower = as.numeric(CI[1, ]),
                      upper = as.numeric(CI[2, ]))
    return(res)
}
