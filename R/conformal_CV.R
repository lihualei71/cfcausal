## Generic weighted CV+. See ?conformal
conformalCV <- function(X, Y,
                        type, side,
                        quantiles,
                        outfun, outparams,
                        wtfun,
                        nfolds, idlist){
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
    sapply(wtfun, check_wtfun)
    if (is.null(wtfun0)){
        wtfun0 <- wtfun
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

    n <- length(Y)
    if (is.null(idlist)){
        idlist <- gen_cv_ids(n, nfolds)
    }
    if (!is.list(idlist) || length(idlist) != nfolds){
        stop("idlist needs to a list of length 'nfolds'")
    }

    outparams0 <- outparams
    info <- list()
    for (k in 1:nfolds){
        testid <- idlist[[k]]
        Xtrain <- X[-testid, ,drop=FALSE]
        Ytrain <- Y[-testid]
        Xval <- X[testid, ,drop=FALSE]
        Yval <- Y[testid]

        outparams <- c(list(Y = Ytrain, X = Xtrain), outparams0)
        Ymodel <- function(X){
            do.call(outfun, c(outparams, list(Xtest = X)))
        }
        Yhat <- Ymodel(Xval)
        Yscore <- conformalScore(Yval, Yhat, type, side)
        wt <- wtfun[[k]](Xval)

        obj <- list(Yscore = Yscore,
                    wt = wt,
                    Ymodel = Ymodel)
        info[[k]] <- obj
    }

    res <- list(info = info,
                wtfun = wtfun0,
                type = type,
                side = side,
                quantiles = quantiles,
                nfolds = nfolds,
                idlist = idlist)
    class(res) <- "conformalCV"
    return(res)
}

#' Predict Method for conformalCV objects
#'
#' Obtains predictive intervals on a testing dataset based on a \code{conformalCV} object
#' from \code{\link{conformal}} with \code{useCV = TRUE}.
#'
#' Given a testing set \eqn{X_1, X_2, \ldots, X_n} and a weight function \eqn{w(x)}, the
#' weight of the weighted distribution \eqn{p_j = w(X_j) / (w(X_1) + \cdots + w(X_n))}.
#' In cases where some of \eqn{p_j} are extreme, we truncate \eqn{p_j} at level \code{wthigh}
#' and \code{wtlow} to ensure stability. If \code{wthigh = Inf, wtlow = 0}, no truncation
#' is being used.
#'
#' @param object an object of class \code{conformalCV}; see \code{\link{conformal}}.
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
#' \code{\link{predict.conformalSplit}}, \code{\link{conformal}}.
#'
#' @export
predict.conformalCV <- function(object, Xtest,
                                alpha = 0.1,
                                wthigh = 20, wtlow = 0.05,
                                ...){
    type <- object$type
    side <- object$side
    nfolds <- object$nfolds
    info <- object$info

    for (k in 1:nfolds){
        info[[k]]$Yhat_test <- info[[k]]$Ymodel(Xtest)
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

    if (type == "CQR" && side == "two"){
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
    } else if (type == "mean" && side == "two"){
        CI <- sapply(1:length(qt), function(i){
            Ylo <- lapply(info, function(x){
                x$Yhat_test[i] - x$Yscore
            })
            Ylo <- do.call(c, Ylo)
            Ylo <- -weightedConformalCutoff(-Ylo, wt, qt[i])
            Yup <- lapply(info, function(x){
                x$Yhat_test[i] + x$Yscore
            })
            Yup <- do.call(c, Yup)
            Yup <- weightedConformalCutoff(Yup, wt, qt[i])
            c(Ylo, Yup)
        })
    } else if (side == "above"){
        CI <- sapply(1:length(qt), function(i){
            Ylo <- -Inf
            Yup <- lapply(info, function(x){
                x$Yhat_test[i] + x$Yscore
            })
            Yup <- do.call(c, Yup)
            Yup <- weightedConformalCutoff(Yup, wt, qt[i])
            c(Ylo, Yup)
        })
    } else if (side == "below"){
        CI <- sapply(1:length(qt), function(i){
            Ylo <- lapply(info, function(x){
                x$Yhat_test[i] - x$Yscore
            })
            Ylo <- do.call(c, Ylo)
            Ylo <- -weightedConformalCutoff(-Ylo, wt, qt[i])
            Yup <- Inf
            c(Ylo, Yup)
        })
    }

    res <- data.frame(lower = as.numeric(CI[1, ]),
                      upper = as.numeric(CI[2, ]))
    return(res)
}
