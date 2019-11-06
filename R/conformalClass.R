#' @export
conformalClass <- function(X, Y,
                           Xtest = NULL,
                           alpha = 0.05,
                           type = c("weighted", "unweighted"),
                           clsp = FALSE,
                           fun = NULL,
                           params = list(),
                           wtfun = NULL,
                           trainprop = 0.75,
                           trainid = NULL){
    n <- length(Y)
    if (class(Y) != "factor"){
        stop("Y must be a factor")
    }
    if (is.null(fun)){
        fun <- RF
    }
    if (is.null(wtfun)){
        wtfun <- function(Xtest){
            rep(1, nrow(Xtest))
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

    prmodel <- function(Xtest){
        do.call(fun, c(params, list(Y = Ytrain, X = Xtrain, Xtest = Xtest)))
    }
    phat <- prmodel(Xval)
    if (type == "weighted"){
        mpr <- freq(Ytrain)
    } else {
        mpr <- NULL
    }
    prscore <- conformalScoreClass(Yval, phat, type, mpr)
    if (clsp){
        prscore <- lapply(levels(Y), function(lv){
            prscore[Yval == lv]
        })
    }
    wt <- wtfun(Xval)

    if (!is.null(Xtest)){
        phat_test <- prmodel(Xtest)
        wt_test <- wtfun(Xtest)
    } else{
        phat_test <- NULL
        wt_test <- NULL
    }

    obj <- list(prscore = prscore, wt = wt,
                phat = phat,
                prmodel = prmodel, wtfun = wtfun,
                phat_test = phat_test,
                wt_test = wt_test,
                type = type,
                clsp = clsp,
                mpr = mpr,
                labels = levels(Y))
    class(obj) <- "conformalClass"
    return(obj)
}

#' @export
predict.conformalClass <- function(object, Xtest = NULL,
                                   alpha = 0.05,
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
    if (!is.null(Xtest)){
        phat_test <- object$prmodel(Xtest)
        wt_test <- object$wtfun(Xtest)
    } else {
        phat_test <- object$phat_test
        wt_test <- object$wt_test
    }
    type <- object$type
    if (type == "weighted"){
        mpr <- object$mpr
        phat_test <- row_quo(phat_test, mpr)
    }

    wt <- censoring(object$wt, wthigh, wtlow)
    wt_test <- censoring(wt_test, wthigh, wtlow)
    if (!clsp){
        cutoff <- weightedConformal(object$prscore, wt, wt_test, 1 - alpha)
        cutoff <- matrix(rep(cutoff, nclass), ncol = nclass)
    } else {
        cutoff <- sapply(object$prscore, function(score){
            weightedConformal(score, wt, wt_test, 1 - alpha)
        })
    }

    res <- phat_test >= cutoff
    return(res)

}
