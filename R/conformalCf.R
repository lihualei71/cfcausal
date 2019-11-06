#' @export
conformalCf <- function(X, Y, 
                        Xtest = NULL,
                        type = c("CQR", "mean"),
                        side = c("two", "above", "below"),
                        alpha = 0.05,
                        outfun = NULL,
                        outparams = list(),
                        psfun = NULL,
                        psparams = list(),
                        trainprop = 0.75,
                        trainid = NULL){
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
    
    T <- !is.na(Y)
    psparams <- c(list(Y = as.numeric(T), X = X), psparams)
    wtfun <- function(Xtest){
        1 / do.call(psfun, c(psparams, list(Xtest = Xtest)))
    }
    X <- X[T, ]
    Y <- Y[T]

    conformal(X, Y, Xtest, type, side, alpha,
              outfun, outparams, wtfun,
              trainprop, trainid)
}
