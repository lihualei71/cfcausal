## Split conformal inference for counterfactuals. See ?conformalCf
conformalCf_split <- function(X, Y,
                              estimand,
                              type, side,
                              quantiles,
                              outfun, outparams,
                              psfun, psparams,
                              trainprop){
    T <- as.numeric(!is.na(Y))
    inds1 <- which(T == 1)
    inds0 <- which(T == 0)
    n1 <- length(inds1)
    n0 <- length(inds0)
    trainid1 <- sample(n1, floor(n1 * trainprop))
    trainid0 <- sample(n0, floor(n0 * trainprop))
    trainid <- c(inds1[trainid1], inds0[trainid0])
    Xtrain <- X[trainid, ,drop=FALSE]
    Ttrain <- T[trainid]

    psparams0 <- psparams
    if (estimand == "unconditional"){
        psparams <- c(list(Y = Ttrain, X = Xtrain), psparams0)
        wtfun <- function(X){
            ps <- do.call(psfun, c(list(Xtest = X), psparams))
            1 / ps
        }
    } else if (estimand == "nonmissing"){
        wtfun <- function(X){
            rep(1, nrow(X))
        }
    } else if (estimand == "missing"){
        psparams <- c(list(Y = Ttrain, X = Xtrain), psparams0)
        wtfun <- function(X){
            ps <- do.call(psfun, c(list(Xtest = X), psparams))
            (1 - ps) / ps
        }
    }

    X <- X[inds1, ,drop=FALSE]
    Y <- Y[inds1]
    res <- conformalSplit(X, Y,
                          type, side,
                          quantiles,
                          outfun, outparams,
                          wtfun,
                          trainprop, trainid1)
    return(res)
}
