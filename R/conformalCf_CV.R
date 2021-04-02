## CV+ for counterfactuals. See ?conformalCf
conformalCf_CV <- function(X, Y,
                           estimand,
                           type, side,
                           quantiles,
                           outfun, outparams,
                           psfun, psparams,
                           nfolds){
    T <- as.numeric(!is.na(Y))
    inds1 <- which(T == 1)
    inds0 <- which(T == 0)
    n1 <- length(inds1)
    n0 <- length(inds0)
    if (n1 < nfolds){
        stop("Insufficient non-missing data")
    }
    idlist1 <- gen_cv_ids(n1, nfolds, offset = 0)
    idlist0 <- gen_cv_ids(n0, nfolds, offset = 0)
    idlist <- lapply(1:nfolds, function(k){
        c(inds1[idlist1[[k]]], inds0[idlist0[[k]]])
    })

    psparams0 <- psparams
    if (estimand == "unconditional"){
        wtfun <- lapply(1:nfolds, function(k){
            testid <- idlist[[k]]
            Xtrain <- X[-testid, ,drop=FALSE]
            Ttrain <- T[-testid]
            psparams <- c(list(Y = Ttrain, X = Xtrain), psparams0)
            function(X){
                ps <- do.call(psfun, c(psparams, list(Xtest = X)))
                1 / ps
            }
        })
        psparams <- c(list(Y = T, X = X), psparams0)
        wtfun_test <- function(X){
            ps <- do.call(psfun, c(psparams, list(Xtest = X)))
            1 / ps
        }
    } else if (estimand == "nonmissing"){
        wtfun_test <- function(X){
            rep(1, nrow(X))
        }
        wtfun <- lapply(1:nfolds, function(k){
            wtfun_test
        })
    } else if (estimand == "missing"){
        wtfun <- lapply(1:nfolds, function(k){
            testid <- idlist[[k]]
            Xtrain <- X[-testid, ,drop=FALSE]
            Ttrain <- T[-testid]
            psparams <- c(list(Y = Ttrain, X = Xtrain), psparams0)
            function(X){
                ps <- do.call(psfun, c(psparams, list(Xtest = X)))
                (1 - ps) / ps
            }
        })
        psparams <- c(list(Y = T, X = X), psparams0)
        wtfun_test <- function(X){
            ps <- do.call(psfun, c(psparams, list(Xtest = X)))
            (1 - ps) / ps
        }
    }

    X <- X[inds1, ,drop=FALSE]
    Y <- Y[inds1]
    res <- conformalCV(X, Y,
                       type, side,
                       quantiles,
                       outfun, outparams,
                       wtfun,
                       nfolds, idlist1)
    res$wtfun <- wtfun_test
    return(res)
}
