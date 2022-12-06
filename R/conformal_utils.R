## weightedConformalCutoff <- function(score, weight, qt){
##     quantile(score, qt)
## }

find_inds <- function(a, b){
    n <- length(a)
    b <- b - 1e-12
    ## n + 1 - rank(-c(a, b), ties.method = "first")[-(1:n)] + rank(-b, ties.method = "first")
    rank(c(a, b), ties.method = "first")[-(1:n)] - rank(b, ties.method = "first") + 1
}

## weightedConformalCutoff0 <- function(score, weight, qt){
##     ord <- order(score)
##     weight <- weight[ord]
##     score <- score[ord]
##     cw <- cumsum(weight)
##     cutoff <- sapply(qt, function(q){
##         ind <- min(which(cw >= q))
##         score[ind]
##     })
##     return(cutoff)
## }

weightedConformalCutoff <- function(score, weight, qt, useInf){
    ord <- order(score)
    weight <- weight[ord]
    score <- score[ord]
    cw <- cumsum(weight)
    inds <- find_inds(cw, pmin(qt, 1))
    cutoff <- score[inds]
    if (useInf){
        cutoff[qt >= 1] <- Inf
    }
    return(cutoff)
}

conformalScore <- function(Y, Yhat, type, side){
    if (is.vector(Y) || (is.matrix(Y) && ncol(Y) == 1)){
        if (type == "CQR" && side == "two"){
            score <- pmax(Yhat[, 1] - Y, Y - Yhat[, 2])
        } else if (type == "mean" && side == "two"){
            score <- abs(Yhat - Y)
        } else if (side == "above"){
            score <- Y - Yhat
        } else if (side == "below"){
            score <- Yhat - Y
        }
    } else if (is.matrix(Y) && ncol(Y) == 2){
        score <- pmax(Yhat[, 1] - Y[, 1], Y[, 2] - Yhat[, 2])
    }
    return(score)
}

conformalScoreClass <- function(Y, phat, type, wt){
    ncl <- nlevels(Y)
    if (type == "weighted"){
        phat <- row_quo(phat, wt)
    }
    Yid <- as.numeric(Y)
    score <- phat[cbind(1:length(Y), Yid)]
    return(score)
}

