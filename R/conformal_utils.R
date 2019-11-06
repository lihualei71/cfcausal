weightedConformal <- function(score, weight, nws, alpha){
    totw <- sum(weight)
    qt <- (totw + nws) / totw * (1 - alpha)
    qt <- pmin(qt, 1)
    cutoff <- quantile(score, qt, type = 1)
    as.numeric(cutoff)
}

conformalScore <- function(Y, Yhat, type, side){
    if (type == "CQR" && side == "two"){
        score <- pmax(Yhat[, 1] - Y, Y - Yhat[, ncol(Yhat)])
    } else if (type == "mean" && side == "two"){
        score <- abs(Yhat - Y)
    } else if (side == "above"){
        score <- Y - Yhat
    } else if (side == "below"){
        score <- Yhat - Y
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

