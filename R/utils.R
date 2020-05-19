censoring <- function(x, high = 20, low = 0.05){
    pmin(pmax(x, low), high)
}

guessClass <- function(x){
    if (length(unique(x)) == 2 || is.logical(x) || (is.factor(x) && nlevels(x) == 2)){
        dist <- "bernoulli"
    } else if (is.factor(x) && nlevels(x) > 2){
        dist <- "multinomial"
    } else if (is.numeric(x)){
        dist <- "gaussian"
    }
    return(dist)
}

freq <- function(x){
    if (!is.factor(x)){
        x <- as.factor(x)
    }
    as.numeric(table(x)) / length(x)
}

row_quo <- function(A, b){
    t(t(A) / b)
}

gen_cv_ids <- function(n, nfolds, offset = 0){
    ids <- sample(n, n)
    quo <- floor(n / nfolds)
    if (quo == 0){
        idlist <- lapply(1:nfolds, function(i){
            if (i <= n){
                i
            } else {
                numeric(0)
            }
        })
    } else {
        resid <- n - quo * nfolds
        idlist <- lapply(1:nfolds, function(i){
            tmp <- (i - 1) * quo + 1:quo
            if (i <= resid){
                tmp <- c(tmp, quo * nfolds + i)
            }
            return(ids[tmp] + offset)
        })
    }
    return(idlist)
}
