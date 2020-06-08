## Truncate a sequence from both sides
censoring <- function(x, high = 20, low = 0.05){
    pmin(pmax(x, low), high)
}

## Guess the type of a variable
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

## Frequency of each level in a factor
freq <- function(x){
    if (!is.factor(x)){
        x <- as.factor(x)
    }
    as.numeric(table(x)) / length(x)
}

## Divide each row of A by a vector b
row_quo <- function(A, b){
    t(t(A) / b)
}

## Generate a list of indices for cross-validation
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

## Convert a valid outfun string to the function
str_outfun <- function(method){
    if (method == "RF"){
        if (!requireNamespace("randomForest")){
            stop("randomForest package needs to be installed")
        }
        return(RF)
    } else if (method == "quantRF"){
        if (!requireNamespace("grf")){
            stop("grf package needs to be installed")
        }
        return(quantRF)
    } else if (method == "Boosting"){
        if (!requireNamespace("gbm")){
            stop("gbm package needs to be installed")
        }
        return(Boosting)
    } else if (method == "quantBoosting"){
        if (!requireNamespace("gbm")){
            stop("gbm package needs to be installed")
        }
        return(quantBoosting)
    } else if (method == "BART"){
        if (!requireNamespace("bartMachine")){
            stop("bartMachine package needs to be installed")
        }
        return(BART)
    } else if (method == "quantBART"){
        if (!requireNamespace("bartMachine")){
            stop("bartMachine package needs to be installed")
        }
        return(quantBART)
    } else {
        stop(paste0(method, " is not supported. Please input a valid string or a function that meets the minimal requirements described in the man page"))
    }
}

## Convert a valid psfun string to the function
str_psfun <- function(method){
    if (method == "RF"){
        if (!requireNamespace("randomForest")){
            stop("randomForest package needs to be installed")
        }
        return(RF)
    } else if (method == "Boosting"){
        if (!requireNamespace("gbm")){
            stop("gbm package needs to be installed")
        }
        return(Boosting)
    } else {
        stop(paste0(method, " is not supported. Please input a valid string or a function that meets the minimal requirements described in the man page"))
    }
}

## Check if the required inputs are included in outfun
check_outfun <- function(fun, type){
    args <- methods::formalArgs(fun)
    if (type == "CQR"){
        res <- all(c("Y", "X", "Xtest", "quantiles") %in% args)
        if (!res){
            stop("outfun should include 'Y', 'X', 'Xtest' and 'quantiles' as inputs when type = 'CQR'")
        }
    } else {
        res <- all(c("Y", "X", "Xtest") %in% args)
        if (!res){
            stop("outfun should include 'Y', 'X' and 'Xtest' as inputs when type = 'mean'")
        }
    }
}

## Check if the required inputs are included in wtfun
check_wtfun <- function(fun){
    args <- formalArgs(fun)
    if (!("X" %in% args)){
        stop("wtfun should include 'X' as inputs")
    }
}

## Check if the required inputs are included in psfun
check_psfun <- function(fun){
    args <- formalArgs(fun)
    if (!all(c("Y", "X", "Xtest") %in% args)){
        stop("psfun should include 'Y', 'X' and 'Xtest' as inputs")
    }
}
