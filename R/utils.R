censoring <- function(x, high = 20, low = 0.05){
    pmin(pmax(x, high), low)
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
