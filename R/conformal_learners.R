## quantile random forest. grf package needed
quantRF <- function(Y, X, Xtest, quantiles, ...){
    fit <- grf::quantile_forest(X, Y, quantiles = quantiles, ...)
    res <- predict(fit, Xtest, quantiles = quantiles)
    if (is.list(res) && !is.data.frame(res)){
        res <- res$predictions
    }
    if (length(quantiles) == 1){
        res <- as.numeric(res)
    } else {
        res <- as.matrix(res)
    }
    return(res)
}

## random forest. randomForest package needed
RF <- function(Y, X, Xtest, ...){
    dist <- guessClass(Y)
    if (dist == "gaussian"){
        fit <- randomForest::randomForest(x = X, y = Y, ...)
        res <- predict(fit, newdata = Xtest)
        res <- as.numeric(res)
    } else if (dist == "bernoulli"){
        if (!is.factor(Y)){
            Y <- as.factor(Y)
        }
        fit <- randomForest::randomForest(x = X, y = Y, ...)
        res <- predict(fit, newdata = Xtest, type = "prob")
        res <- as.numeric(res[, 2])
    } else if (dist == "multinomial"){
        if (!is.factor(Y)){
            Y <- as.factor(Y)
        }
        fit <- randomForest::randomForest(x = X, y = Y, ...)
        res <- predict(fit, newdata = Xtest, type = "prob")
        res <- as.matrix(res)
    }
    return(res)
}

## quantile gradient boosting. gbm package needed
quantBoosting <- function(Y, X, Xtest, quantiles, n.trees = 100, ...){
    if (class(X)[1] != "data.frame"){
        X <- as.data.frame(X)
        Xtest <- as.data.frame(Xtest)
        names(Xtest) <- names(X)
    }
    data <- data.frame(Y = Y, X)
    fit <- gbm::gbm(Y ~ ., distribution = list(name = "quantile", alpha = quantiles[1]), data = data, n.trees = n.trees, ...)
    res <- predict(fit, Xtest, type = "response", n.trees = n.trees)
    if (length(quantiles) == 2){
        fit2 <- gbm::gbm(Y ~ ., distribution = list(name = "quantile", alpha = quantiles[2]), data = data, n.trees = n.trees, ...)
        res2 <- predict(fit2, Xtest, type = "response", n.trees = n.trees)
        res <- cbind(res, res2)
    }
    return(res)
}

## gradient boosting. gbm package needed
Boosting <- function(Y, X, Xtest, n.trees = 100, ...){
    if (class(X)[1] != "data.frame"){
        X <- as.data.frame(X)
        Xtest <- as.data.frame(Xtest)
        names(Xtest) <- names(X)
    }
    distribution <- guessClass(Y)
    if (distribution == "bernoulli" && is.factor(Y)){
        Y <- as.numeric(Y) - 1
    }
    data <- data.frame(Y = Y, X)
    fit <- gbm::gbm(Y ~ ., distribution = distribution, data = data, n.trees = n.trees, ...)
    res <- predict(fit, Xtest, type = "response", n.trees = n.trees)
    if (distribution == "multinomial"){
        res <- matrix(res, nrow = nrow(Xtest))
    }
    return(res)
}

## posterior quantiles of BART. bartMachine package needed
quantBART <- function(Y, X, Xtest, quantiles,
                      ndpost = 100, ...){
    if (class(X)[1] != "data.frame"){
        X <- as.data.frame(X)
        Xtest <- as.data.frame(Xtest)
        names(Xtest) <- names(X)
    }
    y <- Y
    fit <- bartMachine::bartMachine(X, y, verbose = FALSE)
    if (length(quantiles) == 2){
        if (sum(quantiles) != 1){
            warning("Two quantiles should sum up to 1.")
        }
        ci_conf <- quantiles[2] - quantiles[1]
        res <- bartMachine::calc_prediction_intervals(
                                fit, new_data = Xtest,
                                pi_conf = 0.95)$interval
        res <- as.matrix(res)
    } else if (length(quantiles) == 1){
        if (quantiles[1] > 0.5){
            ci_conf <- 2 * quantiles[1]
            res <- bartMachine::calc_prediction_intervals(
                                    fit, new_data = Xtest,
                                    pi_conf = 0.95)$interval[, 2]
            res <- as.numeric(res)
        }  else{
            ci_conf <- 2 * (1 - quantiles[1])
            res <- bartMachine::calc_prediction_intervals(
                                    fit, new_data = Xtest,
                                    pi_conf = 0.95)$interval[, 1]
            res <- as.numeric(res)
        }
    }
    return(res)
}

## BART. bartMachine package needed
BART <- function(Y, X, Xtest,
                 ndpost = 100, ...){
    if (class(X)[1] != "data.frame"){
        X <- as.data.frame(X)
        Xtest <- as.data.frame(Xtest)
        names(Xtest) <- names(X)
    }
    y <- Y
    fit <- bartMachine::bartMachine(X, y, verbose = FALSE)
    res <- predict(fit, Xtest)
    return(res)
}
