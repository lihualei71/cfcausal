quantRF <- function(Y, X, Xtest, quantiles, ...){
    fit <- grf::quantile_forest(X, Y, quantiles, ...)
    res <- predict(fit, Xtest, quantiles = quantiles)
    if (length(quantiles) == 1){
        res <- as.numeric(res)
    } else {
        res <- as.matrix(res)
    }
    return(res)
}

RF <- function(Y, X, Xtest, ...){
    fit <- randomForest::randomForest(x = X, y = Y, ...)
    dist <- guessClass(Y)
    if (dist == "gaussian"){
        res <- predict(fit, newdata = Xtest)
        res <- as.numeric(res)
    } else if (dist == "bernoulli"){
        res <- predict(fit, newdata = Xtest, type = "prob")
        res <- as.numeric(res)
    } else if (dist == "multinomial"){
        res <- predict(fit, newdata = Xtest, type = "prob")
        res <- as.matrix(res)
    }
    return(res)
}

quantBoosting <- function(Y, X, Xtest, quantile, n.trees = 100, ...){
    if (class(X) != "data.frame"){
        X <- as.data.frame(X)
        Xtest <- as.data.frame(Xtest)
        names(Xtest) <- names(X)
    }
    data <- data.frame(Y = Y, X)
    fit <- gbm::gbm(Y ~ ., distribution = list(name = "quantile", alpha = quantile), data = data, n.trees = n.trees, ...)
    res <- predict(fit, Xtest, type = "response", n.trees = n.trees)
    return(res)
}

Boosting <- function(Y, X, Xtest, n.trees = 100, ...){
    if (class(X) != "data.frame"){
        X <- as.data.frame(X)
        Xtest <- as.data.frame(Xtest)
        names(Xtest) <- names(X)
    }
    data <- data.frame(Y = Y, X)
    distribution <- guessClass(Y)
    fit <- gbm::gbm(Y ~ ., distribution = distribution, data = data, n.trees = n.trees, ...)
    res <- predict(fit, Xtest, type = "response", n.trees = n.trees)
    if (distribution == "multinomial"){
        res <- matrix(res, nrow = length(Y))
    }
    return(res)
}
