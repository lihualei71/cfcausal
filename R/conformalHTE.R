conformalHTE_CI <- function(X, Y, T,
                            Xtau, Ytau, Ttau,
                            type, side,
                            alpha,
                            outfun, outparams,
                            wtfun1, wtfun0,
                            trainprop, trainid,
                            wthigh, wtlow){
    X1 <- X[T == 1, ]
    Y1 <- Y[T == 1]
    X0 <- X[T == 0, ]
    Y0 <- Y[T == 0]

    n <- nrow(Xtau)
    if (side == "two"){
        side1 <- side0 <- "two"
    } else if (side == "above"){
        side1 <- "above"
        side0 <- "below"
    } else if (side == "below"){
        side1 <- "below"
        side0 <- "above"
    }
    obj1 <- conformal(X1, Y1, Xtau, type, side1, alpha,
                      outfun, outparams, wtfun1,
                      trainprop, trainid)
    obj0 <- conformal(X0, Y0, Xtau, type, side0, alpha,
                      outfun, outparams, wtfun0,
                      trainprop, trainid)
    CI1 <- predict(obj1, NULL, alpha, wthigh = wthigh, wtlow = wtlow)
    CI0 <- predict(obj0, NULL, alpha, wthigh = wthigh, wtlow = wtlow)
    if (!is.null(Ytau) && !is.null(Ttau)){
        tau_high <- tau_low <- rep(0, n)
        tau_high[Ttau == 1] <- Ytau[Ttau == 1] - CI0$Ylow[Ttau == 1]
        tau_high[Ttau == 0] <- CI1$Yhigh[Ttau == 0] - Ytau[Ttau == 0]
        tau_low[Ttau == 1] <- Ytau[Ttau == 1] - CI0$Yhigh[Ttau == 1]
        tau_low[Ttau == 0] <- CI1$Ylow[Ttau == 0] - Ytau[Ttau == 0]
        CItau <- data.frame(tau_low = tau_low,
                            tau_high = tau_high)
    } else {
        CItau <- NULL
    }
    return(list(CItau = CItau, CIY1 = CI1, CIY0 = CI0,
                Y1model = obj1, Y0model = obj0))
}

#' @export
conformalHTE <- function(X, Y, T,
                         Xtau = NULL,
                         Ytau = NULL,
                         Ttau = NULL,
                         Xtest = NULL,
                         estimand = c("ATT", "ATE"),
                         type = c("CQR", "mean"),
                         side = c("two", "above", "below"),
                         alpha = 0.05, alpha2 = 0.5,
                         outfun = NULL,
                         outparams = list(),
                         psfun = NULL,
                         psparams = list(),
                         wthigh = 20, wtlow = 0.05,
                         HTEtype = c("None", "double", "single"),
                         useCV = FALSE,
                         nfolds = 5, foldid = NULL,
                         trainprop = 0.75, trainid = NULL){
    if (useCV && (!is.null(Xtau) || !is.null(Ytau) || !is.null(Ttau))){
        warning("Xtau, Ytau, and Ttau are dropped when CV is used")
    }
    if (!useCV && (is.null(Xtau) || is.null(Ytau) || is.null(Ttau))){
        stop("Xtau, Ytau, and Ttau must be given when CV is not used")
    }
    if (!is.null(Xtest) && (HTEtype == "None")){
        warning("Xtest is dropped when HTEtype = None")
    }
    if (HTEtype == "single" && is.null(Xtest)){
        stop("Xtest must be given when HTEtype = single")
    }
    if (is.null(outfun)){
        outfun <- switch(type,
                         CQR = quantRF,
                         mean = RF)
    }
    if (is.null(psfun)){
        psfun <- Boosting
    }
    estimand <- estimand[1]
    stopifnot(estimand %in% c("ATT", "ATE"))
    type <- type[1]
    stopifnot(type %in% c("CQR", "mean"))
    side <- side[1]
    stopifnot(side %in% c("two", "above", "below"))

    psparams <- c(list(Y = as.numeric(T), X = X), psparams)
    psmodel <- function(Xtest){
        do.call(psfun, c(psparams, list(Xtest = Xtest)))
    }
    if (estimand == "ATE"){
        wtfun1 <- function(Xtest){
            1 / psmodel(Xtest)
        }
        wtfun0 <- function(Xtest){
            1 / (1 - psmodel(Xtest))
        }
    } else if (estimand == "ATT"){
        wtfun1 <- function(Xtest){
            rep(1, nrow(Xtest))
        }
        wtfun0 <- function(Xtest){
            tmp <- psmodel(Xtest)
            tmp / (1 - tmp)
        }
    }

    if (!useCV){
        CItau <- conformalHTE_CI(
            X, Y, T,
            Xtau, Ytau, Ttau,
            type, side,
            alpha / 2,
            outfun, outparams,
            wtfun1, wtfun0,
            trainprop, trainid,
            wthigh, wtlow)$CItau
        CItau <- cbind(CItau, Xtau)
    } else {
        id1 <- sample(which(T == 1), sum(T == 1))
        id0 <- sample(which(T == 0), sum(T == 0))
        n1 <- length(id1)
        n0 <- length(id0)
        n <- n1 + n0
        bsize1 <- floor(n1 / nfolds)
        bsize0 <- floor(n0 / nfolds)
        ids <- lapply(1:nfolds, function(i){
            start1 <- (i - 1) * bsize1 + 1
            start0 <- (i - 1) * bsize0 + 1
            end1 <- ifelse(i == nfolds, n1, i * bsize1)
            end0 <- ifelse(i == nfolds, n0, i * bsize0)
            c(id1[start1:end1], id0[start0:end0])
        })
        CItau <- data.frame(tau_low = rep(NA, n),
                            tau_high = rep(NA, n))
        for (k in 1:nfolds){
            id <- ids[[k]]
            CItau[id, ] <- conformalHTE_CI(
                X[-id, ], Y[-id], T[-id],
                X[id, ], Y[id], T[id],
                type, side,
                alpha / 2,
                outfun, outparams,
                wtfun1, wtfun0,
                trainprop, trainid,
                wthigh, wtlow)$CItau
        }
        CItau <- cbind(CItau, X)
    }

    if (HTEtype == "None"){
        HTE_model <- NULL
        HTE_test <- NULL
    } else if (HTEtype == "single"){
        if (!useCV){
            X <- cbind(X, Xtau)
            Y <- c(Y, Ytau)
            T <- c(T, Ttau)
        }
        res <- conformalHTE_CI(
            X, Y, T,
            Xtest, NULL, NULL,
            type, side,
            alpha / 2,
            outfun, outparams,
            wtfun1, wtfun0,
            trainprop, trainid,
            wthigh, wtlow)
        tau_high <- res$CIY1$Yhigh - res$CIY0$Ylow
        tau_low <- res$CIY1$Ylow - res$CIY0$Yhigh
        HTE_test <- data.frame(tau_low = tau_low,
                               tau_high = tau_high)
        HTE_model <- function(Xtest){
            CI1 <- predict(res$Y1model, Xtest, alpha / 2, wthigh, wtlow)
            CI0 <- predict(res$Y0model, Xtest, alpha / 2, wthigh, wtlow)
            tau_high <- CI1$Yhigh - CI0$Ylow
            tau_low <- CI1$Ylow - CI0$Yhigh
            data.frame(tau_low = tau_low,
                       tau_high = tau_high)
        }
    } else if (HTEtype == "double"){
        X <- as.matrix(CItau[, -c(1, 2)])
        if (side != "below"){
            tau_high <- CItau$tau_high
            alpha_high <- ifelse(side == "two", alpha2 / 2, alpha2)
            obj_high <- conformal(
                X, tau_high, NULL,
                type, "above",
                alpha_high,
                outfun, outparams,
                NULL,
                trainprop, NULL)
        }
        if (side != "above"){
            tau_low <- CItau$tau_low
            alpha_low <- ifelse(side == "two", alpha2 / 2, alpha2)
            obj_low <- conformal(
                X, tau_low, NULL,
                type, "below",
                alpha_low,
                outfun, outparams,
                NULL,
                trainprop, NULL)
        }

        if (side == "two"){
            HTE_model <- function(Xtest){
                CI_high <- predict(obj_high, Xtest, alpha2 / 2, wthigh, wtlow)$Yhigh
                CI_low <- predict(obj_low, Xtest, alpha2 / 2, wthigh, wtlow)$Ylow
                data.frame(low = CI_low, high = CI_high)
            }
        } else if (side == "above"){
            HTE_model <- function(Xtest){
                CI_high <- predict(obj_high, Xtest, alpha2, wthigh, wtlow)$Yhigh
                data.frame(low = -Inf, high = CI_high)
            }
        } else if (side == "below"){
            HTE_model <- function(Xtest){
                CI_low <- predict(obj_low, Xtest, alpha2, wthigh, wtlow)$Ylow
                data.frame(low = CI_low, high = Inf)
            }
        }

        if (!is.null(Xtest)){
            HTE_test <- HTE_model(Xtest)
        } else {
            HTE_test <- NULL
        }
    }

    return(list(CItau = CItau, HTE_model = HTE_model, HTE_test = HTE_test))
}
