#' @export
conformalIte <- function(X, Y, T,
                         alpha = 0.1,
                         algo = c("nest",
                                  "naive",
                                  "counterfactual"),
                         exact = FALSE,
                         type = c("CQR", "mean"),
                         side = c("two", "above", "below"),
                         quantiles = NULL,
                         outfun = NULL,
                         outparams = list(),
                         psfun = NULL,
                         psparams = list(),
                         cfprop = 0.5,
                         citype = c("CQR", "mean"),
                         lofun = NULL,
                         loquantile = 0.5,
                         loparams = list(),
                         upfun = NULL,
                         upquantile = 0.5,
                         upparams = list(),
                         useCV = FALSE,
                         trainprop = 0.75,
                         nfolds = 10,
                         wthigh = 20, wtlow = 0.05){
    algo <- algo[1]
    stopifnot(algo %in% c("nest", "naive", "counterfactual"))
    if (algo == "nest"){
        obj <- conformalIteNest(X, Y, T,
                                alpha,
                                type, side,
                                quantiles,
                                outfun, outparams,
                                psfun, psparams,
                                exact,
                                cfprop, citype,
                                lofun, loquantile, loparams,
                                upfun, upquantile, upparams,
                                useCV,
                                trainprop, nfolds,
                                wthigh, wtlow)
        function(X){
            predict(obj, X)
        }
    } else if (algo == "naive"){
        obj <- conformalIteNaive(X, Y, T,
                                 type, side,
                                 quantiles,
                                 outfun, outparams,
                                 psfun, psparams,
                                 useCV,
                                 trainprop, nfolds)
        function(X){
            predict(obj, X, alpha, wthigh, wtlow)$Ite
        }
    } else if (algo == "counterfactual"){
        obj <- conformalIteCf(X, Y, T,
                              type, side,
                              quantiles,
                              outfun, outparams,
                              psfun, psparams,
                              useCV,
                              trainprop, nfolds)
        function(X, Y, T){
            predict(obj, X, Y, T, alpha, wthigh, wtlow)
        }
    }
}
