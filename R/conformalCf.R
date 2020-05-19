#' @export
conformalCf <- function(X, Y, 
                        estimand = c("unconditional",
                                     "nonmissing",
                                     "missing"),
                        type = c("CQR", "mean"),
                        side = c("two", "above", "below"),
                        quantiles = NULL,
                        outfun = NULL,
                        outparams = list(),
                        psfun = NULL,
                        psparams = list(),
                        useCV = FALSE,
                        trainprop = 0.75,
                        nfolds = 10){
    if (!useCV){
        conformalCf_split(X, Y,
                          estimand,
                          type, side,
                          quantiles,
                          outfun, outparams,
                          psfun, psparams,
                          trainprop)
    } else {
        conformalCf_CV(X, Y,
                       estimand,
                       type, side,
                       quantiles,
                       outfun, outparams,
                       psfun, psparams,
                       nfolds)
    }
}
