#' @export
conformal <- function(X, Y,
                      type = c("CQR", "mean"),
                      side = c("two", "above", "below"),
                      quantiles = NULL,
                      outfun = NULL,
                      outparams = list(),
                      wtfun = NULL,
                      useCV = FALSE,
                      trainprop = 0.75, trainid = NULL,
                      nfolds = 10, idlist = NULL){
    if (!useCV){
        conformalSplit(X, Y,
                       type, side,
                       quantiles,
                       outfun, outparams,
                       wtfun,
                       trainprop, trainid)
    } else {
        conformalCV(X, Y,
                    type, side,
                    quantiles,
                    outfun, outparams,
                    wtfun,
                    nfolds, idlist)
    }
}
