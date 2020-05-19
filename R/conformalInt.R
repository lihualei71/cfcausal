#' @export
conformalInt <- function(X, Y,
                         type = c("CQR", "mean"),
                         lofun = NULL,
                         loquantile = 0.5,
                         loparams = list(),
                         upfun = NULL,
                         upquantile = 0.5,
                         upparams = list(),
                         wtfun = NULL,
                         useCV = FALSE,
                         trainprop = 0.75, trainid = NULL,
                         nfolds = 10, idlist = NULL){
    if (!useCV){
        conformalIntSplit(X, Y,
                          type,
                          lofun, loquantile, loparams,
                          upfun, upquantile, upparams,
                          wtfun,
                          trainprop, trainid)
    } else {
        conformalIntCV(X, Y,
                       type,
                       lofun, loquantile, loparams,
                       upfun, upquantile, upparams,
                       wtfun,
                       nfolds, idlist)
    }
}
