#' @export
conformalClass <- function(X, Y,
                           type = c("weighted", "unweighted"),
                           clsp = FALSE,
                           outfun = NULL,
                           outparams = list(),
                           wtfun = NULL,
                           useCV = FALSE,
                           trainprop = 0.75, trainid = NULL,
                           nfolds = 10, idlist = NULL){
    if (!useCV){
        conformalClassSplit(X, Y,
                            type, clsp,
                            outfun, outparams,
                            wtfun,
                            trainprop, trainid)
    } else {
        conformalClassCV(X, Y,
                         type, clsp,
                         outfun, outparams,
                         wtfun,
                         nfolds, idlist)
    }
}
