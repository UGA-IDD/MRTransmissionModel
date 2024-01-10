#' Get the Boulianne vaccine success object (rubella)
#'
#' @return xxx
#' @export
#'

get.boulianne.vsucc <- function() {

  #rubella with MMR
  vdata.MMR.boulianne1995 <- data.frame(age=c(12,13,14,15), successes=c(62,74,50,47),
                                        N=c(64,79,50,48))

  return(fit.vsucc.Logistic(vdata.MMR.boulianne1995,.97))
}
