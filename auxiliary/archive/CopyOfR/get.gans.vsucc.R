#' Get the gans vauccine success object
#'
#' @return xxx
#' @export
#'

get.gans.vsucc <- function() {

  #Data sets on vaccine efficacy
  vdata.gans1998 <- data.frame(age=c(6,9,12), successes=c(10,17,21),
                               N=c(23,20,22))

  return(fit.vsucc.Logistic(vdata.gans1998,.97))
}
