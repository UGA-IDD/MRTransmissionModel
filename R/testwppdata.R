#' Test data calls to wpp package
#'
#' @return mean of 2025 data
#' @export
#'
testwppdata <- function(){

  data("popproj", package = "wpp2017", envir = environment()) # wpp2017
  return(mean(popproj$`2025`))

}
