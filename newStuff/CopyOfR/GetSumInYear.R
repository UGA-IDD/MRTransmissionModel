#' Sums over no.gens.in.year to get number per year
#'
#' @param mat matrix with ncol equal to no.gens.year*year
#' @param no.gens.in.year xxx
#'
#' @return xxx
#' @export
#'

GetSumInYear <- function(mat, no.gens.in.year){
  max.year <- ncol(mat)/no.gens.in.year
  out <- matrix(0, nrow(mat), max.year)
  for (y in 1:max.year) {
    for (u in 0:(no.gens.in.year-1)){
      #out[,y] <- out[,y]+mat[,(no.gens.in.year*y+1)-u] #leaves out the first time point, ok through b/c only starting population
      out[,y] <- out[,y]+mat[,(no.gens.in.year*y)-u] #CORRECTION NEED TO MAKE, OTHERWISE VACCINE TIMING BEING MOVED UP A YEAR
    }
  }
  return(out)
}
