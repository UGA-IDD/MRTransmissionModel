#' Function to get age-specific population size at mid-year
#'
#' @param res result data matrix
#' @param trans transition object
#' @param no.gens.in.year xxx
#'
#' @return matrix of pop size by age (in years) per year
#' @export
#'

GetPop.per.Year <- function(res, trans, no.gens.in.year){

  mid.year.index <- seq(13,ncol(res)-1, no.gens.in.year) #mid-year every year
  pop.year <- matrix(NA, 101, length(mid.year.index))
  for (i in 1:length(mid.year.index)){
    t <- mid.year.index[i]
    state <- res[,t]
    pop.all <- GetNumber.per.AgeGroup(state=state, trans=trans)
    pop.year[,i] <- GetNumber.per.AgeYear(vec=pop.all, age.classes=trans@age.class)
  }

  return(pop.year)
}
