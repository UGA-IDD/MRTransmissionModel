#' Function to get age-specific rubella cases size per year
#'
#' @param res result data matrix
#' @param trans transition object
#' @param epi.state epidemiological state indexes
#' @param no.gens.in.year xxx
#'
#' @return matrix of pop size by age (in years) per year
#' @export
#'

GetRubellaCases.per.Year <- function(res, trans, epi.state, no.gens.in.year){

  #get number of cases per age in year
  inf <- (res[epi.state==3,])  #cases only
  inf.age <- matrix(NA, 101, ncol(inf))
  for (t in 1:ncol(inf)){
    inf.age[,t] <- GetNumber.per.AgeYear(vec=inf[,t], age.classes=trans@age.class)
  }

  #sum over each year
  inf.pop <- GetSumInYear(inf.age, no.gens.in.year=no.gens.in.year)

  return(inf.pop)
}
