  #' Function to get age-specific rubella seroprevalence per year
#'
#' @param res result data matrix
#' @param trans transition object
#' @param epi.state epidemiological state indexes
#' @param no.gens.in.year xxx
#'
#' @return matrix of seroprevalence by age (in years) per year
#' @export
#'

GetRubellaSeroprevalence.per.Year <- function(res, trans, epi.state, no.gens.in.year){

  #get number of cases per age in year
  sus <- (res[epi.state==2,])  #susceptibles only
  sus.age <- matrix(NA, floor(max(trans@age.class)/12), ncol(sus))
  for (t in 1:ncol(sus)){
    sus.age[,t] <- GetNumber.per.AgeYear(vec=sus[,t], age.classes=trans@age.class)
    sus.age[1,t] <- sum(sus[(9:12),t])
  }

  #sum over each year
  sus.pop <- sus.age[,seq(floor(no.gens.in.year/2),ncol(res),no.gens.in.year)] #number of susceptibles at mid-year
  pop <- GetPop.per.Year(res, trans, no.gens.in.year) #population at mid-year

  seroprev.age <- 1-(sus.pop/pop)
  seroprev.age[which(pop==0)] <- 1 #replace seroprev with 1 if no population in that age group

  return(seroprev.age)
}
