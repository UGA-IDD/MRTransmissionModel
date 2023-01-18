#' Function to get number of individuals per age group
#'
#' @param trans transition object
#' @param state vector of state of the population at one time point
#' @param tval xxx
#'
#' @return vector size n.age.classes with the number of individuals in each age
#' @export
#'

GetNumber.per.AgeGroup <- function(trans, state, tval=1) {
  epi.class <- trans@epi.class
  pop.per.agegroup <- rep(0, length(trans@age.class))
  for (k in 1:trans@n.epi.class) {
    pop.per.agegroup <- pop.per.agegroup+state[epi.class==k]
  }
  return(pop.per.agegroup)
}
