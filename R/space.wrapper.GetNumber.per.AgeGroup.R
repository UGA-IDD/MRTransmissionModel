#' Spatial wrapper Function to get number of individuals per age group and subpopulation
#'
#' @param trans transition object
#' @param state vector of state of the population at one time point
#'
#' @return matrix, size n.subpops (rows) by n.ages.class (columns) with the number of individuals in each age and space
#' @export
#'
space.wrapper.GetNumber.per.AgeGroup <- function(trans, state){

  for (s in 1:trans@n.subpops){
    trans.tmp <- trans
    trans.tmp@epi.class <- trans.tmp@epi.class[trans@subpop.class.label==s]
    if (is.vector(state)) tmp <- GetNumber.per.AgeGroup(trans.tmp, state[trans@subpop.class.label==s])
    if (is.matrix(state)) tmp <- GetNumber.per.AgeGroup(trans.tmp, state[trans@subpop.class.label==s,1])
    if (s==1) out <- tmp
    if (s!=1) out <- rbind(out, tmp)
  }
  return(out)

}
