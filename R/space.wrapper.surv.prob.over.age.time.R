#' Wrapper function to add space to create.surv.prob.over.age.time()
#'
#' @param age.classes vector; a set of age classes for which tran is being built
#' @param generation.time numeric; generation time
#' @param space.nMx space.nMx object; age and space specific death rates over time
#' @param years.interpolate vector; years to interpolate age specific death rates
#' @param check boolean; whether or not to plot the death rates over age
#'
#' @return returns the age and spatial profile of survivorship (in units of the generation time) (rows) by time (cols)
#' @export
#'
space.wrapper.surv.prob.over.age.time <- function(age.classes, generation.time,
                                                  space.nMx=NULL, years.interpolate = seq(1980,(1980+121),1),
                                                  check=FALSE){

  n.subpops <- space.nMx@n.subpops

  for (s in 1:n.subpops){
    n.age.classes <- length(space.nMx@mid.age)
    nMx.tmp <- new("nMx",
                   rate.years = space.nMx@rate.years,
                   rates = space.nMx@rates[((s*n.age.classes)-n.age.classes+1):(s*n.age.classes),],
                   mid.age = space.nMx@mid.age)
    surv.tmp <- create.surv.prob.over.age.time(age.classes, generation.time, nMx = nMx.tmp, nMx.years = years.interpolate, check=check)
    if (s==1) surv <- surv.tmp
    if (s!=1) surv <- rbind(surv, surv.tmp)
  }

  return(surv)

}
