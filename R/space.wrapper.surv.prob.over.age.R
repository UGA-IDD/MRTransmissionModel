#' Wrapper function to add space to create.surv.prob.over.age
#'
#' @param age.classes vector; a set of age classes for which tran is being built
#' @param generation.time numeric; generation time
#' @param space.nMx space.nMx object; age and space specific death rates over time
#' @param year year to pull age specific death rates that will be used to run out the transients (should coincide with DFE year)
#' @param check logical;
#'
#' @return returns the age profile (cols) of survivorship (in units of the generation time) by spatial unit (rows)
#' @export
#'
space.wrapper.surv.prob.over.age <- function(age.classes, generation.time,
                                             space.nMx=NULL, year = 1980, check=FALSE){

  n.subpops <- space.nMx@n.subpops

  for (s in 1:n.subpops){
    n.age.classes <- length(space.nMx@mid.age)
    nMx.tmp <- new("nMx",
                   rate.years = space.nMx@rate.years,
                   rates = space.nMx@rates[((s*n.age.classes)-n.age.classes+1):(s*n.age.classes),],
                   mid.age = space.nMx@mid.age)
    surv.tmp <- create.surv.prob.over.age(age.classes, generation.time, nMx = nMx.tmp, year, check)
    if (s==1) surv <- surv.tmp
    if (s!=1) surv <- rbind(surv, surv.tmp)
  }

  return(surv)

}
