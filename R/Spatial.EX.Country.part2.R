#' Function to Run the SPATIAL Simulations/Experiments
#'
#' @param uncode UN country code
#' @param generation.time generation time in months
#' @param pop.rescale numeric - population by which you want to rescale at pop.time
#' @param pop.time numeric - time in YEARS that you want to rescale the population by
#' @param is.stochastic logical
#' @param get.births vector of births with length of 5*n.age.classes
#' @param t.max numeric - time in year that plan to run the experiment
#' @param rescale.WAIFW logical
#' @param yr.births.per.1000.acrossyears.bysupop vector of length(t.max) - crude birth rate per 1000 per year
#' @param space.asdr.object space nMx object - with rate and mid-age (age specific death rates)
#' @param year numeric - start year
#' @param EXt0 object from part1
#' @param time.specific.MR1cov - numeric matrix
#' @param time.specific.MR2cov - numeric matrix
#' @param time.specific.SIAcov - numeric matrix
#' @param time.specific.schoolvacc.cov - numeric matrix
#' @param time.specific.min.age.MR1 - numeric vector
#' @param time.specific.max.age.MR1 - numeric vector
#' @param time.specific.min.age.MR2 - numeric vector
#' @param time.specific.max.age.MR2 - numeric vector
#' @param time.specific.min.age.SIA - numeric vector
#' @param time.specific.max.age.SIA - numeric vector
#' @param time.specific.min.age.schoolvacc - numeric vector
#' @param time.specific.max.age.schoolvacc - numeric vector
#' @param obj.vcdf.MR1 - object
#' @param obj.vcdf.MR2 - object
#' @param list.obj.vcdf.schoolenroll - list of objects
#' @param obj.prob.vsucc - object
#' @param sia.timing.in.year - numeric scalar
#' @param schoolvacc.timing.in.year - numeric scalar
#' @param MR1MR2correlation - NULL or not
#' @param MR1SIAcorrelation - NULL or not
#' @param MR2SIAcorrelation - NULL or not
#' @param SIAinacc - xxx
#' @param prop.inacc - xxx
#' @param SIAinefficient - xxx
#'
#' @return experiment result
#' @export
#'
Spatial.EX.Country.part2 <- function(uncode,
                                     generation.time = 0.5, #generation time in months
                                     pop.rescale=NULL,
                                     pop.time=4,
                                     is.stochastic=FALSE,
                                     get.births=NULL,
                                     t.max = 40,
                                     rescale.WAIFW=F,
                                     yr.births.per.1000.acrossyears.bysupop,
                                     space.asdr.object,
                                     year=1990,
                                     EXt0=EXt0,
                                     time.specific.MR1cov,
                                     time.specific.MR2cov,
                                     time.specific.SIAcov,
                                     time.specific.schoolvacc.cov=NULL,
                                     time.specific.min.age.MR1,
                                     time.specific.max.age.MR1,
                                     time.specific.min.age.MR2,
                                     time.specific.max.age.MR2,
                                     time.specific.min.age.SIA,
                                     time.specific.max.age.SIA,
                                     time.specific.min.age.schoolvacc,
                                     time.specific.max.age.schoolvacc,
                                     obj.vcdf.MR1 = get.vcdf.normal(6, 12),
                                     obj.vcdf.MR2 = get.vcdf.normal(15, 21),
                                     list.obj.vcdf.schoolenroll,
                                     obj.prob.vsucc = pvacsuccess(1:(20*12), get.boulianne.vsucc()),
                                     sia.timing.in.year = (3/12),
                                     schoolvacc.timing.in.year = (9/12),
                                     MR1MR2correlation = NULL,
                                     MR1SIAcorrelation = NULL,
                                     MR2SIAcorrelation = NULL,
                                     SIAinacc = NULL,
                                     prop.inacc = NULL,
                                     SIAinefficient = NULL){


  ## Changing experiment type
  if (is.null(time.specific.schoolvacc.cov)) EX <- new("experiment.updatedemog.vaccinationchange.spatial") #hard coded MR1 and MR2 dependent
  if (!is.null(time.specific.schoolvacc.cov)) EX <- new("experiment.updatedemog.vaccinationchange.school.spatial") #hard coded MR1 and MR2 dependent
  #if (!is.null(MR1MR2correlation)) EX <- new("experiment.updatedemog.vaccinationchange.vaccinationlimitations.spatial") #this experiment is not set up yet

  ## Country name
  name <- countrycode::countrycode(uncode, origin="un", destination="country.name")

  ## Update new experiment
  EX@trans <- EXt0@trans
  EX@state.t0 <- EXt0@state.t0
  EX@maternal.obj <- EXt0@maternal.obj
  EX@name <- paste("Vaccination Experiment:", name, year, "to", (year+t.max), sep=" ")
  EX@description <- "Frequency Dependent Stochastic"
  EX@t.min <- 0
  EX@t0.doy = 0
  EX@step.size <- EXt0@step.size #1/no.gens.per.year
  EX@season.obj <- EXt0@season.obj

  ## Add in the new features
  # Time horizon
  EX@t.max <- t.max

  ## Get time steps information
  no.time.steps.in.experiment <- round((EX@t.max-EX@t.min)/EX@step.size)+1
  no.gens.in.year <- (no.time.steps.in.experiment-1)/t.max

  # Make stochastic if so desired
  EX@trans@is.stochastic <- is.stochastic
  if (is.stochastic)  EX@state.t0[,1] <- round(EX@state.t0[,1])

  # Population by which to rescale and when, if it is called for
  EX@pop.rescale.each.timestep <- matrix(NaN, EX@trans@n.subpops, no.time.steps.in.experiment)
  if (!is.null(pop.rescale)) { #if pop.rescale is not NULL fill in the pops to rescale at mid-year
    for (p in 1:ncol(pop.rescale)){
      EX@pop.rescale.each.timestep[,(pop.time[p]*no.gens.in.year-1)] <- pop.rescale[,p] #rescalling at end of previous year
    }
  }

  # Setting up changing birth rates
  EX@births.per.1000.each.timestep <- cbind(yr.births.per.1000.acrossyears.bysupop[,1],
                                            matrix(rep(t(yr.births.per.1000.acrossyears.bysupop), each = no.gens.in.year), byrow=T, nrow=EX@trans@n.subpops))*generation.time/12

  # Setting up survival
  EX@surv.each.timestep <- space.wrapper.surv.prob.over.age.time(EX@trans@age.class, generation.time, space.nMx=space.asdr.object, years.interpolate =seq(year,(year+t.max-1),1)) #xxamy minus 1 to get correct number of time steps

  # Setting up vaccination pieces
  EX@time.specific.MR1cov =  time.specific.MR1cov #matrix annual coverage values by space (rows) and years 1980-2100 (cols)
  EX@time.specific.MR2cov = time.specific.MR2cov #matrix annual coverage values by space (rows) and years 1980-2100 (cols)
  EX@time.specific.SIAcov =   time.specific.SIAcov #matrix annual coverage values by space (rows) and years 1980-2100 (cols)
  EX@time.specific.min.age.MR1 = time.specific.min.age.MR1
  EX@time.specific.max.age.MR1 = time.specific.max.age.MR1
  EX@time.specific.min.age.MR2 = time.specific.min.age.MR2
  EX@time.specific.max.age.MR2 = time.specific.max.age.MR2
  EX@time.specific.min.age.SIA = time.specific.min.age.SIA
  EX@time.specific.max.age.SIA = time.specific.max.age.SIA
  EX@obj.vcdf.MR1 = obj.vcdf.MR1
  EX@obj.vcdf.MR2 = obj.vcdf.MR2
  EX@obj.prob.vsucc = obj.prob.vsucc
  EX@sia.timing.in.year = sia.timing.in.year
  if (!is.null(MR1MR2correlation)) {
    EX@MR1MR2correlation = MR1MR2correlation
    EX@MR1SIAcorrelation = MR1SIAcorrelation
    EX@MR2SIAcorrelation = MR2SIAcorrelation
    EX@SIAinacc = SIAinacc
    #Setting up population always inaccessible to vaccine
    EX@prop.inacc <- c(rep(prop.inacc, each=no.gens.in.year),  prop.inacc[length(prop.inacc)]) #by time step
    #EX@prop.inacc = prop.inacc #by year
    EX@SIAinefficient = SIAinefficient
  }
  if (!is.null(time.specific.schoolvacc.cov)) {
    EX@time.specific.schoolvacc.cov =  time.specific.schoolvacc.cov #matrix annual coverage values by space (rows) and years 1980-2100 (cols)
    EX@time.specific.min.age.schoolvacc = time.specific.min.age.schoolvacc
    EX@time.specific.max.age.schoolvacc = time.specific.max.age.schoolvacc
    EX@list.obj.vcdf.schoolenroll <- list.obj.vcdf.schoolenroll
    EX@schoolvacc.timing.in.year <- schoolvacc.timing.in.year
  }

  return(run(EX, rescale.WAIFW=rescale.WAIFW))

}
