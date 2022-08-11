#' Function to set up experiment and run the transients out
#'
#' @param uncode UN country code
#' @param generation.time the desired generation time in months
#' @param age.classes vector; the upper limit of the age classes we are interested in months
#' @param maternal.decay.rt rate of maternal decay of immunity to rubella
#' @param exponent numeric; exponent for the infected
#' @param frequency.dep xxx
#' @param tot.pop numeric; population you want to scale you whole experiment by at the beginning, NULL or numeric
#' @param yr.births.per.1000.acrossyears vector of length t.max; crude birth rate per 1000 per year
#' @param intro.rate numeric; rate at which new infections are introduced
#' @param targeted.intro logical; FALSE means space introduction out over all age classes, TRUE means to concentrate introductions
#' @param R0 numeric; basic reproductive number of rubella assumed, used to rescale the waifw
#' @param t.max numeric; time in year that plan to run out transients
#' @param get.births vector of births with length of 5*n.age.classes
#' @param seasonal.amp numeric; amplitude of the seasonal cosine function
#' @param flat.WAIFW logical; if FALSE then default waifw is polymod GB
#' @param country.specific.WAIFW logical; if true then overrides polymod or flat waifw
#' @param vynnycky.waifw logical; true if want emilia's waifw <13 and 15-50 (beta_young, beta_old, and 0.7*beta_old)
#' @param vynnycky.waifw.betas vector of two betas (beta_young and beta_old)
#' @param asdr.object nMx object with rate and mid-age (age specific death rates)
#' @param year numeric; year to pull DFE demography (population and age structure)
#' @param use_montagu_demog logical; T = use montagu demography or F = use UNPD demography
#' @param routine.vac xxx
#' @param routine.vac.age.index xxx
#' @param demog_data optional demography data for use if demog_data = NULL
#'
#' @include setClasses.R
#' @importFrom methods new
#'
#' @return experiment result
#' @export

EX.Country.part1 <- function(uncode,
                             generation.time = 0.5, #generation time in months
                             age.classes = c(1:240, seq(241,720,12)),
                             maternal.decay.rt=0.95, #based on Metcalf, 2012 paper
                             exponent = 0.97,
                             frequency.dep=TRUE,
                             tot.pop=NULL,
                             yr.births.per.1000.acrossyears,
                             intro.rate=1e-10,
                             targeted.intro=FALSE,
                             R0=5,
                             t.max = 20,
                             get.births=NULL,
                             seasonal.amp=0.35,#metcalf et al 2012 used 0.35
                             flat.WAIFW=FALSE,
                             country.specific.WAIFW=FALSE,
                             vynnycky.waifw=FALSE,
                             vynnycky.waifw.betas = c(2,1),
                             asdr.object,
                             year=1990,
                             use_montagu_demog=FALSE,
                             routine.vac=0,
                             routine.vac.age.index=12,
                             demog_data = NULL) {

  # build.R
  tmp <- Get.CountryX.Starting.Pop.MSIRV(uncode=uncode,
                                         generation.time=generation.time,
                                         age.classes=age.classes,
                                         maternal.decay.rt=maternal.decay.rt,
                                         exponent=exponent,
                                         frequency.dep = frequency.dep,
                                         is.stochastic = FALSE, # Make NOT stochastic for the transient ridding stage
                                         tot.pop=tot.pop,
                                         yr.births.per.1000 = yr.births.per.1000.acrossyears[1],
                                         intro.rate=intro.rate,
                                         flat.WAIFW=flat.WAIFW,
                                         asdr.object=asdr.object,
                                         targeted.intro=targeted.intro,
                                         year=year,
                                         get.births=get.births,
                                         use_montagu_demog=use_montagu_demog,
                                         routine.vac=routine.vac,
                                         routine.vac.age.index=routine.vac.age.index,
                                         demog_data = demog_data)

  EX <- new("experiment.updatedemog")
  EX@name <- paste(uncode, "Frequency Dependent", sep=" ")
  EX@state.t0 <- tmp$state
  EX@trans <- tmp$tran
  EX@t.min <- 0
  EX@t0.doy <- 0
  EX@step.size <- 1/(12*1/generation.time) #1/no.gens.per.year
  EX@season.obj <- new("seasonal.cosine", amplitude = seasonal.amp)
  EX@R0 <- R0
  EX@t.max <-t.max
  EX@maternal.obj <- new("maternal.exp.decay", decay.rt=maternal.decay.rt)

  #change the WAIFW to country-specific if requested (based on prem et al.)

  # function "get.prem.WAIFW()" in build.R
  if (country.specific.WAIFW & use_montagu_demog){
    EX@trans@waifw <- get.prem.WAIFW(age.class.boundries=age.classes/12,
                                     uncode, other.contact.matrix=F,
                                     bandwidth=c(3,3),
                                     adjustment_start_year=T, year=1980)
  }

  if (country.specific.WAIFW & !use_montagu_demog){
    EX@trans@waifw <- get.prem.WAIFW(age.class.boundries=age.classes/12,
                                     uncode, other.contact.matrix=F,
                                     bandwidth=c(3,3),
                                     adjustment_start_year=F, year=1980)
  }

  if (vynnycky.waifw) {
    EX@trans@waifw <- get.vynnycky.WAIFW(age.classes/12,
                                         beta_young=vynnycky.waifw.betas[1],
                                         beta_old=vynnycky.waifw.betas[2])
  }

  # To run out the transients, keep the birth rate constant across time
  EX@births.per.1000.each.timestep <- c(yr.births.per.1000.acrossyears[1],
                                        rep(yr.births.per.1000.acrossyears,
                                            each=12*1/generation.time))*generation.time/12
  # Setting up survival over time as NULL because not changing over time
  EX@surv.each.timestep <-  matrix(NA,1,1)
  # Rescale the pop size (take the popuation distibution from state.t0 and multiply by POP) if you want to
  if (!is.null(tot.pop)) {
    EX@state.t0[,1] <- tot.pop*EX@state.t0[,1]/sum(EX@state.t0[,1])
  }
  # seed an infected and run out to equilibrium
  EX@state.t0[3+5*11,1] <- sum(EX@state.t0[,1])*0.00001 # seeding an infected person into i.inds because ind is 3rd epi class out of 5 therefore every 8 is an infected class - and want it intereted at age 11 or 88th row
  EX@state.t0[2+5*11,1] <- EX@state.t0[2+5*11,1]-sum(EX@state.t0[,1])*0.00001 # to keep population the same size we are pulling this 100 from susceptibles at age 11
  # This condition was initiating a bunch of vaccinated - change to recovered
  EX@state.t0[EX@trans@r.inds,1][EX@trans@age.class>60] <-EX@state.t0[EX@trans@s.inds,1][EX@trans@age.class>60]
  EX@state.t0[EX@trans@s.inds,1][EX@trans@age.class>60] <-0
  # Get rid of negatives
  EX@state.t0[which(EX@state.t0[,1]<0),1] <- 0
  # Run out transients - plot if you want to be sure...

  res <- run(EX, rescale.WAIFW=T)
  #plot(res@result)

  # Reset

  # function "GetNumber.per.AgeGroup()" in build.R
  age.struct.1991 <- GetNumber.per.AgeGroup(state=res@experiment.def@state.t0, trans=EX@trans)
  age.struc.sim <- GetNumber.per.AgeGroup(state=res@result@.Data[,ncol(res@result@.Data)], trans=EX@trans)
  prop.struc.sim <- res@result[,ncol(res@result)]/rep(age.struc.sim, each=5)
  new.state <- rep(age.struct.1991, each=5)*prop.struc.sim
  EX@state.t0[,1] <- new.state
  #EX@state.t0[,1] <- sum(res@result[,1])*res@result[,ncol(res@result)]/sum(res@result[,ncol(res@result)])
  EX@trans@waifw <- res@experiment.def@trans@waifw

  res <- run(EX, rescale.WAIFW=T)
  #plot(res@result)

  # Reset
  age.struct.1991 <- GetNumber.per.AgeGroup(state=res@experiment.def@state.t0, trans=EX@trans)
  age.struc.sim <- GetNumber.per.AgeGroup(state=res@result@.Data[,ncol(res@result@.Data)], trans=EX@trans)
  prop.struc.sim <- res@result[,ncol(res@result)]/rep(age.struc.sim, each=5)
  new.state <- rep(age.struct.1991, each=5)*prop.struc.sim
  EX@state.t0[,1] <- new.state
  #EX@state.t0[,1] <- sum(res@result[,1])*res@result[,ncol(res@result)]/sum(res@result[,ncol(res@result)])
  EX@trans@waifw <- res@experiment.def@trans@waifw

  res <- run(EX, rescale.WAIFW=T)
  #plot(res@result)

  # Reset
  age.struct.1991 <- GetNumber.per.AgeGroup(state=res@experiment.def@state.t0, trans=EX@trans)
  age.struc.sim <- GetNumber.per.AgeGroup(state=res@result@.Data[,ncol(res@result@.Data)], trans=EX@trans)
  prop.struc.sim <- res@result[,ncol(res@result)]/rep(age.struc.sim, each=5)
  new.state <- rep(age.struct.1991, each=5)*prop.struc.sim
  EX@state.t0[,1] <- new.state
  #EX@state.t0[,1] <- sum(res@result[,1])*res@result[,ncol(res@result)]/sum(res@result[,ncol(res@result)])
  EX@trans@waifw <- res@experiment.def@trans@waifw

  return(EX)
}
