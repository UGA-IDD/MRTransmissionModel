#' Function to set up SPATIAL experiment and run the transients out
#'
#' @param uncode UN country code
#' @param iso3code character iso3 code to matches uncode
#' @param generation.time the desired generation time in months
#' @param age.classes vector - the upper limit of the age classes we are interested in months
#' @param maternal.decay.rt rate of maternal decay of immunity to rubella
#' @param exponent exponent for the infected
#' @param frequency.dep xxx
#' @param tot.subpop population you want to scale you whole experiment by at the beginning, NULL or numeric
#' @param yr.births.per.1000.bysubpop crude birth rate per 1000 for each subpopulation
#' @param intro.rate numeric - rate at which new infections are introduced
#' @param targeted.intro boolean - FALSE means space introduction out over all age classes, TRUE means to concentrate introductions
#' @param R0 numeric - basic reproductive number of rubella assumed, used to rescale the waifw
#' @param t.max numeric - time in year that plan to run out transients
#' @param get.births vector of births with length of 5*n.age.classes
#' @param seasonal.amp numeric - amplitude of the seasonal cosine function
#' @param flat.WAIFW boolean - if false then default waifw is polymod GB
#' @param country.specific.WAIFW boolean - if true then overrides polymod or flat waifw
#' @param vynnycky.waifw boolean - true if want emilia's waifw <13 and 15-50 (beta_young, beta_old, and 0.7*beta_old)
#' @param vynnycky.waifw.betas vector of two betas (beta_young and beta_old)
#' @param space.asdr.object nMx object with rate and mid-age (age specific death rates)
#' @param year numeric - year to pull DFE demography (population and age structure)
#' @param routine.vac xxx
#' @param routine.vac.age.index xxx
#' @param n.subpops number of subpopulations
#' @param coupling matrix of connectivity between subpopulations
#' @param starting.prop.immune vector - starting proportion immune, length is length(starting.prop.immune.ages.in.months)*n.subpops
#' @param starting.prop.immune.ages.in.months vector - the ages (in months) associated with starting.prop.immune
#'
#' @return tran object and state.t0 after transients run out
#' @export
#'

Spatial.EX.Country.part1 <- function(uncode,
                                     iso3code,
                                     generation.time = 0.5, #generation time in months
                                     age.classes = c(1:240, seq(241,720,12)),
                                     maternal.decay.rt=0.95, #based on Metcalf, 2012 paper
                                     exponent = 0.97,
                                     frequency.dep=TRUE,
                                     tot.subpop=NULL,
                                     yr.births.per.1000.bysubpop,
                                     intro.rate=1e-10,
                                     targeted.intro=FALSE,
                                     R0=5,
                                     t.max = 20,
                                     get.births=NULL,
                                     seasonal.amp=0.35, #metcalf et al 2012 used 0.35
                                     flat.WAIFW=FALSE,
                                     country.specific.WAIFW=FALSE,
                                     vynnycky.waifw=FALSE,
                                     vynnycky.waifw.betas = c(2,1),
                                     space.asdr.object,
                                     year=1990,
                                     routine.vac=0,
                                     routine.vac.age.index=12,
                                     n.subpops=2,
                                     coupling = matrix(1, nrow=n.subpops, ncol=n.subpops),
                                     starting.prop.immune = c(0.75, 0.75, 0.75, 0.75),
                                     starting.prop.immune.ages.in.months = c(9, 12, 24, 240),
                                     max.immunity = 0.98) {

  tmp <- Space.Get.CountryX.Starting.Pop.MSIRV(uncode=uncode,
                                               iso3code=iso3code,
                                               generation.time=generation.time,
                                               age.classes=age.classes,
                                               maternal.decay.rt=maternal.decay.rt,
                                               exponent=exponent,
                                               frequency.dep = frequency.dep,
                                               is.stochastic = FALSE, # Make NOT stochastic for the transient ridding stage
                                               tot.subpop=tot.subpop,
                                               yr.births.per.1000.bysubpop = yr.births.per.1000.bysubpop,
                                               intro.rate=intro.rate,
                                               flat.WAIFW=flat.WAIFW,
                                               space.asdr.object=space.asdr.object,
                                               targeted.intro=targeted.intro,
                                               year=year,
                                               get.births=get.births,
                                               routine.vac=routine.vac,
                                               routine.vac.age.index=routine.vac.age.index,
                                               n.subpops=n.subpops,
                                               coupling=coupling)


  EX <- new("experiment.updatedemog.spatial")
  EX@name <- paste(uncode, "Frequency Dependent", sep=" ")
  EX@state.t0 = tmp$state
  EX@trans = tmp$tran
  EX@t.min = 0
  EX@t0.doy = 0
  EX@step.size = 1/(12*1/generation.time) #1/no.gens.per.year
  EX@season.obj = new("seasonal.cosine", amplitude = seasonal.amp)
  EX@R0 <- R0
  EX@t.max <-t.max
  EX@maternal.obj <- new("maternal.exp.decay", decay.rt=maternal.decay.rt)

  #change the WAIFW to country-specific if requested (based on prem et al.)
  if (country.specific.WAIFW) EX@trans@waifw <- get.prem.WAIFW(age.class.boundries=age.classes/12,
                                                               uncode, other.contact.matrix=FALSE,
                                                               bandwidth=c(3,3),
                                                               adjustment_start_year=FALSE, year=1980)

  if (vynnycky.waifw) EX@trans@waifw <- get.vynnycky.WAIFW(age.classes/12,
                                                           beta_young=vynnycky.waifw.betas[1],
                                                           beta_old=vynnycky.waifw.betas[2])

  # To run out the transients, keep the birth rate constant across time
  EX@births.per.1000.each.timestep <- matrix(rep(yr.births.per.1000.bysubpop, each = (t.max+1)*12*1/generation.time), byrow=T, nrow=n.subpops)*generation.time/12

  # Setting up survival over time as NA in order to use the original surv.matrix set up for trans object and keep constant over time
  EX@surv.each.timestep <-  matrix(NA,1,1)

  # Re-scale the pop size (take the population distribution from state.t0 and multiply by POP) if you want to
  if (!is.null(tot.subpop)) {
    for (s in 1:n.subpops){
      EX@state.t0[EX@trans@subpop.class.label==s,1] <- tot.subpop[s]*EX@state.t0[EX@trans@subpop.class.label==s,1]/sum(EX@state.t0[EX@trans@subpop.class.label==s,1])
    }
  }

  for (s in 1:n.subpops){
    # seed an infected and run out to equilibrium
    EX@state.t0[which(EX@trans@subpop.class.label==s)[(3+5*60)],1] <- sum(EX@state.t0[EX@trans@subpop.class.label==s,1])*0.00001 # seeding infected 60 month olds into i.inds because ind is 3rd epi class out of 5
    EX@state.t0[which(EX@trans@subpop.class.label==s)[(2+5*60)],1] <- EX@state.t0[which(EX@trans@subpop.class.label==s)[(2+5*60)],1]-sum(EX@state.t0[EX@trans@subpop.class.label==s,1])*0.00001 # to keep population the same size we are pulling this 100 from susceptibles at age 60 months old

    # This condition was initiating a bunch of vaccinated - change to recovered
    index.loc <-  EX@trans@subpop.class.label[EX@trans@i.inds]
    this.loc <- index.loc==s
    EX@state.t0[EX@trans@r.inds[this.loc],1][EX@trans@age.class>60] <- EX@state.t0[EX@trans@s.inds[this.loc],1][EX@trans@age.class>60]
    EX@state.t0[EX@trans@s.inds[this.loc],1][EX@trans@age.class>60] <-0
  }

  # Get rid of negatives
  EX@state.t0[which(EX@state.t0[,1]<0),1] <- 0

  # Don't rescale the population
  no.time.steps.in.experiment <- round((EX@t.max-EX@t.min)/EX@step.size)+1
  EX@pop.rescale.each.timestep <- matrix(NaN, EX@trans@n.subpops, no.time.steps.in.experiment)

  # Run out transients - plot if you want to be sure...
  res <- run(EX, rescale.WAIFW=T)
  #plot(res@result)

  # Reset
  age.struc.t0 <- space.wrapper.GetNumber.per.AgeGroup(state=res@experiment.def@state.t0, trans=EX@trans)
  age.struc.tmax <- space.wrapper.GetNumber.per.AgeGroup(state=matrix(res@result@.Data[,ncol(res@result@.Data)],ncol=1), trans=EX@trans)
  for (s in 1:n.subpops){
    prop.struc.tmax <- res@result[res@experiment.def@trans@subpop.class.label==s,ncol(res@result)]/rep(age.struc.tmax[s,], each=5)
    tmp <- rep(age.struc.t0[s,], each=5)*prop.struc.tmax
    if (s==1) new.state <- tmp
    if (s!=1) new.state <- c(new.state, tmp)
  }
  EX@state.t0[,1] <- new.state
  EX@trans@waifw <- res@experiment.def@trans@waifw

  res <- run(EX, rescale.WAIFW=T)
  #plot(res@result)

  # Reset
  age.struc.t0 <- space.wrapper.GetNumber.per.AgeGroup(state=res@experiment.def@state.t0, trans=EX@trans)
  age.struc.tmax <- space.wrapper.GetNumber.per.AgeGroup(state=matrix(res@result@.Data[,ncol(res@result@.Data)],ncol=1), trans=EX@trans)
  for (s in 1:n.subpops){
    prop.struc.tmax <- res@result[res@experiment.def@trans@subpop.class.label==s,ncol(res@result)]/rep(age.struc.tmax[s,], each=5)
    tmp <- rep(age.struc.t0[s,], each=5)*prop.struc.tmax
    if (s==1) new.state <- tmp
    if (s!=1) new.state <- c(new.state, tmp)
  }
  EX@state.t0[,1] <- new.state
  EX@trans@waifw <- res@experiment.def@trans@waifw

  res <- run(EX, rescale.WAIFW=T)
  #plot(res@result)

  # Fix the starting susceptible population
  #reset based on tmax
  age.struc.t0 <- space.wrapper.GetNumber.per.AgeGroup(state=res@experiment.def@state.t0, trans=EX@trans)
  age.struc.tmax <- space.wrapper.GetNumber.per.AgeGroup(state=res@result@.Data[,ncol(res@result@.Data)], trans=EX@trans)
  for (s in 1:n.subpops){
    prop.struc.tmax <- res@result[res@experiment.def@trans@subpop.class.label==s,ncol(res@result)]/rep(age.struc.tmax[s,], each=5)
    tmp <- rep(age.struc.t0[s,], each=5)*prop.struc.tmax
    if (s==1) new.state <- tmp
    if (s!=1) new.state <- c(new.state, tmp)
  }
  #determine which age indexes match the argument starting.prop.sus.ages.in.months
  indexes <- which(findInterval(age.classes, starting.prop.immune.ages.in.months[1])==1 &
                     findInterval(age.classes, starting.prop.immune.ages.in.months[length(starting.prop.immune.ages.in.months)])==0)
  #interpolate ages based on age.classes for each subpop
  for (s in 1:n.subpops){
    f <- smooth.spline(starting.prop.immune.ages.in.months, starting.prop.immune[,s])
    pred.prop.immune <- predict(f,age.classes[indexes])$y
    pred.prop.immune <- sapply(pred.prop.immune, function(x) min(x, max.immunity))
    #replace age-specific profiles IF proportion susceptible given for that age group
    new.state[EX@trans@r.inds[indexes]] <- age.struc.t0[s,indexes]*pred.prop.immune
    new.state[EX@trans@s.inds[indexes]] <- age.struc.t0[s,indexes]*(1-pred.prop.immune)
    new.state[EX@trans@m.inds[indexes]] <- 0
    new.state[EX@trans@i.inds[indexes]] <- 0
    new.state[EX@trans@v.inds[indexes]] <- 0
  }

  EX@state.t0[,1] <- new.state

  return(EX)

}
