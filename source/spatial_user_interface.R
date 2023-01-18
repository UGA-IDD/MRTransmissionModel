## Wrapper functions to run SPATIAL experiments
## Authors: Amy Winter


#Function to set up SPATIAL experiment on disease free population
#
#Parameters -
#     uncode - UN country code
#     generation.time - the desired generation time in months
#     age.classes - vector - the upper limit of the age classes we are interested in months
#     maternal.decay.rt -rate of maternal decay of immunity to rubella
#     exponent - numeric - exponent for the infected
#     tot.subpop - numeric - population you want to scale you whole experiment by at the beginning, NULL or numeric
#     yr.births.per.1000.bysubpop - crude birth rate per 1000 for each subpopulation
#     coverage - proportion - vaccination coverage
#     targeted.intro - boolean - FALSE means space introduction out over all age classes, TRUE means to concentrate introductions
#     intro.rate - numeric - rate at which new infections are introduced
#     R0 - numeric - basic reproductive number of rubella assumed, used to rescale the waifw
#     t.max - numeric - time in year that plan to run out transients
#     get.births - vector of births with length of 5*n.age.classes
#     seasonal.amp - numeric - amplitude of the seasonal cosine function
#     flat.WAIFW - boolean - if false then default waifw is polymod GB
#     country.specific.WAIFW - boolean - if true then overrides polymod or flat waifw
#     vynnycky.waifw - boolean - true if want emilia's waifw <13 and 15-50 (beta_young, beta_old, and 0.7*beta_old)
#     vynnycky.waifw.betas - vector of two betas (beta_young and beta_old)
#     space.asdr.object - nMx object with rate and mid-age (age specific death rates)
#     year - numeric - year to pull DFE demography (population and age structure)
#     n.subpops - number of subpopulations
#     subpop.class.numeric.label - numeric label xxxamy 22dec2022
#     subpop.class.character.label - character label associted with the numeric label xxxamy 22dec2022
#     coupling - matrix of connectivity between subpopulations
#
#Returns -
#     tran object and state.t0 after transients run out
Spatial.EX.DiseaseFree.part1 <- function(uncode,
                                     generation.time = 0.5, #generation time in months
                                     age.classes = c(1:240, seq(241,720,12)),
                                     maternal.decay.rt=0.95, #based on Metcalf, 2012 paper
                                     exponent = 0.97,
                                     frequency.dep=T,
                                     tot.subpop=NULL,
                                     yr.births.per.1000.bysubpop,
                                     intro.rate=1e-10,
                                     targeted.intro=F,
                                     R0=5,
                                     t.max = 20,
                                     get.births=NULL,
                                     seasonal.amp=0.35, #metcalf et al 2012 used 0.35
                                     flat.WAIFW=F,
                                     country.specific.WAIFW=F,
                                     vynnycky.waifw=F,
                                     vynnycky.waifw.betas = c(2,1),
                                     space.asdr.object,
                                     year=1990,
                                     routine.vac=0,
                                     routine.vac.age.index=12,
                                     n.subpops=2,
                                     subpop.class.numeric.label, #xxxamy 22dec2022
                                     subpop.class.character.label, #xxxamy 22dec2022
                                     coupling = matrix(1, nrow=n.subpops, ncol=n.subpops)) { 
  
  tmp <- Space.Get.CountryX.Starting.Pop.MSIRV(uncode=uncode, 
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
                                               subpop.class.label = subpop.class.numeric.label, #xxxamy 22dec2022
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
  EX@subpop.class.character.label <- subpop.class.character.label #xxxamy 22dec2022
  
  #change the WAIFW to country-specific if requested (based on prem et al.)
  if (country.specific.WAIFW) EX@trans@waifw <- get.prem.WAIFW(age.class.boundries=age.classes/12, 
                                                               uncode, other.contact.matrix=F,
                                                               bandwidth=c(3,3), 
                                                               adjustment_start_year=F, year=1980)
  
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
  
  # Rescale the WAIFW if a specific R0 is specified
  EX@trans@waifw <- scale.WAIFW.space(EX@R0,
                                     EX@state.t0,EX@trans@waifw,
                                     frequency.dep=EX@trans@frequency.dep,
                                     suscept.state=EX@trans@s.inds[1])
  
  return(EX)
  
}


#Function to set up SPATIAL experiment and run the transients out
#
#Parameters -
#     uncode - UN country code
#     generation.time - the desired generation time in months
#     age.classes - vector - the upper limit of the age classes we are interested in months
#     maternal.decay.rt -rate of maternal decay of immunity to rubella
#     exponent - numeric - exponent for the infected
#     tot.subpop - numeric - population you want to scale you whole experiment by at the beginning, NULL or numeric
#     yr.births.per.1000.bysubpop - crude birth rate per 1000 for each subpopulation
#     coverage - proportion - vaccination coverage
#     targeted.intro - boolean - FALSE means space introduction out over all age classes, TRUE means to concentrate introductions
#     intro.rate - numeric - rate at which new infections are introduced
#     R0 - numeric - basic reproductive number of rubella assumed, used to rescale the waifw
#     t.max - numeric - time in year that plan to run out transients
#     get.births - vector of births with length of 5*n.age.classes
#     seasonal.amp - numeric - amplitude of the seasonal cosine function
#     flat.WAIFW - boolean - if false then default waifw is polymod GB
#     country.specific.WAIFW - boolean - if true then overrides polymod or flat waifw
#     vynnycky.waifw - boolean - true if want emilia's waifw <13 and 15-50 (beta_young, beta_old, and 0.7*beta_old)
#     vynnycky.waifw.betas - vector of two betas (beta_young and beta_old)
#     space.asdr.object - nMx object with rate and mid-age (age specific death rates)
#     year - numeric - year to pull DFE demography (population and age structure)
#     n.subpops - number of subpopulations
#     subpop.class.numeric.label - numeric label xxxamy 22dec2022
#     subpop.class.character.label - character label associted with the numeric label xxxamy 22dec2022
#     coupling - matrix of connectivity between subpopulations
#
#Returns -
#     tran object and state.t0 after transients run out
Spatial.EX.Country.part1 <- function(uncode,
                             generation.time = 0.5, #generation time in months
                             age.classes = c(1:240, seq(241,720,12)),
                             maternal.decay.rt=0.95, #based on Metcalf, 2012 paper
                             exponent = 0.97,
                             frequency.dep=T,
                             tot.subpop=NULL,
                             yr.births.per.1000.bysubpop,
                             intro.rate=1e-10,
                             targeted.intro=F,
                             R0=5,
                             t.max = 20,
                             get.births=NULL,
                             seasonal.amp=0.35, #metcalf et al 2012 used 0.35
                             flat.WAIFW=F,
                             country.specific.WAIFW=F,
                             vynnycky.waifw=F,
                             vynnycky.waifw.betas = c(2,1),
                             space.asdr.object,
                             year=1990,
                             routine.vac=0,
                             routine.vac.age.index=12,
                             n.subpops=2,
                             subpop.class.numeric.label, #xxxamy 22dec2022
                             subpop.class.character.label, #xxxamy 22dec2022
                             coupling = matrix(1, nrow=n.subpops, ncol=n.subpops)) { 
  
  tmp <- Space.Get.CountryX.Starting.Pop.MSIRV(uncode=uncode, 
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
                                               subpop.class.label = subpop.class.numeric.label, #xxxamy 22dec2022
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
  EX@subpop.class.character.label <- subpop.class.character.label #xxxamy 22dec2022
  
  #change the WAIFW to country-specific if requested (based on prem et al.)
  if (country.specific.WAIFW) EX@trans@waifw <- get.prem.WAIFW(age.class.boundries=age.classes/12, 
                                                               uncode, other.contact.matrix=F,
                                                               bandwidth=c(3,3), 
                                                               adjustment_start_year=F, year=1980)
  
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
  age.struct.1991 <- space.wrapper.GetNumber.per.AgeGroup(state=res@experiment.def@state.t0, trans=EX@trans)
  age.struc.sim <- space.wrapper.GetNumber.per.AgeGroup(state=matrix(res@result@.Data[,ncol(res@result@.Data)],ncol=1), trans=EX@trans)
  for (s in 1:n.subpops){
    prop.struc.sim <- res@result[res@experiment.def@trans@subpop.class.label==s,ncol(res@result)]/rep(age.struc.sim[s,], each=5)
    tmp <- rep(age.struct.1991[s,], each=5)*prop.struc.sim
    if (s==1) new.state <- tmp
    if (s!=1) new.state <- c(new.state, tmp)
  }
  EX@state.t0[,1] <- new.state
  EX@trans@waifw <- res@experiment.def@trans@waifw
  
  res <- run(EX, rescale.WAIFW=T)
  #plot(res@result)
  
  # Reset
  age.struct.1991 <- space.wrapper.GetNumber.per.AgeGroup(state=res@experiment.def@state.t0, trans=EX@trans)
  age.struc.sim <- space.wrapper.GetNumber.per.AgeGroup(state=matrix(res@result@.Data[,ncol(res@result@.Data)],ncol=1), trans=EX@trans)
  for (s in 1:n.subpops){
    prop.struc.sim <- res@result[res@experiment.def@trans@subpop.class.label==s,ncol(res@result)]/rep(age.struc.sim[s,], each=5)
    tmp <- rep(age.struct.1991[s,], each=5)*prop.struc.sim
    if (s==1) new.state <- tmp
    if (s!=1) new.state <- c(new.state, tmp)
  }
  EX@state.t0[,1] <- new.state
  EX@trans@waifw <- res@experiment.def@trans@waifw
  
  res <- run(EX, rescale.WAIFW=T)
  plot(res@result)
  
  # Reset
  age.struct.1991 <- space.wrapper.GetNumber.per.AgeGroup(state=res@experiment.def@state.t0, trans=EX@trans)
  age.struc.sim <- space.wrapper.GetNumber.per.AgeGroup(state=matrix(res@result@.Data[,ncol(res@result@.Data)],ncol=1), trans=EX@trans)
  for (s in 1:n.subpops){
    prop.struc.sim <- res@result[res@experiment.def@trans@subpop.class.label==s,ncol(res@result)]/rep(age.struc.sim[s,], each=5)
    tmp <- rep(age.struct.1991[s,], each=5)*prop.struc.sim
    if (s==1) new.state <- tmp
    if (s!=1) new.state <- c(new.state, tmp)
  }
  EX@state.t0[,1] <- new.state
  EX@trans@waifw <- res@experiment.def@trans@waifw
  
  return(EX)
  
}


#Function to Run the SPATIAL Simulations/Experiments
#
#Parameters -
#     uncode - UN country code
#     generation.time - generation time in months
#     coverage.sia - numeric - proportion of children that will be reached by SIA
#     age.range.sia - numeric - age range of the SIA in months
#     sia.biweek - numeric - the time the SIAs are to begin in biweeks
#     pop.rescale - numeric - population by which you want to rescale at pop.time 
#     pop.time - numeric - time in YEARS that you want to rescale the population by
#     is.stochastic - logical
#     get.births - vector of births with length of 5*n.age.classes
#     t.max - numeric - time in year that plan to run the experiment
#     RCV1.max.age - numeric - age by which RCV1 vaccination is achieved in months
#     RCV1.vaccination - vector - RCV1 coverage over time, length(RCV1.vaccination) must match length(RCV2.vaccination)
#     RCV2.age.range - vector - min and max age or RCV2, in months
#     RCV2.vaccination - vector - RCV2 coverage over time, length(RCV2.vaccination) must match length(RCV1.vaccination)
#     time.specific.routine - matrix - proportion routine vaccination, nrow=length(setup$RCV1.coverage.1980to2100) ncol=length(EX@trans@age.class))
#     index.routine.to.time.steps - vector - time steps corresponding the RCV1.vaccination and RCV2.vaccination
#     prop.vacc.inaccessible - numeric - prop. inaccessible to vacccination  (but can revise to be per year)
#     rescale.WAIFW - logical  
#     yr.births.per.1000.acrossyears.bysupop - vector of length(t.max) - crude birth rate per 1000 per year
#     space.asdr.object - space nMx object - with rate and mid-age (age specific death rates)
#     year - numeric - year that pulled DFE demography (population and age structure) from part1
#     EX - experiment object - based on part1()
#
#Returns -
#     experiment results
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
                             time.specific.min.age.MR1,
                             time.specific.max.age.MR1,
                             time.specific.min.age.MR2,
                             time.specific.max.age.MR2,
                             time.specific.min.age.SIA,
                             time.specific.max.age.SIA,
                             obj.vcdf.MR1 = get.vcdf.normal(6, 12),
                             obj.vcdf.MR2 = get.vcdf.normal(15, 21),
                             obj.prob.vsucc = pvacsuccess(1:(20*12), get.boulianne.vsucc()),
                             sia.timing.in.year = (3/12),
                             MR1MR2correlation = NULL,
                             MR1SIAcorrelation = NULL,
                             MR2SIAcorrelation = NULL,
                             SIAinacc = NULL,
                             prop.inacc = NULL,
                             SIAinefficient = NULL){   
  
  
  ## Changing experiment type
  if (is.null(MR1MR2correlation)) EX <- new("experiment.updatedemog.vaccinationchange.spatial") #default
  if (!is.null(MR1MR2correlation)) EX <- new("experiment.updatedemog.vaccinationchange.vaccinationlimitations.spatial")
  
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
      EX@pop.rescale.each.timestep[,pop.time[p]*no.gens.in.year-1] <- pop.rescale[,p] #rescalling at end of previous year
    }
  }
  
  # Setting up changing birth rates
  EX@births.per.1000.each.timestep <- cbind(yr.births.per.1000.acrossyears.bysupop[,1], 
                                            matrix(rep(yr.births.per.1000.acrossyears.bysupop, each = no.gens.in.year), byrow=T, nrow=EX@trans@n.subpops))*generation.time/12
  
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
  
  return(run(EX, rescale.WAIFW=rescale.WAIFW))
  
}

