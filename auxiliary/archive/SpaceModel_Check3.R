## 20 January 2023
## Attempting non-spatial run of subpopulation 1
## Looking for 3,3 non-zero in trans matrix
## Looking to see if running out the transients helps (i.e., different starting point)
## Looking to see if non-spatial model differs in other ways

load("/Users/winter/Library/CloudStorage/GoogleDrive-amykwinter@gmail.com/My Drive/Gavi_ZeroD_Personal/SchoolEntry/Subnational_Model_Sims/output/Zambia.Spatial.Setup.Object.RData")


year=2016
uncode=setup$uncode
iso3code="ZMB"
generation.time = 0.5
age.classes = c(1:240, seq(252,960,12))
maternal.decay.rt = 0.45
exponent = 0.97
frequency.dep=TRUE
tot.pop=NULL
yr.births.per.1000.acrossyears <- setup$cbr.1950.2100[1,year-1950+1]
intro.rate=1/24/320
targeted.intro=FALSE
R0=18
t.max = 20
get.births=setup$get.births.here
seasonal.amp = 0.15
flat.WAIFW=FALSE
country.specific.WAIFW=TRUE
vynnycky.waifw=FALSE
vynnycky.waifw.betas = c(2,1)
asdr.object <- new("nMx",
                   rate.years = setup$asdr.object@rate.years,
                   rates = setup$asdr.object@rates[1:101,],
                   mid.age = setup$asdr.object@mid.age)
use_montagu_demog=FALSE
routine.vac=0
routine.vac.age.index=12

## EX.Country.part1 ####

## Get.CountryX.Starting.Pop.MSIRV ####

## Calculate the aging rate using the age classes and the generation time
age.lows <- c(0,age.classes[2:length(age.classes)-1]) # makes it 0 to 697 rather than 1 to 709, so lowest age in wach age range
ac.sz <- age.classes - age.lows # the size in months of each age range (1 month till age 20, then 12 months to age 59)
aging.rate <- generation.time/ac.sz # converting generation time units into year/month time units based on size of age class
aging.rate[length(aging.rate)] <- 0 # forcing the last aging range to 0

## Returns the age profile of survivorship in units of the generation time
survs <- create.surv.prob.over.age(age.classes=age.classes,
                                   generation.time=generation.time,
                                   nMx=asdr.object, year=year)

## Create WAIFW matrix age.classes X age.classes dimensions, default is polymod basedo on great britain
if (flat.WAIFW) {
  waifw <- get.flat.WAIFW(age.classes/12)
}else{
  waifw <- get.polymod.WAIFW(age.classes/12)
}

## Set up maternal immunity
maternal.obj = new("maternal.exp.decay", decay.rt=maternal.decay.rt)

## Create the transition object
yr.births.per.1000 = yr.births.per.1000.acrossyears[1]
tran <- create.ID.transition.MSIRV(n.age.class = length(age.classes),
                                   aging.rate = aging.rate,
                                   survival.rate = survs,
                                   waifw = waifw,
                                   routine.vac=routine.vac,
                                   routine.vac.age.index=routine.vac.age.index,
                                   maternal.obj = maternal.obj,
                                   time.step = generation.time,
                                   age.class = age.classes,
                                   birth.rate =
                                     (yr.births.per.1000*tot.pop/1000)*
                                     (generation.time/12), #the expected number of total births for each time step
                                   exponent = exponent,
                                   frequency.dep=frequency.dep,
                                   is.stochastic = FALSE,
                                   get.births=get.births)


## space.create.country.x.DFE.ID.state.matrix ####
epi.class.label = c("M","S","I","R","V")

  #make template
  rc <- space.create.ID.state.matrix(tran@n.age.class,
                                     tran@n.epi.class,
                                     epi.class.label = epi.class.label,
                                     1)

  #fill in the template
  demog <- space.getDemography(tran@age.class, iso3code=iso3code)
  demog.n.ages <- nrow(demog$pop.age.byageclasses.1950.2100)/116

  #loop through the number of sub-populations
  for (s in 1){
    pop.struct <- demog$pop.age.byageclasses.1950.2100[(demog.n.ages*s-demog.n.ages+1):(demog.n.ages*s),(year-1950+1)]
    age <- as.numeric(tran@age.class)

    tot.pop <- demog$pop.total.1950.2100[s,(year-1950+1)]

    #Turn this into prop of desired age classes
    #note that NUMBER is imposed externally (since we might not want full countries... because interested in stochasticity)
    sp <- smooth.spline(age, pop.struct)  # gives "slope" between age (range 0-100) and population size per age (therefore 101 items)
    pred <- predict(sp, tran@age.class/12)$y # predicts population size per age class (total 280 age classes)
    #plot(age,pop.struct)
    #points(tran@age.class/12,pred,type="l",col=2)

    #Adjust for varying bin width
    pred <- pred*diff(c(0,tran@age.class))
    #added a 0 at the beginning of age classes
    #took the difference between each class n from n-1.   therefore difference is 1 (or 1 month) for item 1-241 and then 12 for items 242-281
    #then weighing each predicted population size per age class by the size of the age class

    #Fill in
    rc[tran@s.inds,1] <- tot.pop*pred/sum(pred)
    #if pop=NULL then assuming that the population w/in each state with missing ages are distributed within the population based on the pop structure that we do have ages for
    #filling in the susceptibles per age class
  }

state <- rc

## END of space.create.country.x.DFE.ID.state.matrix ####

## Update births now because need sum(state) = pop
## the expected number of total births for each time step
tran@birth.rate <- (yr.births.per.1000*sum(state)/1000)*generation.time/12

## Input transition introduction of infection rate
if (targeted.intro) {
  tran@introduction.rate <- rep(0, tran@n.age.class)
  tran@introduction.rate[60:71] <- intro.rate
} else {
  tran@introduction.rate <- rep(intro.rate, tran@n.age.class)
}

## Assumes everyone has maternal protection
state[tran@m.inds,1] <- state[tran@s.inds,1] * pmaternal(age.lows, maternal.obj)
state[tran@s.inds,1] <- state[tran@s.inds,1] * (1-pmaternal(age.lows, maternal.obj))

tmp <- list(state = state, tran = tran)

# END of Get.CountryX.Starting.Pop.MSIRV ####

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

# No rescaling of population -
EX@pop.rescale.each.timestep <- 0 #--xxamy - added this lines because when changed pop.rescale.each.timestep to "ANY" is became NULL by default but we need it to be a 0

waifw.prior.nospace <- EX@trans@waifw
# Run out transients - plot if you want to be sure...
res <- run(EX, rescale.WAIFW=T)
waifw.post.nospace <- res@experiment.def@trans@waifw
EX00 <- EX

#plot(res@result)
# Reset

# function "GetNumber.per.AgeGroup()" in build.R
#age.struct.1991 <- GetNumber.per.AgeGroup(state=res@experiment.def@state.t0, trans=EX@trans)
#age.struc.sim <- GetNumber.per.AgeGroup(state=res@result@.Data[,ncol(res@result@.Data)], trans=EX@trans)
#prop.struc.sim <- res@result[,ncol(res@result)]/rep(age.struc.sim, each=5)
#new.state <- rep(age.struct.1991, each=5)*prop.struc.sim
#EX@state.t0[,1] <- new.state
#EX@state.t0[,1] <- sum(res@result[,1])*res@result[,ncol(res@result)]/sum(res@result[,ncol(res@result)])
#EX@trans@waifw <- res@experiment.def@trans@waifw

#res <- run(EX, rescale.WAIFW=T)
#plot(res@result)

#load("/Users/winter/Library/CloudStorage/GoogleDrive-amykwinter@gmail.com/My Drive/Gavi_ZeroD_Personal/SchoolEntry/Subnational_Model_Sims/output/EXt0_hacked.RData")
load("~/Downloads/EX_space_part1.RData")

EX00@state.t0@.Data[1:1500,1] <- EX@state.t0@.Data[1:1500]
EX00@trans@waifw <- EX@trans@waifw

## END of EX.Country.part1 ####
EXt00 <- EX00
save(EXt00, file="~/Downloads/EXt00.RData")



## EX.Country.part2 ####
load("~/Downloads/EXt00.RData")
t.max = 30
year = 2016
uncode = setup$uncode
generation.time = 0.5
pop.rescale = as.numeric(setup$pop.total.1950.2100[1,seq(2026,(year+t.max-10),10)-1950+1])
pop.time = seq(2026,(year+t.max-10),10)-year
is.stochastic = FALSE
get.births = setup$get.births.here
t.max = t.max
rescale.WAIFW = F
yr.births.per.1000.acrossyears = setup$cbr.1950.2100[1,(year-1950+1):((year-1950+1)+t.max)]
asdr.object <- new("nMx",
                   rate.years = setup$asdr.object@rate.years,
                   rates = setup$asdr.object@rates[1:101,],
                   mid.age = setup$asdr.object@mid.age)
EXt0=EXt00
time.specific.MR1cov = setup$MCV1.coverage.1980to2100[1,(year-1980+1):(year-1980+t.max)]
time.specific.MR2cov = setup$MCV2.coverage.1980to2100[1,(year-1980+1):(year-1980+t.max)]
time.specific.SIAcov = setup$measlesSIA.coverage.1980to2100[1,(year-1980+1):(year-1980+t.max)]
time.specific.min.age.MR1 = setup$time.specific.min.age.MR1[(year-1980+1):(year-1980+t.max)]
time.specific.max.age.MR1 = setup$time.specific.max.age.MR1[(year-1980+1):(year-1980+t.max)]
time.specific.min.age.MR2 = setup$time.specific.min.age.MR2[(year-1980+1):(year-1980+t.max)]
time.specific.max.age.MR2 = setup$time.specific.max.age.MR2[(year-1980+1):(year-1980+t.max)]
time.specific.min.age.SIA = setup$time.specific.min.age.SIA[(year-1980+1):(year-1980+t.max)]
time.specific.max.age.SIA = setup$time.specific.max.age.SIA[(year-1980+1):(year-1980+t.max)]
obj.vcdf.MR1 = get.MR1cdf.survival(uncode=setup$uncode, 1, 24)
obj.vcdf.MR2 = get.vcdf.normal(25, 36)
obj.prob.vsucc = pvacsuccess(1:(20*12), get.boulianne.vsucc())
sia.timing.in.year = (10/12)
MR1MR2correlation = FALSE
MR1SIAcorrelation = FALSE
MR2SIAcorrelation = FALSE
SIAinacc = FALSE
prop.inacc = NULL
SIAinefficient = FALSE


## Changing experiment type
# call to new returns a newly allocated object from the class identified by first argument
if (!MR1MR2correlation) EX <- new("experiment.updatedemog.vaccinationchange") #default
if (MR1MR2correlation) EX <- new("experiment.updatedemog.vaccinationchange.vaccinationlimitations")

## Country name
name <- countrycode::countrycode(uncode, origin="un", destination="country.name")

## Update new experiment
# EX is an experiment object
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
EX@pop.rescale.each.timestep <- rep(NaN,no.time.steps.in.experiment)
if (!is.null(pop.rescale)) { #if pop.rescale is not NULL fill in the pops to rescale at mid-year
  for (p in 1:length(pop.rescale)){
    EX@pop.rescale.each.timestep[pop.time[p]*no.gens.in.year-1] <- pop.rescale[p] #rescalling at end of previous year xxamy change from poc
  }
}

# Setting up changing birth rates
EX@births.per.1000.each.timestep <- c(as.numeric(yr.births.per.1000.acrossyears[1]),
                                      rep(as.numeric(yr.births.per.1000.acrossyears), each=no.gens.in.year))*generation.time/12

# Setting up survival
EX@surv.each.timestep <- create.surv.prob.over.age.time(EX@trans@age.class, generation.time, nMx=asdr.object, nMx.years=seq(year,(year+t.max-1),1)) #xxamny minus 1 to get correct number of tie steps

# Setting up vaccination pieces
EX@time.specific.MR1cov =  time.specific.MR1cov #vector of values annual coverage values 1980-2100
EX@time.specific.MR2cov = time.specific.MR2cov #vector of values annual coverage values 1980-2100
EX@time.specific.SIAcov =   time.specific.SIAcov #vector of values annual coverage values 1980-2100
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
if (MR1MR2correlation | MR1SIAcorrelation | MR2SIAcorrelation | SIAinefficient | SIAinacc){
  EX@MR1MR2correlation = MR1MR2correlation
  EX@MR1SIAcorrelation = MR1SIAcorrelation
  EX@MR2SIAcorrelation = MR2SIAcorrelation
  EX@SIAinefficient = SIAinefficient
  EX@SIAinacc = SIAinacc
  if (SIAinacc) {
    #Setting up population always inaccessible to vaccine
    EX@prop.inacc <- c(rep(prop.inacc, each=no.gens.in.year),  prop.inacc[length(prop.inacc)]) #by time step
    #EX@prop.inacc = prop.inacc #by year
  }
}
MR1MR2correlation=TRUE

#return(run(EX, rescale.WAIFW=rescale.WAIFW))
exper <- EX

