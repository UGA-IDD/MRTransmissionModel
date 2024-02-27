## 18 Janurary 2024
## Making sure space model and non-space model (which uses c code) are doing the same thing
## As far as I can tell the production of the trans matrix and multiplication to the state vector
## is exactly the same

## Will look further at introduction rates


load("/Users/winter/Library/CloudStorage/GoogleDrive-amykwinter@gmail.com/My Drive/Gavi_ZeroD_Personal/SchoolEntry/Subnational_Model_Sims/output/Zambia.Spatial.Setup.Object.RData")
load("/Users/winter/Library/CloudStorage/GoogleDrive-amykwinter@gmail.com/My Drive/Gavi_ZeroD_Personal/SchoolEntry/Subnational_Model_Sims/output/EXt0_hacked.RData")

## Space ####
year = 2016
t.max = 2
uncode = setup$uncode
generation.time = 0.5
pop.rescale = NULL
is.stochastic = FALSE
get.births = setup$get.births.here
rescale.WAIFW = F
yr.births.per.1000.acrossyears.bysupop = setup$cbr.1950.2100[,(year-1950+1):((year-1950+1)+t.max)]
space.asdr.object = setup$asdr.object
time.specific.MR1cov = setup$MCV1.coverage.1980to2100[,(year-1980+1):(year-1980+t.max)]
time.specific.MR2cov = setup$MCV2.coverage.1980to2100[,(year-1980+1):(year-1980+t.max)]
time.specific.SIAcov = setup$measlesSIA.coverage.1980to2100[,(year-1980+1):(year-1980+t.max)]
time.specific.schoolvacc.cov=NULL
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
schoolvacc.timing.in.year = (9/12)
MR1MR2correlation=NULL

#part2 ####
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

#return(run(EX, rescale.WAIFW=rescale.WAIFW))
exper <- EX

##run####
#print(rescale.WAIFW)
state <- exper@state.t0

#number of sub-populations
n.subpops <- exper@trans@n.subpops

#rescale the WAIFW if a specific R0 is specified
if (rescale.WAIFW & length(exper@R0)>0) {
  #print("RESCALING!!!")
  exper@trans@waifw <- scaleWAIFW.space(exper@R0,
                                        state,exper@trans@waifw,
                                        frequency.dep=exper@trans@frequency.dep,
                                        suscept.state=exper@trans@s.inds[1])
}

#get the number of time steps in the experiment
T <- round((exper@t.max-exper@t.min)/exper@step.size)+1

#get seasonal mult
if (!is.null(exper@season.obj)) {
  mults <- get.seasonal.mult(exper@t0.doy/365+(1:T-1)*exper@step.size,
                             exper@season.obj)
} else {
  mults <- rep(1,T)
}

#make a temporary transmission object
tmp.trans <- exper@trans

#generate the age and space (rows) and year (columns) specific vaccination matrix
routine <- space.wrapper.get.routine.time.age.specific(time.step= exper@step.size*12,
                                                       age.classes=exper@trans@age.class,
                                                       space.time.specific.MR1cov=exper@time.specific.MR1cov,
                                                       age.min.MR1=exper@time.specific.min.age.MR1,
                                                       age.max.MR1=exper@time.specific.max.age.MR1,
                                                       space.time.specific.MR2cov=exper@time.specific.MR2cov,
                                                       age.min.MR2=exper@time.specific.min.age.MR2,
                                                       age.max.MR2=exper@time.specific.max.age.MR2,
                                                       obj.vcdf.MR1=exper@obj.vcdf.MR1,
                                                       obj.vcdf.MR2=exper@obj.vcdf.MR2,
                                                       obj.prob.vsucc=exper@obj.prob.vsucc,
                                                       MR1MR2correlation=T)

#getting year of routine introductions
routine.intro <- rep(0, T)
if (any(colSums(exper@time.specific.MR1cov)!=0)) routine.intro[min(which(colSums(exper@time.specific.MR1cov)>0))*(1/exper@step.size)+1] <- 1
if (any(colSums(exper@time.specific.MR2cov)!=0)) routine.intro[min(which(colSums(exper@time.specific.MR2cov)>0))*(1/exper@step.size)+1] <- 1

#specifying that each year (column) of routine should be repeated 24 times
index.routine.vacc <- c(1,rep(1:ncol(routine$age.time.specific.routine), each=(T-1)/exper@t.max))

#generate the age and space (rows) and year (columns) specific vaccination matrix - same size as `routine`
SIA <- space.wrapper.get.sia.time.age.specific(age.classes=exper@trans@age.class,
                                               space.time.specific.SIAcov=exper@time.specific.SIAcov,
                                               age.min.sia=exper@time.specific.min.age.SIA,
                                               age.max.sia=exper@time.specific.max.age.SIA,
                                               obj.prob.vsucc=exper@obj.prob.vsucc)

#getting sia.times as a vector of 0 / 1 to represent when the SIA is to take place
index.sia.vacc <- rep(NA,T)
year.sia <- which(colSums(exper@time.specific.SIAcov)!=0)
index.sia.vacc[(year.sia-1)*(T-1)/exper@t.max + round((exper@sia.timing.in.year*(T-1)/exper@t.max))] <-  year.sia #minus 1 because adding the sia.timing
sia.times <- ifelse(!is.na(index.sia.vacc), 1, 0)

#need output vectors for primary vaccination failure over time
MR1.fail.each.timestep <- MR2.fail.each.timestep <- SIA.fail.each.timestep <- matrix(0, n.subpops, T)
MR1.fail.each.timestep[,1] <- routine$prop.fail.MR1[,1]
MR2.fail.each.timestep[,1] <- routine$prop.fail.MR2[,1]
SIA.fail.each.timestep[,1] <- 0

t=2

print(t)

#scale the waifw by seasonality
if (!is.array(mults)){
  tmp.trans@waifw <- exper@trans@waifw*mults[t]
} else {
  tmp.trans@waifw <- exper@trans@waifw*mults[,,t]
}

#if exper@pop.rescale.each.timestep[t]==NaN then tmp.trans@pop.rescale==NaN and will not rescale
#otherwise if any number other than NaN it will rescale
#or if exper@pop.rescale.each.timestep not completed , then numeric(0) and any index [t] is NA
for (s in 1:n.subpops) {
  if (!is.na(exper@pop.rescale.each.timestep[s,t])) {
    last.index.tmp <- (s*exper@trans@n.epi.class*exper@trans@n.age.class)
    first.index.tmp <- last.index.tmp-(exper@trans@n.epi.class*exper@trans@n.age.class)+1
    state[first.index.tmp:last.index.tmp] <- exper@pop.rescale.each.timestep[s,t]*
      (state[first.index.tmp:last.index.tmp]/sum(state[first.index.tmp:last.index.tmp]))
  }
}

#put in the correct birth rate for that time-step, if it varies
for (s in 1:n.subpops) {
  if (ncol(exper@births.per.1000.each.timestep)>1) {
    last.index.tmp <- (s*exper@trans@n.epi.class*exper@trans@n.age.class)
    first.index.tmp <- last.index.tmp-(exper@trans@n.epi.class*exper@trans@n.age.class)+1
    tmp.trans@birth.rate[s] = (exper@births.per.1000.each.timestep[s,t]*
                                 sum(state[first.index.tmp:last.index.tmp])/1000)
  } else {
    tmp.trans@birth.rate[s] <- exper@trans@birth.rate[s]
  }
}

#put in time appropriate survival rate otherwise it uses original surv.matrix set up for trans object and keeps constant over time
#if (!is.na(exper@surv.each.timestep[1,1])) {
#  for (s in 1:n.subpops) {
#    last.index.tmp1 <- (s*exper@trans@n.age.class)
#    first.index.tmp1 <- last.index.tmp1-(exper@trans@n.age.class)+1
#    tmp.age.surv.matrix <- ExtractAgeSpecificSurvivalMatrix(tmp.trans=tmp.trans, maternal.obj=exper@maternal.obj,
#                                                            surv.at.timestep.t=exper@surv.each.timestep[first.index.tmp1:last.index.tmp1,t])
#    last.index.tmp2 <- (s*exper@trans@n.epi.class*exper@trans@n.age.class)
#    first.index.tmp2 <- last.index.tmp2-(exper@trans@n.epi.class*exper@trans@n.age.class)+1
#    tmp.trans@age.surv.matrix[first.index.tmp2:last.index.tmp2,] <- tmp.age.surv.matrix
#  }
#}

##put in correct vaccination coverage
for (s in 1:n.subpops) {
  last.index.tmp <- (s*exper@trans@n.age.class)
  first.index.tmp <- last.index.tmp-(exper@trans@n.age.class)+1

  routine.vacc.prob <- routine$age.time.specific.routine[first.index.tmp:last.index.tmp,]
  sia.vacc.prob <- SIA$age.time.specific.SIA[first.index.tmp:last.index.tmp,]

  if (!is.na(index.sia.vacc[t])){ #if SIA
    tmp.trans@vac.per@pvacc.in.age.class[first.index.tmp:last.index.tmp] <-
      routine.vacc.prob[,index.routine.vacc[t]] +
      sia.vacc.prob[,index.sia.vacc[t]] -
      (routine.vacc.prob[,index.routine.vacc[t]]*sia.vacc.prob[,index.sia.vacc[t]])
    #stow output
    MR1.fail.each.timestep[t] <- routine$prop.fail.MR1[s,index.routine.vacc[t]]
    MR2.fail.each.timestep[t] <- routine$prop.fail.MR2[s,index.routine.vacc[t]]
    SIA.fail.each.timestep[t] <- SIA$prop.fail.SIA[s,index.sia.vacc[t]]

  } else { #if no SIA
    tmp.trans@vac.per@pvacc.in.age.class[first.index.tmp:last.index.tmp] <-
      routine.vacc.prob[,index.routine.vacc[t]]
    #stow output
    MR1.fail.each.timestep[t] <- routine$prop.fail.MR1[s,index.routine.vacc[t]]
    MR2.fail.each.timestep[t] <- routine$prop.fail.MR2[s,index.routine.vacc[t]]
  }
}

ts1 <- tmp.trans@age.surv.matrix[41:45, 41:45]

#run experiment to get the next time step
#state2 <- next.ID.state(state, tran=tmp.trans)
#state1_space <- state[1:1500,1]
#state2_space <- state2[1:1500,1]

state1_space <- state
tran=tmp.trans

##next.ID.state####
#get the transition matrix
tran.matrix <- tran@age.surv.matrix

#get a spatial index on same scale as epi indexes
index.loc <-  tran@subpop.class.label[tran@i.inds]
one.loc <- index.loc==1

n=1
this.loc <- index.loc==n

#The vaccination logic
age.spec.vacc.prob <- tran@vac.per@pvacc.in.age.class[this.loc]

##Updating transition matrix to include vaccination
#M->V transition
tran.matrix[tran@v.inds[this.loc],tran@m.inds[one.loc]] <-  tran@age.surv.matrix[tran@v.inds[this.loc],tran@m.inds[one.loc]] *
  age.spec.vacc.prob
#M->M transition
tran.matrix[tran@m.inds[this.loc],tran@m.inds[one.loc]] <-  tran@age.surv.matrix[tran@m.inds[this.loc],tran@m.inds[one.loc]] *
  (1-age.spec.vacc.prob)
#M->S transition
tran.matrix[tran@s.inds[this.loc],tran@m.inds[one.loc]] <-  tran@age.surv.matrix[tran@s.inds[this.loc],tran@m.inds[one.loc]] *
  (1-age.spec.vacc.prob)
#S->V transition
tran.matrix[tran@v.inds[this.loc],tran@s.inds[one.loc]] <-  tran@age.surv.matrix[tran@v.inds[this.loc],tran@s.inds[one.loc]] *
  age.spec.vacc.prob
#S->S transition
tran.matrix[tran@s.inds[this.loc],tran@s.inds[one.loc]] <-  tran@age.surv.matrix[tran@s.inds[this.loc],tran@s.inds[one.loc]] *
  (1-age.spec.vacc.prob)
#S->I transition
tran.matrix[tran@i.inds[this.loc],tran@s.inds[one.loc]] <-  tran@age.surv.matrix[tran@i.inds[this.loc],tran@s.inds[one.loc]] *
  (1-age.spec.vacc.prob)

ts2 <- tran.matrix[41:45, 41:45]

## Calculate the phi matrix
#define denominator of phi matrix - preventing NAs by setting min to 1
if (tran@frequency.dep) {denom<-max(sum(state[tran@subpop.class.label==n]),1)} else {denom<-1}

#phi <- as.double(tran@waifw)*((state[tran@i.inds[this.loc],]^tran@exponent)/denom)
phi <- (tran@waifw%*%(state[tran@i.inds[this.loc],]^tran@exponent))/denom
phi <- 1 - exp(-phi)
phi <- matrix(phi, nrow=state@n.age.class ,
              ncol=state@n.age.class)
phi <- t(phi)

#print(c("phi",range(phi)))

##Now that phi is calculated, update the tran matrix
#people who stay susceptible
tran.matrix[tran@s.inds[this.loc], tran@s.inds[one.loc]] <-
  tran.matrix[tran@s.inds[this.loc], tran@s.inds[one.loc]] * (1-phi)

#people who become infected
tran.matrix[tran@i.inds[this.loc], tran@s.inds[one.loc]] <-
  tran.matrix[tran@i.inds[this.loc], tran@s.inds[one.loc]] * (phi)

#no one stays infected
tran.matrix[tran@i.inds[this.loc], tran@i.inds[one.loc]] <- 0

ts3 <- tran.matrix[41:45, 41:45]

n=1
here <- which(tran@subpop.class.label==n,arr.ind=TRUE)
state2_space <- tran.matrix[here,]%*%state[here,]



##non-space####

tmp.trans <- exper@trans
state <- exper@state.t0[1:1500]
t=2

#scale the waifw by seasonality
if (!is.array(mults)){
  tmp.trans@waifw <- exper@trans@waifw*mults[t]
} else {
  tmp.trans@waifw <- exper@trans@waifw*mults[,,t]
}


tmp.trans@age.surv.matrix <- ExtractAgeSpecificSurvivalMatrix(tmp.trans=tmp.trans, maternal.obj=exper@maternal.obj,
                                                                surv.at.timestep.t=exper@surv.each.timestep[1:300,t])

#put in the correct birth rate for that time-step, if it varies
if (length(exper@births.per.1000.each.timestep)>1) {
  # Read 'birth rate'
  tmp.trans@birth.rate = (exper@births.per.1000.each.timestep[1,t]*sum(state)/1000)
}

##put in correct vaccination coverage
routine.vacc.prob <- routine$age.time.specific.routine[1:300,]
sia.vacc.prob <- SIA$age.time.specific.SIA[1:300]
#if no SIA
tmp.trans@vac.per@pvacc.in.age.class <- routine.vacc.prob[,1]

tran <- tmp.trans
#state <- next.ID.state(state, tmp.trans)

tn1 <- tran@age.surv.matrix[41:45, 41:45]
tran@age.surv.matrix <- .Call("update_age_surv_MSIRV",
                              tran@age.surv.matrix,
                              sz = nrow(tran@age.surv.matrix),
                              tran@vac.per@pvacc.in.age.class,
                              tran@v.inds[1:300],
                              tran@m.inds[1:300],
                              tran@s.inds[1:300],
                              tran@i.inds[1:300])

tn2 <- tran@age.surv.matrix[41:45, 41:45]

#SIRID####
if (tran@frequency.dep) {denom<-sum(state)} else {denom<-1}
tran.matrix <- .Call("calc_phi_and_update_tran",
                     tran@waifw,
                     state,
                     tran@s.inds[1:300],
                     tran@i.inds[1:300],
                     tran@exponent,
                     denom,
                     tran@age.surv.matrix)

tn3 <- tran@age.surv.matrix[41:45, 41:45]

state2 <- tran.matrix%*%state

state2_nospace <- state2
state1_nospace <- state
