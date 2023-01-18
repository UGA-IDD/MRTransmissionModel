## add these functions


## Revised & New functions
## 29 November 2022

## New Functions - These functions are all new function####

#Class definition for space.nMx object - aka Age Specific Death Rates (ASDR) - by space
# setClass("space.nMx",
#          representation(rate.years="numeric", #year associated with each time-specific rate
#                         rates="data.frame", #age and space specific death rates per 1 (rows) by year (columns)
#                         mid.age="numeric", #mid-age associated with each age-specific rates
#                         n.subpops="numeric" #number of subpopulations
#          ))

#Wrapper function to add space to create.surv.prob.over.age.time()
#
#Paramaters -
#   age.classes - vector - a set of age classes for which tran is being built
#   generation time - numeric - generation time
#   space.nMx - space.nMx object -  age and space specific death rates over time
#   years.interpolate - vector - years to interpolate age specific death rates
#   check - boolean - whether or not to plot the death rates over age
#
#Returns-
#   returns the age and spatial profile of survivorship (in units of the generation time) (rows) by time (cols)
# space.wrapper.surv.prob.over.age.time <- function(age.classes, generation.time, space.nMx=NULL, years.interpolate = seq(1980,(1980+121),1), check=F){
#
#   n.subpops <- space.nMx@n.subpops
#
#   for (s in 1:n.subpops){
#     n.age.classes <- length(space.nMx@mid.age)
#     nMx.tmp <- new("nMx",
#                    rate.years = space.nMx@rate.years,
#                    rates = space.nMx@rates[((s*n.age.classes)-n.age.classes+1):(s*n.age.classes),],
#                    mid.age = space.nMx@mid.age)
#     surv.tmp <- create.surv.prob.over.age.time(age.classes, generation.time, nMx = nMx.tmp, nMx.years = years.interpolate, check=check)
#     if (s==1) surv <- surv.tmp
#     if (s!=1) surv <- rbind(surv, surv.tmp)
#   }
#
#   return(surv)
#
# }


#Wrapper function to add space to create.surv.prob.over.age()
#
#Paramaters -
#   age.classes - vector - a set of age classes for which tran is being built
#   generation time - numeric - generation time
#   space.nMx - space.nMx object -  age and space specific death rates over time
#   year - year to pull age specific death rates that will be used to run out the transients (shoudld coincide with DFE year)
#
#Returns-
#   returns the age profile (cols) of survivorship (in units of the generation time) by spatial unit (rows)
# space.wrapper.surv.prob.over.age <- function(age.classes, generation.time, space.nMx=NULL, year = 1980, check=F){
#
#   n.subpops <- space.nMx@n.subpops
#
#   for (s in 1:n.subpops){
#     n.age.classes <- length(space.nMx@mid.age)
#     nMx.tmp <- new("nMx",
#                    rate.years = space.nMx@rate.years,
#                    rates = space.nMx@rates[((s*n.age.classes)-n.age.classes+1):(s*n.age.classes),],
#                    mid.age = space.nMx@mid.age)
#     surv.tmp <- create.surv.prob.over.age(age.classes, generation.time, nMx = nMx.tmp, year, check=T)
#     if (s==1) surv <- surv.tmp
#     if (s!=1) surv <- rbind(surv, surv.tmp)
#   }
#
#   return(surv)
#
# }


#Function to scale the WAIFW matrix to a particular R0 given a particular subpopulation, but for a spatially structured population
#
#Parameters -
#    R0 - the reproductive rate to scale to
#    state - a population at disease free equilibrium
#    waifw - the matrix to scale
#
#Returns -
#    a scaled versions of the age-specific waifw
# scale.WAIFW.space <- function(R0, state, waifw, frequency.dep=F, suscept.state = 1) {
#
#   #move everyone into susceptible category for DFE
#   DFE.state <- matrix(0,length(state[,1]),1)
#   for (k in  1:state@n.epi.class)
#     DFE.state[state@epi.class==suscept.state,1] <-
#       DFE.state[state@epi.class==suscept.state,1]+state[state@epi.class==k,1]
#
#   #state.here <- rowSums(matrix(DFE.state[state@epi.class==suscept.state,1],
#   #                     state@n.age.class,length(state[state@epi.class==suscept.state,1])/state@n.age.class))
#
#   #just use first site (not whole pop)
#   state.here <- DFE.state[state@epi.class==suscept.state,1][1:state@n.age.class]
#   if (frequency.dep) denom <- sum(state[,1]) else denom <- sum(state[1:(state@n.age.class*state@n.epi.class),1])
#   #plot(state.here,type="l")
#
#   #more correct
#   next.gen <- state.here*(1-exp(-waifw/denom))
#
#   #get the first eigen value
#   cur.R0 <- Re(eigen(next.gen)$value[1])
#
#   #More correct transform
#   R.ratio <- R0/cur.R0; #print(R0); #print(cur.R0); #print(R.ratio)
#   waifw <- -log(1-R.ratio*(1-exp(-waifw/denom)))*denom
#
#   return(waifw)
# }


#Function creates an ID state matrix with the given number of
#epidemiological categories, age classes, and subpopulations
#
#Parameters -
#   n.age.class - the number of age classes
#   n.epi.class - the number of states in the model
#   epi.class.label - the names of the states
#   value - starting values, defaults to all 0
#
#Returns -
#  an ID.state.matrix object
# space.create.ID.state.matrix <- function(n.age.class,
#                                          n.epi.class = 3,
#                                          epi.class.label = c("S","I","R"),
#                                          n.subpops = 1,
#                                          value = rep(0,n.age.class*n.epi.class*n.subpops)) {
#
#   rc <- new("ID.state.matrix",
#             nrow = n.age.class*n.epi.class*n.subpops,
#             ncol = 1,
#             n.epi.class = n.epi.class,
#             epi.class = rep(rep(1:n.epi.class, n.age.class),n.subpops),
#             epi.class.label = epi.class.label,
#             n.age.class = n.age.class)
#
#   rc[,1] <- value
#   return(rc)
# }

#Function to create a state matrix at the disease free equilibrium (DFE)
#from a transition object
#
#Parameters -
#     uncode - numeric- UN country code
#     tot.subpop - vector - sub-populations by which you can rescale if you want
#     tran - the transition object
#     year - numeric - options are 1950 to 2100
#
#Return -
#     a starting DFE vector of S and R individuals
# space.create.country.x.DFE.ID.state.matrix <- function(uncode, tot.subpop=NULL, tran,
#                                                        epi.class.label = c("M","S","I","R","V"),
#                                                        year=1990){
#
#   #make template
#   rc <- space.create.ID.state.matrix(tran@n.age.class,
#                                      tran@n.epi.class,
#                                      epi.class.label = epi.class.label,
#                                      tran@n.subpops)
#
#   #fill in the template
#   demog <- space.getDemography(uncode=uncode)
#   demog.n.ages <- nrow(demog$pop.age.byageclasses.1950.2100)/tran@n.subpops
#
#   #loop through the number of sub-populations
#   for (s in 1:tran@n.subpops){
#
#     pop.struct <- demog$pop.age.byageclasses.1950.2100[(demog.n.ages*s-demog.n.ages+1):(demog.n.ages*s),(year-1950+1)]
#     age <- as.numeric(names(pop.struct))
#
#     if (!is.null(tot.subpop)) {
#       tot.pop <- tot.subpop[s]
#     } else {
#       tot.pop <- demog$pop.total.1950.2100[s,(year-1950+1)]
#     }
#
#     #Turn this into prop of desired age classes
#     #note that NUMBER is imposed externally (since we might not want full countries... because interested in stochasticity)
#     sp <- smooth.spline(age, pop.struct)  # gives "slope" between age (range 0-100) and population size per age (therefore 101 items)
#     pred <- predict(sp, tran@age.class/12)$y # predicts population size per age class (total 280 age classes)
#     #plot(age,pop.struct, ylim=c(0,3000000))
#     #points(tran@age.class/12,pred,type="l",col=2)
#
#     #Adjust for varying bin width
#     pred <- pred*diff(c(0,tran@age.class))
#     #added a 0 at the beginning of age classes
#     #took the difference between each class n from n-1.   therefore difference is 1 (or 1 month) for item 1-241 and then 12 for items 242-281
#     #then weighing each predicted population size per age class by the size of the age class
#
#     #Fill in
#     rc[tran@s.inds[(tran@n.age.class*s-tran@n.age.class+1):(tran@n.age.class*s)],1] <- tot.pop*pred/sum(pred)
#     #if pop=NULL then assuming that the population w/in each state with missing ages are distributed within the population based on the pop structure that we do have ages for
#     #filling in the susceptibles per age class
#
#   }
#
#   return(rc)
#
# }


#Class to represent transitions in the context of vaccination and maternal antibodies with multiple locations
# setClass("ID.transition.MSIRV.space",
#          representation(n.subpops = "numeric", #total number of locations
#                         subpop.class.label = "numeric", #location index
#                         coupling = "matrix"), #matrix of connection between each sub-population
#          contains = "ID.transition.MSIRV.SIA"
# )



# # Method to create an ID.transition.MSIRV.space object
# space.create.ID.transition.MSIRV <- function(n.age.class, #number of age classes (not location specific)
#                                              aging.rate,  #aging rates (includes all locations concatenated)
#                                              survival.rate, #vector with locations survival rates concatenated
#                                              waifw, # assume same across locations
#                                              routine.vac = 0,
#                                              routine.vac.age.index = 12,
#                                              maternal.obj, #maternal antibody object
#                                              time.step = 1/2, #time step in months == generation time
#                                              age.class = numeric(),
#                                              birth.rate.subpop = rep(0, n.subpops),
#                                              sia.vac = 0,
#                                              sia.vsucc =  new("vsucc.constant", success.rate=1),
#                                              introduction.rate = 0,
#                                              exponent = exponent,
#                                              frequency.dep = F,
#                                              is.stochastic = F,
#                                              get.births = function(pop,tran){return(c(tran@birth.rate,rep(0,length(pop)-1)))},
#                                              n.subpops = 1, #number of sub-populations
#                                              subpop.class.label,
#                                              coupling = 1 #matrix describing movement between sub-populations
# ) {
#
#
#   #First do all of the standard logic
#   rc <- new("ID.transition.MSIRV.space",
#             n.age.class = n.age.class,
#             n.epi.class = 5,
#             epi.class.label = c("M","S","I","R","V"),
#             epi.class = rep(rep(1:5, n.age.class),n.subpops),
#             age.class = age.class,
#             aging.rate = aging.rate,
#             survival.rate = survival.rate,
#             birth.rate = birth.rate.subpop,
#             waifw = waifw,
#             sia.vsucc = sia.vsucc,
#             introduction.rate = introduction.rate,
#             exponent = exponent,
#             frequency.dep=frequency.dep,
#             is.stochastic=is.stochastic,
#             get.births = get.births,
#             n.subpops=n.subpops,
#             subpop.class.label=rep(1:n.subpops,each=n.age.class*5),
#             coupling=coupling)
#
#   rc@m.inds = which(rc@epi.class==1)
#   rc@s.inds = which(rc@epi.class==2)
#   rc@i.inds = which(rc@epi.class==3)
#   rc@r.inds = which(rc@epi.class==4)
#   rc@v.inds = which(rc@epi.class==5)
#
#   rc@age.surv.matrix <- matrix(0,nrow = n.age.class*rc@n.epi.class*n.subpops,
#                                ncol = n.age.class*rc@n.epi.class)
#
#   #for zeroing out all the appropriate transitions
#   template.mtrx <- matrix(c(1,1,0,0,1,0,1,1,0,1,0,0,1,1,0,0,0,0,1,0,0,0,0,0,1), nrow=5, ncol=5)
#
#   #useful low age
#   low.age <-rep(c(0, age.class[2:length(age.class)-1]), n.subpops)
#
#   #for each subpop and age class
#   for (j in 1:n.subpops) {
#     for (i in 1:(n.age.class)) {
#       #first fill in the diagonal matrix
#       tmp.template <- template.mtrx
#       tmp.template[2,1] <- 0
#       #move into right spatial location
#       adj <- (j-1)*n.age.class*rc@n.epi.class
#       #index ages
#       inds <- adj+(i-1)*rc@n.epi.class+(1:rc@n.epi.class)
#       rc@age.surv.matrix[inds,inds-adj] <- (1-aging.rate[i])*survival.rate[j,i]*tmp.template
#
#       #now fill in the off diagonal matrix
#       if ((i%%n.age.class)!=0) {
#         #create a modified template matrix to move
#         #people from the M class to the S class
#         tmp.template <- template.mtrx
#
#         #calculate the probability of losing maternal protection during
#         #an age class
#         p.mat.loss <- (pmaternal(low.age[i], maternal.obj)-
#                          pmaternal(age.class[i], maternal.obj))/pmaternal(low.age[i],maternal.obj)
#         p.mat.loss[pmaternal(low.age[i],maternal.obj)==0] <- 1
#
#         tmp.template[1,1] <- 1-p.mat.loss
#         tmp.template[2,1] <- 1-tmp.template[1,1]
#
#         inds2 <- adj + (i)*rc@n.epi.class+(1:rc@n.epi.class)
#         #print(inds2)
#         rc@age.surv.matrix[inds2,inds-adj] <- aging.rate[i]*survival.rate[j,i]*tmp.template
#       }
#     }
#   }
#
#
#   #Next do all of the vaccination logic, concatenating two sites together, with desired coverage
#   #get the %vaccinated in each age group per time step time step must be in years!
#   vac.per <- c()
#   for (l in 1:n.subpops) {
#     vac.per.loc <- new("vacc.per.time.step.by.age",
#                        pvacc.in.age.class = numeric(length=n.age.class))
#     vac.per.loc@pvacc.in.age.class[routine.vac.age.index[l]] <- routine.vac[l]
#     vac.per <- c(vac.per, vac.per.loc@pvacc.in.age.class)
#   }
#
#   #rc@vac.per <- vac.per #non-patch uses this
#   rc@vac.per@pvacc.in.age.class <- vac.per
#   rc@sia.vac <- sia.vac
#
#   return(rc)
# }

#Function to get starting state vector and transition object - built originally from Get.CountryX.Starting.Pop.MSIRV()
#
#Parameters -
#     uncode - UN country code
#     generation.time - the desired generation time in months
#     age.classes - vector - the upper limit of the age classes we are interested in months
#     maternal.decay.rt -rate of maternal decay of immunity to rubella
#     exponent - numeric - exponent for the infected
#     frequency.dep - boolean - TRUE
#     is.stochastic - boolean - FALSE
#     tot.subpop - numeric vector - population you want to scale you whole experiment by at the beginning by spatial unit, NULL or numeric
#     yr.births.per.1000.bysubpop - numeric vector - crude birth rate per year per 1000 by spatial unit
#     intro.rate - numeric -
#     flat.WIAFW - boolean - if F then default waifw is polymod GB
#     space.asdr.object - space.nMx object with rate and mid-age and year and space (age specific death rates)
#     targeted.intro - boolean - FALSE means space introduction out over all age classes, TRUE means to concentrate introductions
#     year - numeric - year to pull country DFE demography (age structure and population size)
#     get.births - vector of births with length of 5*n.age.classes
#     routine.vac - routine vaccination coverage for each subpop
#     routine.vac.age.index - age index based on age classes for each subpop
#     n.subpops - the number of subpopulations
#     coupling - matrix - coupling of subpopulations
#
#Returns -
#     starting state vector and transition object
# Space.Get.CountryX.Starting.Pop.MSIRV <- function(uncode ,
#                                                   generation.time = 0.5,  #generation time in months
#                                                   age.classes = c(1:240, seq(241,720,12)),
#                                                   maternal.decay.rt=0.95, #based on Metcalf, 2012 paper
#                                                   exponent = 0.97,
#                                                   frequency.dep=TRUE,
#                                                   is.stochastic=FALSE,
#                                                   tot.subpop=NULL,
#                                                   yr.births.per.1000.bysubpop,
#                                                   intro.rate=intro.rate,
#                                                   flat.WAIFW=F,
#                                                   space.asdr.object=NULL,
#                                                   targeted.intro=FALSE,
#                                                   year=1980,
#                                                   get.births,
#                                                   routine.vac=0,
#                                                   routine.vac.age.index=12,
#                                                   n.subpops = 1,
#                                                   coupling){
#
#   ## Calculate the aging rate using the age classes and the generation time
#   age.lows <- c(0,age.classes[2:length(age.classes)-1]) # makes it 0 to 697 rather than 1 to 709, so lowest age in wach age range
#   ac.sz <- age.classes - age.lows # the size in months of each age range (1 month till age 20, then 12 months to age 59)
#   aging.rate <- generation.time/ac.sz # converting generation time units into year/month time units based on size of age class
#   aging.rate[length(aging.rate)] <- 0 # forcing the last aging range to 0
#
#   ## Returns the age profile (cols) of survivorship (in units of the generation time) by spatial unit (rows)
#   survs <- space.wrapper.surv.prob.over.age(age.classes, generation.time, space.nMx=space.asdr.object, year)
#
#   ## Create WAIFW matrix age.classes X age.classes dimensions, default is POLYMOD based on Great Britain
#   waifw <- get.polymod.WAIFW(age.classes/12)
#   if (flat.WAIFW) waifw <- get.flat.WAIFW(age.classes/12)
#
#   ## Set up maternal immunity
#   maternal.obj = new("maternal.exp.decay", decay.rt=maternal.decay.rt)
#
#   ## Create the transition object
#   tran <- space.create.ID.transition.MSIRV(n.age.class = length(age.classes),
#                                            aging.rate = aging.rate,
#                                            survival.rate = survs,
#                                            waifw = waifw,
#                                            routine.vac = routine.vac,
#                                            routine.vac.age.index = routine.vac.age.index,
#                                            maternal.obj = maternal.obj,
#                                            time.step = generation.time,
#                                            age.class = age.classes,
#                                            #birth.rate,
#                                            exponent = exponent,
#                                            frequency.dep=frequency.dep,
#                                            is.stochastic=is.stochastic,
#                                            get.births=get.births,
#                                            #sia.vac = 0,
#                                            #sia.vsucc =  new("vsucc.constant", success.rate=1),
#                                            #introduction.rate = 0,
#                                            n.subpops = n.subpops, #number of sub-populations
#                                            coupling = matrix(1, nrow=n.subpops, ncol=n.subpops))
#
#
#   ## Putting in starting state where everyone susceptible
#   state <- space.create.country.x.DFE.ID.state.matrix(uncode=uncode, tot.subpop=tot.subpop, tran=tran,
#                                                       epi.class.label = c("M","S","I","R","V"), year=year)
#
#   ## Update starting births now because need sum(state) = pop for each subpop
#   ## the expected number of total births for each time step
#   for (s in 1:n.subpops){
#     subpop.tmp <- sum(state[tran@s.inds[(tran@n.age.class*s-tran@n.age.class+1):(tran@n.age.class*s)],1])
#     tran@birth.rate[s] <- (yr.births.per.1000.bysubpop[s]/1000*subpop.tmp)*generation.time/12
#   }
#
#   ## Input transition introduction of infection rate
#   if (targeted.intro) {
#     tran@introduction.rate <- rep(0, tran@n.age.class*tran@n.subpops)
#     tran@introduction.rate[(60:71)*(rep(1:tran@n.subpops, each=length(60:71)))] <- intro.rate
#   } else {
#     tran@introduction.rate <- rep(intro.rate, tran@n.age.class*tran@n.subpops)
#   }
#
#   ## Assumes everyone has maternal protection
#   state[tran@m.inds,1] <- state[tran@s.inds,1] * pmaternal(age.lows, maternal.obj)
#   state[tran@s.inds,1] <- state[tran@s.inds,1] * (1-pmaternal(age.lows, maternal.obj))
#
#   return(list(state = state, tran = tran))
# }


#' #Wrapper function to add space to get.routine.time.age.specific()
#'
#' @param time.step numeric, time step in months
#' @param age.classes vector, age classes in months, same as tran object
#' @param space.time.specific.MR1cov matrix, space (rows) and year (columns) specific routine MR1 coverage as a proportion
#' @param age.min.MR1 vector, time specific minimum age eligible for MR1
#' @param age.max.MR1 vector, time specific maximum age eligible for MR1
#' @param space.time.specific.MR2cov matrix, space (rows) and year (columns) specific routine MR2 coverage as a proportion
#' @param age.min.MR2 vector, time specific minimum age eligible for MR2
#' @param age.max.MR2 vector, time specific maximum age eligible for MR2
#' @param obj.vcdf.MR1 vaccine.cdf.byage object, age distribution of routine MR1
#' @param obj.vcdf.MR2 vaccine.cdf.byage object, age distribution of routine MR2
#' @param obj.prob.vsucc prob.vsucc.byage object, VE by age relevant to both routine vacccine doses
#' @param MR1MR2correlation boolean, if MR1 and MR2 are dependent (TRUE) or independent (FALSE)
#'
#' @return list, features specific to routine vaccination, including age.time.specific.routine matrix which is age and space (rows) by year (cols) of proportion routine vaccination
#' @export
#'
#' @examples
# space.wrapper.get.routine.time.age.specific <- function(time.step=0.5, age.classes = c(1:240, seq(252,1212,12)),
#                                                         space.time.specific.MR1cov = rbind(rep(0.7, 50),rep(0.6, 50)), age.min.MR1=rep(12, 50), age.max.MR1=rep(23,50),
#                                                         space.time.specific.MR2cov = rbind(rep(0.4, 50),rep(0.3, 50)), age.min.MR2=rep(24, 50), age.max.MR2=rep(35, 50),
#                                                         obj.vcdf.MR1=get.vcdf.uniform(12, 23), obj.vcdf.MR2=get.vcdf.uniform(24, 35),
#                                                         obj.prob.vsucc = pvacsuccess(1:(14*12), new("vsucc.constant", success.rate=0.8)),
#                                                         MR1MR2correlation=F){
#   n.subpops <- nrow(space.time.specific.MR1cov)
#   for (s in 1:n.subpops){
#     tmp <- get.routine.time.age.specific(time.step, age.classes,
#                                          time.specific.MR1cov=space.time.specific.MR1cov[s,], age.min.MR1, age.max.MR1,
#                                          time.specific.MR2cov=space.time.specific.MR2cov[s,], age.min.MR2, age.max.MR2,
#                                          obj.vcdf.MR1, obj.vcdf.MR2,
#                                          obj.prob.vsucc, MR1MR2correlation)
#
#     #we need to transpose these matrices so that age (rows) and time (cols)
#     tmp$age.time.specific.routine <- t(tmp$age.time.specific.routine)
#     tmp$prop.fail.MR1.byage <- t(tmp$prop.fail.MR1.byage)
#     tmp$prop.fail.MR2.byage <- t(tmp$prop.fail.MR2.byage)
#
#     if (s==1) {
#       out <- tmp
#     }
#
#     # row bind so that each new sub-population adds length(age.classes) number of rows to the bottom of the matrix
#     if (s!=1){
#       out$age.time.specific.routine <- rbind(out$age.time.specific.routine, tmp$age.time.specific.routine)
#       out$prop.fail.MR1.byage <- rbind(out$prop.fail.MR1.byage, tmp$prop.fail.MR1.byage)
#       out$prop.fail.MR2.byage <- rbind(out$prop.fail.MR2.byage, tmp$prop.fail.MR2.byage)
#       out$one.minus.ve1 <- rbind(out$one.minus.ve1, tmp$one.minus.ve1)
#       out$one.minus.ve2 <- rbind(out$one.minus.ve2, tmp$one.minus.ve2)
#       out$prop.fail.MR1 <- rbind(out$prop.fail.MR1, tmp$prop.fail.MR1)
#       out$prop.fail.MR2 <- rbind(out$prop.fail.MR2, tmp$prop.fail.MR2)
#     }
#   }
#   return(out)
# }

#' Wrapper function to add space to get.sia.time.age.specific()
#'
#' @param age.classes vector, age classes in months, same as tran object
#' @param space.time.specific.SIAcov matrix, space (rows) and year (columns) specific campaign MR coverage as a proportion
#' @param age.min.sia vector, time specific minimum age eligible for MR campaign dose
#' @param age.max.sia vector, time specific maximum age eligible for MR campaign dose
#' @param obj.prob.vsucc prob.vsucc.byage object, VE by age relevant to campaign vacccine doses
#'
#' @return list, features specific to campaign vaccination, including age.time.specific.SIA matrix which is age and space (rows) by year (cols) of proportion campaign vaccination
#' @export
#'
#' @examples
# space.wrapper.get.sia.time.age.specific <- function(age.classes=c(1:240, seq(252,1212,12)),
#                                                     space.time.specific.SIAcov=rbind(c(rep(0,29),0.5,rep(0,19),0.7),c(rep(0,29),0.7,rep(0,19),0.8)),
#                                                     age.min.sia=rep(12,50),
#                                                     age.max.sia=rep(60,50),
#                                                     obj.prob.vsucc=pvacsuccess(1:240, new("vsucc.constant", success.rate=0.85))){
#   n.subpops <- nrow(space.time.specific.SIAcov)
#   for (s in 1:n.subpops){
#     tmp <- get.sia.time.age.specific(age.classes=age.classes,
#                                      time.specific.SIAcov=space.time.specific.SIAcov[s,],
#                                      age.min.sia=age.min.sia,
#                                      age.max.sia=age.max.sia,
#                                      obj.prob.vsucc=obj.prob.vsucc)
#
#     #we need to transpose these matrices so that age (rows) and time (cols)
#     tmp$age.time.specific.SIA <- t(tmp$age.time.specific.SIA)
#
#     if (s==1) {
#       out <- tmp
#     }
#
#     # row bind so that each new sub-population adds length(age.classes) number of rows to the bottom of the matrix
#     if (s!=1){
#       out$age.time.specific.SIA <- rbind(out$age.time.specific.SIA, tmp$age.time.specific.SIA)
#       out$prop.fail.SIA <- rbind(out$prop.fail.SIA, tmp$prop.fail.SIA)
#     }
#   }
#   return(out)
# }

#Class for a spatial experiment updating demography and with changing vaccination coverage over time
# setClass("experiment.updatedemog.vaccinationchange.spatial",
#          representation(trans = "ID.transition.MSIRV.space"),
#          contains=c("experiment.updatedemog.vaccinationchange"))

# # Run method for experiment.updatedemog.vaccinationchange.spatial objects
# setMethod("run",
#           "experiment.updatedemog.vaccinationchange.spatial",
#           function(exper, rescale.WAIFW=T, ...) {
#
#             #print(rescale.WAIFW)
#             state <- exper@state.t0
#
#             #number of sub-populations
#             n.subpops <- exper@trans@n.subpops
#
#             #rescale the WAIFW if a specific R0 is specified
#             if (rescale.WAIFW & length(exper@R0)>0) {
#               #print("RESCALING!!!")
#               exper@trans@waifw <- scale.WAIFW.space(exper@R0,
#                                                      state,exper@trans@waifw,
#                                                      frequency.dep=exper@trans@frequency.dep,
#                                                      suscept.state=exper@trans@s.inds[1])
#             }
#
#             #get the number of time steps in the experiment
#             T <- round((exper@t.max-exper@t.min)/exper@step.size)+1
#
#             #get seasonal mult
#             if (!is.null(exper@season.obj)) {
#               mults <- get.seasonal.mult(exper@t0.doy/365+(1:T-1)*exper@step.size,
#                                          exper@season.obj)
#             } else {
#               mults <- rep(1,T)
#             }
#
#             #make a temporary transmission object
#             tmp.trans <- exper@trans
#
#             #hold the states as we walk through
#             rc <- matrix(ncol = T, nrow = nrow(state))
#             rc[,1] <- state
#
#             #need output vector for births over time
#             births.each.timestep <- growth.rate.each.timestep <- N0 <- matrix(NA, n.subpops, T)
#             births.each.timestep[,1] <- tmp.trans@birth.rate
#
#             #generate the age and space (rows) and year (columns) specific vaccination matrix
#             routine <- space.wrapper.get.routine.time.age.specific(time.step= exper@step.size*12,
#                                                                    age.classes=exper@trans@age.class,
#                                                                    space.time.specific.MR1cov=exper@time.specific.MR1cov,
#                                                                    age.min.MR1=exper@time.specific.min.age.MR1,
#                                                                    age.max.MR1=exper@time.specific.max.age.MR1,
#                                                                    space.time.specific.MR2cov=exper@time.specific.MR2cov,
#                                                                    age.min.MR2=exper@time.specific.min.age.MR2,
#                                                                    age.max.MR2=exper@time.specific.max.age.MR2,
#                                                                    obj.vcdf.MR1=exper@obj.vcdf.MR1,
#                                                                    obj.vcdf.MR2=exper@obj.vcdf.MR2,
#                                                                    obj.prob.vsucc=exper@obj.prob.vsucc,
#                                                                    MR1MR2correlation=F)
#
#             #getting year of routine introductions
#             routine.intro <- rep(0, T)
#             if (any(colSums(exper@time.specific.MR1cov)!=0)) routine.intro[min(which(colSums(exper@time.specific.MR1cov)>0))*(1/exper@step.size)+1] <- 1
#             if (any(colSums(exper@time.specific.MR2cov)!=0)) routine.intro[min(which(colSums(exper@time.specific.MR2cov)>0))*(1/exper@step.size)+1] <- 1
#
#             #specifying that each year (column) of routine should be repeated 24 times
#             index.routine.vacc <- c(1,rep(1:ncol(routine$age.time.specific.routine), each=(T-1)/exper@t.max))
#
#             #generate the age and space (rows) and year (columns) specific vaccination matrix - same size as `routine`
#             SIA <- space.wrapper.get.sia.time.age.specific(age.classes=exper@trans@age.class,
#                                                            space.time.specific.SIAcov=exper@time.specific.SIAcov,
#                                                            age.min.sia=exper@time.specific.min.age.SIA,
#                                                            age.max.sia=exper@time.specific.max.age.SIA,
#                                                            obj.prob.vsucc=exper@obj.prob.vsucc)
#
#             #getting sia.times as a vector of 0 / 1 to represent when the SIA is to take place
#             index.sia.vacc <- rep(NA,T)
#             year.sia <- which(colSums(exper@time.specific.SIAcov)!=0)
#             index.sia.vacc[(year.sia-1)*(T-1)/exper@t.max + round((exper@sia.timing.in.year*(T-1)/exper@t.max))] <-  year.sia #minus 1 because adding the sia.timing
#             sia.times <- ifelse(!is.na(index.sia.vacc), 1, 0)
#
#             #need output vectors for primary vaccination failure over time
#             MR1.fail.each.timestep <- MR2.fail.each.timestep <- SIA.fail.each.timestep <- matrix(0, n.subpops, T)
#             MR1.fail.each.timestep[,1] <- routine$prop.fail.MR1[,1]
#             MR2.fail.each.timestep[,1] <- routine$prop.fail.MR2[,1]
#             SIA.fail.each.timestep[,1] <- 0
#
#             for (t in 2:T) {
#
#               print(t)
#
#               #scale the waifw by seasonality
#               if (!is.array(mults)){
#                 tmp.trans@waifw <- exper@trans@waifw*mults[t]
#               } else {
#                 tmp.trans@waifw <- exper@trans@waifw*mults[,,t]
#               }
#
#               #if exper@pop.rescale.each.timestep[t]==NaN then tmp.trans@pop.rescale==NaN and will not rescale
#               #otherwise if any number other than NaN it will rescale
#               #or if exper@pop.rescale.each.timestep not completed , then numeric(0) and any index [t] is NA
#               for (s in 1:n.subpops) {
#                   if (!is.na(exper@pop.rescale.each.timestep[s,t])) {
#                   last.index.tmp <- (s*exper@trans@n.epi.class*exper@trans@n.age.class)
#                   first.index.tmp <- last.index.tmp-(exper@trans@n.epi.class*exper@trans@n.age.class)+1
#                   state[first.index.tmp:last.index.tmp] <- exper@pop.rescale.each.timestep[s,t]*
#                     (state[first.index.tmp:last.index.tmp]/sum(state[first.index.tmp:last.index.tmp]))
#                 }
#               }
#
#               #put in the correct birth rate for that time-step, if it varies
#               for (s in 1:n.subpops) {
#                 if (ncol(exper@births.per.1000.each.timestep)>1) {
#                   last.index.tmp <- (s*exper@trans@n.epi.class*exper@trans@n.age.class)
#                   first.index.tmp <- last.index.tmp-(exper@trans@n.epi.class*exper@trans@n.age.class)+1
#                   tmp.trans@birth.rate[s] = (exper@births.per.1000.each.timestep[s,t]*
#                                                sum(state[first.index.tmp:last.index.tmp])/1000)
#                 } else {
#                   tmp.trans@birth.rate[s] <- exper@trans@birth.rate[s]
#                 }
#               }
#               births.each.timestep[,t] <-  tmp.trans@birth.rate
#
#               #put in time appropriate survival rate otherwise it uses original surv.matrix set up for trans object and keeps constant over time
#               if (!is.na(exper@surv.each.timestep[1,1])) {
#                 for (s in 1:n.subpops) {
#                   last.index.tmp1 <- (s*exper@trans@n.age.class)
#                   first.index.tmp1 <- last.index.tmp1-(exper@trans@n.age.class)+1
#                   tmp.age.surv.matrix <- ExtractAgeSpecificSurvivalMatrix(tmp.trans=tmp.trans, maternal.obj=exper@maternal.obj,
#                                                                           surv.at.timestep.t=exper@surv.each.timestep[first.index.tmp1:last.index.tmp1,t])
#
#                   last.index.tmp2 <- (s*exper@trans@n.epi.class*exper@trans@n.age.class)
#                   first.index.tmp2 <- last.index.tmp2-(exper@trans@n.epi.class*exper@trans@n.age.class)+1
#                   tmp.trans@age.surv.matrix[first.index.tmp2:last.index.tmp2,] <- tmp.age.surv.matrix
#                 }
#               }
#
#               ##put in correct vaccination coverage
#               for (s in 1:n.subpops) {
#                 last.index.tmp <- (s*exper@trans@n.age.class)
#                 first.index.tmp <- last.index.tmp-(exper@trans@n.age.class)+1
#
#                 routine.vacc.prob <- routine$age.time.specific.routine[first.index.tmp:last.index.tmp,]
#                 sia.vacc.prob <- SIA$age.time.specific.SIA[first.index.tmp:last.index.tmp,]
#
#                 if (!is.na(index.sia.vacc[t])){ #if SIA
#                   tmp.trans@vac.per@pvacc.in.age.class[first.index.tmp:last.index.tmp] <-
#                     routine.vacc.prob[,index.routine.vacc[t]] +
#                     sia.vacc.prob[,index.sia.vacc[t]] -
#                     (routine.vacc.prob[,index.routine.vacc[t]]*sia.vacc.prob[,index.sia.vacc[t]])
#                   #stow output
#                   MR1.fail.each.timestep[t] <- routine$prop.fail.MR1[s,index.routine.vacc[t]]
#                   MR2.fail.each.timestep[t] <- routine$prop.fail.MR2[s,index.routine.vacc[t]]
#                   SIA.fail.each.timestep[t] <- SIA$prop.fail.SIA[s,index.sia.vacc[t]]
#
#                 } else { #if no SIA
#                   tmp.trans@vac.per@pvacc.in.age.class[first.index.tmp:last.index.tmp] <-
#                     routine.vacc.prob[,index.routine.vacc[t]]
#                   #stow output
#                   MR1.fail.each.timestep[t] <- routine$prop.fail.MR1[s,index.routine.vacc[t]]
#                   MR2.fail.each.timestep[t] <- routine$prop.fail.MR2[s,index.routine.vacc[t]]
#                 }
#               }
#
#               #stow the previous population size
#               for (s in 1:n.subpops) {
#                 last.index.tmp <- (s*exper@trans@n.epi.class*exper@trans@n.age.class)
#                 first.index.tmp <- last.index.tmp-(exper@trans@n.epi.class*exper@trans@n.age.class)+1
#                 N0[s,t] <- sum(state[first.index.tmp:last.index.tmp])
#               }
#
#               #run experiment to get the next time step
#               state <- next.ID.state(state, tran=tmp.trans)
#
#               #stow the previous population size
#               for (s in 1:n.subpops) {
#                 last.index.tmp <- (s*exper@trans@n.epi.class*exper@trans@n.age.class)
#                 first.index.tmp <- last.index.tmp-(exper@trans@n.epi.class*exper@trans@n.age.class)+1
#                 NT <- sum(state[first.index.tmp:last.index.tmp])
#                 growth.rate.each.timestep[s,t] <- log(NT/N0[s,t]) #instantaneous biweekly growth rate
#               }
#
#               #print(t)
#               #print(dim(state))
#               rc [,t] <- state
#             }
#
#
#             result <- new("sim.results.MSIRV.update.demog.vaccine.change.space",
#                       data=rc,
#                       m.inds = exper@trans@m.inds,
#                       s.inds = exper@trans@s.inds,
#                       i.inds = exper@trans@i.inds,
#                       r.inds = exper@trans@r.inds,
#                       v.inds = exper@trans@v.inds,
#                       t = exper@t.min+(1:T-1)* exper@step.size,
#                       age.class = exper@trans@age.class,
#                       births.each.timestep = births.each.timestep,
#                       growth.rate.each.timestep = growth.rate.each.timestep,
#                       MR1.fail.each.timestep = MR1.fail.each.timestep,
#                       MR2.fail.each.timestep = MR2.fail.each.timestep,
#                       SIA.fail.each.timestep = SIA.fail.each.timestep,
#                       routine.intro = routine.intro,
#                       sia.times = sia.times,
#                       n.subpops = exper@trans@n.subpops,
#                       subpop.class.label = exper@trans@subpop.class.label)
#
#             rc <- new("experiment.result",
#                       experiment.def = exper,
#                       result = result)
#
#
#             return(rc)
#           }
# )


# #Class for a spatial experiment updating demography over time
# setClass("experiment.updatedemog.spatial",
#          representation(trans = "ID.transition.MSIRV.space"),
#          contains=c("experiment.updatedemog"))

# Run method for experiment.updatedemog.spatial objects
# setMethod("run",
#           "experiment.updatedemog.spatial",
#           function(exper, rescale.WAIFW=T, ...){
#
#             #print(rescale.WAIFW)
#             state <- exper@state.t0
#
#             #number of sub-populations
#             n.subpops <- exper@trans@n.subpops
#
#             #rescale the WAIFW if a specific R0 is specified
#             if (rescale.WAIFW & length(exper@R0)>0) {
#               #print("RESCALING!!!")
#               exper@trans@waifw <- scale.WAIFW.space(exper@R0,
#                                                      state,exper@trans@waifw,
#                                                      frequency.dep=exper@trans@frequency.dep,
#                                                      suscept.state=exper@trans@s.inds[1])
#             }
#
#             #get the number of time steps in the experiment
#             T <- round((exper@t.max-exper@t.min)/exper@step.size)+1
#
#             #get seasonal mult
#             if (!is.null(exper@season.obj)) {
#               mults <- get.seasonal.mult(exper@t0.doy/365+(1:T-1)*exper@step.size,
#                                          exper@season.obj)
#             } else {
#               mults <- rep(1,T)
#             }
#
#             #make a temporary transmission object
#             tmp.trans <- exper@trans
#
#             #hold the states as we walk through
#             rc <- matrix(ncol = T, nrow = nrow(state))
#             rc[,1] <- state
#
#             #need output vector for births over time
#             births.each.timestep <- growth.rate.each.timestep <- N0 <- matrix(NA, n.subpops, T)
#             births.each.timestep[,1] <- tmp.trans@birth.rate
#
#             #output vector for when SIAs were administered - default in this experiment is 0
#             sia.times <- routine.intro <- rep(0,T)
#
#
#             for (t in 2:T) {
#
#               print(t)
#
#               #scale the waifw by seasonality
#               if (!is.array(mults)){
#                 tmp.trans@waifw <- exper@trans@waifw*mults[t]
#               } else {
#                 tmp.trans@waifw <- exper@trans@waifw*mults[,,t]
#               }
#
#               #if exper@pop.rescale.each.timestep[t]==NaN then tmp.trans@pop.rescale==NaN and will not rescale
#               #otherwise if any number other than NaN it will rescale
#               #or if exper@pop.rescale.each.timestep not completed , then numeric(0) and any index [t] is NA
#               for (s in 1:n.subpops) {
#                 if (!is.na(exper@pop.rescale.each.timestep[s,t])) {
#                   last.index.tmp <- (s*exper@trans@n.epi.class*exper@trans@n.age.class)
#                   first.index.tmp <- last.index.tmp-(exper@trans@n.epi.class*exper@trans@n.age.class)+1
#                   state[first.index.tmp:last.index.tmp] <- exper@pop.rescale.each.timestep[s,t]*
#                     (state[first.index.tmp:last.index.tmp]/sum(state[first.index.tmp:last.index.tmp]))
#                 }
#               }
#
#               #put in the correct birth rate for that time-step, if it varies
#               for (s in 1:n.subpops) {
#                 if (ncol(exper@births.per.1000.each.timestep)>1) {
#                   last.index.tmp <- (s*exper@trans@n.epi.class*exper@trans@n.age.class)
#                   first.index.tmp <- last.index.tmp-(exper@trans@n.epi.class*exper@trans@n.age.class)+1
#                   tmp.trans@birth.rate[s] = (exper@births.per.1000.each.timestep[s,t]*
#                                                sum(state[first.index.tmp:last.index.tmp])/1000)
#                 } else {
#                   tmp.trans@birth.rate[s] <- exper@trans@birth.rate[s]
#                 }
#               }
#               births.each.timestep[,t] <-  tmp.trans@birth.rate
#
#               #put in time appropriate survival rate otherwise it uses original surv.matrix set up for trans object and keeps constant over time
#               if (!is.na(exper@surv.each.timestep[1,1])) {
#                 for (s in 1:n.subpops) {
#                   last.index.tmp1 <- (s*exper@trans@n.age.class)
#                   first.index.tmp1 <- last.index.tmp1-(exper@trans@n.age.class)+1
#                   tmp.age.surv.matrix <- ExtractAgeSpecificSurvivalMatrix(tmp.trans=tmp.trans, maternal.obj=exper@maternal.obj,
#                                                                           surv.at.timestep.t=exper@surv.each.timestep[first.index.tmp1:last.index.tmp1,t])
#
#                   last.index.tmp2 <- (s*exper@trans@n.epi.class*exper@trans@n.age.class)
#                   first.index.tmp2 <- last.index.tmp2-(exper@trans@n.epi.class*exper@trans@n.age.class)+1
#                   tmp.trans@age.surv.matrix[first.index.tmp2:last.index.tmp2,] <- tmp.age.surv.matrix
#                 }
#               }
#
#               #stow the previous population size
#               for (s in 1:n.subpops) {
#                 last.index.tmp <- (s*exper@trans@n.epi.class*exper@trans@n.age.class)
#                 first.index.tmp <- last.index.tmp-(exper@trans@n.epi.class*exper@trans@n.age.class)+1
#                 N0[s,t] <- sum(state[first.index.tmp:last.index.tmp])
#               }
#
#               #run experiment to get the next time step
#               state <- next.ID.state(state, tran=tmp.trans)
#
#               #stow the previous population size
#               for (s in 1:n.subpops) {
#                 last.index.tmp <- (s*exper@trans@n.epi.class*exper@trans@n.age.class)
#                 first.index.tmp <- last.index.tmp-(exper@trans@n.epi.class*exper@trans@n.age.class)+1
#                 NT <- sum(state[first.index.tmp:last.index.tmp])
#                 growth.rate.each.timestep[s,t] <- log(NT/N0[s,t]) #instantaneous biweekly growth rate
#               }
#
#               #print(t)
#               #print(dim(state))
#               rc [,t] <- state
#             }
#
#             result <- new("sim.results.MSIRV.update.demog.space",
#                           data=rc,
#                           m.inds = exper@trans@m.inds,
#                           s.inds = exper@trans@s.inds,
#                           i.inds = exper@trans@i.inds,
#                           r.inds = exper@trans@r.inds,
#                           v.inds = exper@trans@v.inds,
#                           t = exper@t.min+(1:T-1)* exper@step.size,
#                           age.class = exper@trans@age.class,
#                           births.each.timestep = births.each.timestep,
#                           growth.rate.each.timestep = growth.rate.each.timestep,
#                           routine.intro = routine.intro,
#                           sia.times = sia.times,
#                           n.subpops = exper@trans@n.subpops,
#                           subpop.class.label = exper@trans@subpop.class.label)
#
#             rc <- new("experiment.result",
#                       experiment.def = exper,
#                       result = result)
#
#
#             return(rc)
#           }
# )

# #Class for a spatial vaccination experiment updating demography and with changing vaccination coverage over time, with vaccination limitations
# setClass("experiment.updatedemog.vaccinationchange.vaccinationlimitations.spatial",
#          representation(trans = "ID.transition.MSIRV.space"),
#          contains=c("experiment.updatedemog.vaccinationchange.vaccinationlimitations"))

# # Run method for experiment.updatedemog.vaccinationchange.vaccinationlimitations.spatial objects
# setMethod("run",
#           "experiment.updatedemog.vaccinationchange.vaccinationlimitations.spatial",
#           function(exper, rescale.WAIFW=T, ...){
#             print("THIS EXPERIMENT RUN METHOD HAS NOT YET BEEN CODED")
#           }
# )


# setMethod("getAverageInfectionAge",
#           "sim.results.MSIRV",
#           function(sim.res) {
#             age.mids <- (sim.res@age.class +
#                            c(0,sim.res@age.class[2:length(sim.res@age.class)-1]))/2
#             tmp <- sim.res[sim.res@i.inds,]*age.mids
#
#             rc <- colSums(tmp)/
#               colSums(sim.res[sim.res@i.inds,])
#             return(rc)
#
#           })


# # next ID state for MSIRV.space class
# setMethod("next.ID.state",
#           c("ID.state.matrix", "ID.transition.MSIRV.space"),
#           function (state, tran) {
#
#             #get the transition matrix
#             tran.matrix <- tran@age.surv.matrix
#
#             #get a spatial index on same scale as epi indexes
#             index.loc <-  tran@subpop.class.label[tran@i.inds]
#             one.loc <- index.loc==1
#
#             #define the phi matrix - assume for now same matrix
#             for (n in 1:tran@n.subpops) {
#               this.loc <- index.loc==n
#
#               #The vaccination logic
#               age.spec.vacc.prob <- tran@vac.per@pvacc.in.age.class[this.loc]
#
#               #M->V transition
#               tran.matrix[tran@v.inds[this.loc],tran@m.inds[one.loc]] <-  tran@age.surv.matrix[tran@v.inds[this.loc],tran@m.inds[one.loc]] *
#                 age.spec.vacc.prob
#               #M->M transition
#               tran.matrix[tran@m.inds[this.loc],tran@m.inds[one.loc]] <-  tran@age.surv.matrix[tran@m.inds[this.loc],tran@m.inds[one.loc]] *
#                 (1-age.spec.vacc.prob)
#               #M->S transition
#               tran.matrix[tran@s.inds[this.loc],tran@m.inds[one.loc]] <-  tran@age.surv.matrix[tran@s.inds[this.loc],tran@m.inds[one.loc]] *
#                 (1-age.spec.vacc.prob)
#               #S->V transition
#               tran.matrix[tran@v.inds[this.loc],tran@s.inds[one.loc]] <-  tran@age.surv.matrix[tran@v.inds[this.loc],tran@s.inds[one.loc]] *
#                 age.spec.vacc.prob
#               #S->S transition
#               tran.matrix[tran@s.inds[this.loc],tran@s.inds[one.loc]] <-  tran@age.surv.matrix[tran@s.inds[this.loc],tran@s.inds[one.loc]] *
#                 (1-age.spec.vacc.prob)
#               #S->I transition
#               tran.matrix[tran@i.inds[this.loc],tran@s.inds[one.loc]] <-  tran@age.surv.matrix[tran@i.inds[this.loc],tran@s.inds[one.loc]] *
#                 (1-age.spec.vacc.prob)
#
#
#               #define denominator of phi matrix - preventing NAs by setting min to 1
#               if (tran@frequency.dep) {denom<-max(sum(state[tran@subpop.class.label==n]),1)} else {denom<-1}
#
#               phi <- (tran@waifw%*%(state[tran@i.inds[this.loc],]^tran@exponent))/denom
#               phi <- 1 - exp(-phi)
#               phi <- matrix(phi, nrow=state@n.age.class ,
#                             ncol=state@n.age.class)
#               phi <- t(phi)
#
#               #print(c("phi",range(phi)))
#
#               #make susceptible part of matrix
#               tran.matrix[tran@s.inds[this.loc], tran@s.inds[one.loc]] <-
#                 tran.matrix[tran@s.inds[this.loc], tran@s.inds[one.loc]] * (1-phi)
#
#               #make infected part of matrix
#               tran.matrix[tran@i.inds[this.loc], tran@s.inds[one.loc]] <-
#                 tran.matrix[tran@i.inds[this.loc], tran@s.inds[one.loc]] * (phi)
#
#               #parametric fit
#               #prop local susceptibility (age specific)
#               age.struct.here <- (state[tran@m.inds[this.loc]]+state[tran@s.inds[this.loc]]+
#                                     state[tran@i.inds[this.loc]]+state[tran@r.inds[this.loc]]+
#                                     state[tran@v.inds[this.loc]])
#               xtj <- (state[tran@s.inds[this.loc]]/pmax(age.struct.here,1))
#               xtj <- xtj/length(tran@introduction.rate[this.loc]) #rescale since have multiple age classes
#               #print(c("xtj",range(xtj)))
#               #print(range(age.struct.here))
#
#               #prop non-local infectious - sum over age classes - broken down by location
#               #since need to multiply by each coupling par
#               yt <- sapply(split(state[tran@i.inds[!this.loc]],index.loc[!this.loc]),sum)/
#                 pmax(sapply(split(state[tran@m.inds[!this.loc]]+state[tran@s.inds[!this.loc]]+
#                                     state[tran@i.inds[!this.loc]]+state[tran@r.inds[!this.loc]]+
#                                     state[tran@v.inds[!this.loc]],index.loc[!this.loc]),sum),1)
#
#               #ADD NEW GUYS TO THE INTRODUCTION.RATE FOR THIS SITE - at the moment distributed evenly over age
#               if (!tran@is.stochastic) {
#                 tran@introduction.rate[this.loc] <-  tran@introduction.rate[this.loc] +
#                   (1-exp(-sum(tran@coupling[n,-n]*yt)*xtj))
#               } else {
#                 immigs <- rbinom(length(xtj),1,pmax(1-exp(-sum(tran@coupling[n,-n]*yt)*xtj),0))
#                 #print(yt); print(xtj); print(tran@coupling[n,-n])
#                 #print(immigs)
#                 if (sum(immigs)>0) {
#                   tran@introduction.rate[this.loc] <-
#                     rbinom(length(tran@i.inds[this.loc]),1,tran@introduction.rate[this.loc])+immigs
#
#                 } else {
#                   tran@introduction.rate[this.loc] <-
#                     rbinom(length(tran@i.inds[this.loc]),1,tran@introduction.rate[this.loc])
#                 }
#               }
#             }
#
#             #no one stays infected...
#             tran.matrix[tran@i.inds, tran@i.inds[one.loc]] <- 0
#
#             #print(c(unique(tran.matrix)))
#
#             if (!tran@is.stochastic) {
#
#               #loop over sites and multiply matrices
#               for (n in 1:tran@n.subpops) {
#                 here <- which(tran@subpop.class.label==n,arr.ind=TRUE)
#                 state[here,] <- tran.matrix[here,]%*%state[here,]
#               }
#
#               #add in the births to 1,1 for now, assuming that is correct.
#               #might go back on this later.
#               state[!duplicated(tran@subpop.class.label),1] <-  state[!duplicated(tran@subpop.class.label),1] +
#                 tran@birth.rate #A BIT OF A HACK
#
#               #crude way of handling introductions and emigrations
#               state[tran@i.inds,1] <- state[tran@i.inds,1] + tran@introduction.rate
#
#             } else {
#
#               #mortality probability in each category of the n.age.class * no classes
#               mort <- 1-rep(tran@survival.rate, each=state@n.epi.class)
#
#
#               #loop over and distribute via a multinomial
#               newstate <- rep(0,1+length(tran.matrix[,1]))
#
#               for (n in 1:tran@n.subpops) {
#
#                 here <- which(tran@subpop.class.label==n,arr.ind=TRUE)
#                 here.next <- c(here,length(newstate))
#
#
#                 #  print(here)
#
#                 for (k in 1:length(tran.matrix[1,])) {
#
#                   #print("living, then dying")
#                   #print(sum(tran.matrix[here,k]))
#                   #print(mort[((n-1)*tran@n.age.class)+k])
#
#                   newstate[here.next] <- newstate[here.next] +
#                     rmultinom(1,state[here[k],1],c(tran.matrix[here,k],mort[((n-1)*tran@n.age.class)+k]))
#                 }
#               }
#
#               #print("prop dead")
#               #print(newstate[length(newstate)]/sum(newstate))
#
#               state[,1] <- newstate[1:length(tran.matrix[,1])]
#
#               #print(unique(state[,1]))
#
#               #print("before births")
#               #print(range(c(state@.Data)))
#
#               #add in the births to 1,1 for now, assuming that is correct.
#               state[!duplicated(tran@subpop.class.label),1] <-  state[!duplicated(tran@subpop.class.label),1] +
#                 rpois(tran@n.subpops,tran@birth.rate)
#
#
#               #print("after births")
#               #print(range(c(state@.Data)))
#
#               #print(unique(tran@introduction.rate))
#
#
#               #the stoch dist introduced above for efficiency
#               state[tran@i.inds,1] <- state[tran@i.inds,1] + tran@introduction.rate
#
#
#
#             }
#
#
#             return(state)
#           })


# Holds the results of a simulation with a MSIRV.space object
# setClass("sim.results.MSIRV.space",
#          representation(n.subpops = "numeric",
#                         subpop.class.label = "numeric"),
#          contains= "sim.results.MSIRV")


# Holds the results of a simulation with an experiment.updatedemog object
# setClass("sim.results.MSIRV.update.demog.space",
#          representation(births.each.timestep = "ANY", #the number of births per time step as output from simulation
#                         growth.rate.each.timestep = "ANY"
#          ),
#          contains="sim.results.MSIRV.space")

# # Holds the results output for a simulation with an experiment.updatedemog.vaccinationchange object
# setClass("sim.results.MSIRV.update.demog.vaccine.change.space",
#          representation(MR1.fail.each.timestep = "ANY", #the number of births per time step as output from simulation
#                         MR2.fail.each.timestep = "ANY",
#                         SIA.fail.each.timestep = "ANY"
#          ),
#          contains="sim.results.MSIRV.update.demog.space")


# #Plot a sim.results.MSIRV.space object
# #does a 4 panel plot that shows M+S,I,R+V and average age for each subpopulation
# #Parameters -
# #     x - the sim.results.MSIRV.space object
# #     y - ignored
# #
# setMethod("plot",
#           "sim.results.MSIRV.space",
#           function (x,y,from=0, to=max(x@t), low.age = 0,
#                     high.age = max(x@age.class), proportions = FALSE,
#                     plot.events = TRUE,...) {
#
#             result.all <- x
#             orig.t <- x@t
#
#             for (s in 1:result.all@n.subpops){
#
#               #limit the data to only the subpopulation
#               x <- result.all
#               here <- which(x@subpop.class.label==s,arr.ind=TRUE)
#               x@.Data <- x[here,x@t>=from & x@t<=to]
#
#               #limit the time to the time of interest
#               x@t <- x@t[x@t>=from & x@t<=to]
#
#               #limit to the age classes of interest and the subpopulation of interest
#               #the age classes will be the same regardless of the subpopulation, so I didn't specify the indeces for age.class vector
#               age.classes <- c()
#               for (i in 1:length(x@age.class)) {
#                 age.classes <- c(age.classes, rep(x@age.class[i],5))
#               }
#               #similarly the msirv indeces will be the same once we have shrunken the .Data object to the correct rows of the subpopulation
#               last.index <- length(x@age.class)
#               first.index <- last.index-(length(x@age.class))+1
#               #make everything be in the age classes of interest
#               x@.Data <- x[age.classes>=low.age & age.classes<=high.age,]
#               x@m.inds <- x@m.inds[(first.index:last.index)[x@age.class>=low.age & x@age.class<=high.age]]
#               x@s.inds <- x@s.inds[(first.index:last.index)[x@age.class>=low.age & x@age.class<=high.age]]
#               x@i.inds <- x@i.inds[(first.index:last.index)[x@age.class>=low.age & x@age.class<=high.age]]
#               x@r.inds <- x@r.inds[(first.index:last.index)[x@age.class>=low.age & x@age.class<=high.age]]
#               x@v.inds <- x@v.inds[(first.index:last.index)[x@age.class>=low.age & x@age.class<=high.age]]
#               x@age.class <- x@age.class[x@age.class>=low.age & x@age.class<=high.age]
#
#
#               tots <- getCompartmentTotals(x)
#
#               if (proportions) {
#                 tots <- toProportions(tots)
#               }
#
#
#
#               plt.events <- function() {
#                 if (!plot.events) return()
#
#                 sia.times <- orig.t[which(x@sia.times>0)]
#                 routine.times <- orig.t[which(x@routine.intro>0)]
#
#                 if (length(sia.times)>0) {
#                   for (i in 1:length(sia.times)) {
#                     lines(rep(sia.times[i],2),c(0, 10^10),
#                           col=rgb(t(col2rgb("cornflowerblue")/256), alpha=.45),
#                           lty=3)
#
#                   }
#                 }
#
#
#                 if (length(routine.times)>0) {
#                   for (i in 1:length(routine.times)) {
#                     lines(rep(routine.times[i],2),c(0, 10^10),
#                           col=rgb(t(col2rgb("coral")/256), alpha=.45),
#                           lty=3)
#                   }
#                 }
#
#               }
#
#               par(mfcol = c(4,1), mar=c(3,4,2,2), oma=c(0,0,2,0))
#               plot(tots@t, tots[tots@s.inds,]+tots[tots@m.inds,], type="b", xlab="t",
#                    ylab="M+S", cex=.75, ylim=c(0, max(tots[tots@s.inds,]+tots[tots@m.inds,])))
#               lines(tots@t, tots[tots@s.inds,], lty=2, col="green")
#               lines(tots@t, tots[tots@m.inds,], lty=3, col="blue")
#               plt.events()
#               plot(tots@t,tots[tots@i.inds,], type="b", xlab="t", ylab="I", cex=.75)
#               plt.events()
#               plot(tots@t,tots[tots@r.inds,]+tots[tots@v.inds,], type="b", xlab="t",
#                    ylab="R+V", cex=.75,ylim=c(0, max(tots[tots@r.inds,]+tots[tots@v.inds,])))
#               lines(tots@t, tots[tots@r.inds,], lty=2, col="green")
#               lines(tots@t, tots[tots@v.inds,], lty=3, col="blue")
#               plt.events()
#               plot(x@t,getAverageInfectionAge(x) , type="b", xlab="t", ylab="Age", cex=.75)
#               plt.events()
#               mtext(paste("Subpopulation",s), line=0, side=3, outer=TRUE, cex=2)
#             }
#
#           })


#'
#' #' Spatial wrapper Function to get number of individuals per age group and subpopulation
#' #'
#' #' @param trans transition object
#' #' @param state vector of state of the population at one time point
#' #'
#' #' @return matrix, size n.subpops (rows) by n.ages.class (columns) with the number of individuals in each age and space
#' #' @export
#' #'
#' #' @examples
#' space.wrapper.GetNumber.per.AgeGroup <- function(trans, state){
#'
#'   for (s in 1:trans@n.subpops){
#'     trans.tmp <- trans
#'     trans.tmp@epi.class <- trans.tmp@epi.class[trans@subpop.class.label==s]
#'     tmp <- GetNumber.per.AgeGroup(trans.tmp, state[trans@subpop.class.label==s,1])
#'     if (s==1) out <- tmp
#'     if (s!=1) out <- rbind(out, tmp)
#'   }
#'   return(out)
#'
#' }

#' #' Function to get spatial demography - this function is temporary and will need to be totally revamped
#' #'
#' #' @param uncode
#' #'
#' #' @return list, demographic data
#' #' @export
#' #'
#' #' @examples
#' space.getDemography <- function(uncode){
#'
#'   library(logspline)
#'   setOldClass("oldlogspline")
#'   library(KernSmooth)
#'   library(scales) #for alpha function
#'   library(zoo) #to fill in NAs from the dpt1 estimates
#'   library(survival) #for the PDF of age of vaccination
#'   library(countrycode) #to transfer b/w uncode and iso3codes and country names
#'   library(readxl)
#'   library(dplyr)
#'
#'   dyn.load("./source/MRModel-funcs.so") #run "R CMD SHLIB source/MRModel-funcs.c" in terminal to compile
#'   source("./source/build.R")
#'   source("./source/base.R")
#'   source("./source/user_interface.R")
#'   source("./source/who_un_inputs.R")
#'   source("./source/new_functions.R")
#'
#'   library(wpp2019) #UN and WHO data model inputs require wpp2019
#'   setup <- setupCountry.Dec2021(country="Zambia")
#'   year <- 1980
#'   t.max <- length(year:2100)
#'   generation.time <- 0.5
#'   age.classes <- c(1:240, seq(252,1212,12))
#'
#'   pop.total.1950.2100 <- rbind(setup$pop.total.1950.2100, setup$pop.total.1950.2100*0.5)
#'   pop.age.byageclasses.1950.2100 <- rbind(setup$pop.age.byageclasses.1950.2100, setup$pop.age.byageclasses.1950.2100*0.5)
#'   tfr.1950.2100 <- rbind(setup$tfr.1950.2100, setup$tfr.1950.2100)
#'   e0.1950.2100 <- rbind(setup$e0.1950.2100, setup$e0.1950.2100)
#'   asfr.1950.2100 <- rbind(setup$asfr.1950.2100, setup$asfr.1950.2100)
#'   repro.age.sex.dist.1950.2100 <- rbind(setup$repro.age.sex.dist.1950.2100, setup$repro.age.sex.dist.1950.2100)
#'   births.1950.2100 <- rbind(setup$births.1950.2100, setup$births.1950.2100)
#'   cbr.1950.2100 <- rbind(setup$cbr.1950.2100, setup$cbr.1950.2100)
#'   yr.agespecificbirths.per.1000 <- rbind(setup$yr.agespecificbirths.per.1000, setup$yr.agespecificbirths.per.1000)
#'   asdr.1950.2100.by5 <- rbind(setup$asdr.1950.2100.by5, setup$asdr.1950.2100.by5)
#'   asdr.object <- new("space.nMx",
#'                      rate.years = setup$asdr.object@rate.years,
#'                      rates = rbind(setup$asdr.object@rates, setup$asdr.object@rates),
#'                      mid.age = setup$asdr.object@mid.age,
#'                      n.subpops = 2)
#'
#'   return(list(pop.total.1950.2100=pop.total.1950.2100,
#'               pop.age.byageclasses.1950.2100=pop.age.byageclasses.1950.2100,
#'               tfr.1950.2100=tfr.1950.2100,
#'               e0.1950.2100=e0.1950.2100,
#'               asfr.1950.2100=asfr.1950.2100,
#'               repro.age.sex.dist.1950.2100=repro.age.sex.dist.1950.2100,
#'               births.1950.2100=births.1950.2100,
#'               cbr.1950.2100=cbr.1950.2100,
#'               asdr.1950.2100.by5=asdr.1950.2100.by5,
#'               asdr.object=asdr.object))
#' }
#'
#'










