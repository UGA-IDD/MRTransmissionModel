# all set class set methods functions

#### Vaccination #### ----------------------------------------------------------

#Class that holds the information on how many people should be
#vaccinated during each time step in each age class
#in the TSIR framework.
setClass("vacc.per.time.step.by.age",
         slots = list(pvacc.in.age.class = "numeric")
         # vector of how many are vaccinated in the age class per time step
)


## object for the CDF of vaccination by age (non scaled)
setClass("vaccine.cdf.byage",
         slots = list(cdf = "numeric", #vector of length ages
                      ages = "numeric" #vector of length cdf
         )
)

#Make the generic vaccine success object
setClass("vsucc")

#make a constants vaccine success object
setClass("vsucc.constant",
         slots = list(success.rate="numeric"),
         contains="vsucc")

#Make an object to hold a logistic function
#for vaccine success
setClass("vsucc.Logistic",
         slots = list(intercept = "numeric",
                      mo.eff = "numeric",
                      full.efficacy = "numeric" #the max reachable efficacy
         ),
         contains = "vsucc"
)

#### Transition Objects  #### --------------------------------------------------

# (object that holds everything needed to make transitions)

#Class definition for transition matrix. This matrix knows how to shift the
#ID.state.matrix from one time step to another, and contains all the fun stuff
#we need to do the transition (i.e., the matrix)
setClass("ID.transition.SIR",
         slots = list(n.epi.class="numeric", #number of different epi classes
                      epi.class="numeric", #the epi class of each row, between 1 and n.epi.classes
                      epi.class.label = "character", #epi class labels
                      s.inds = "numeric", #indexes of the susceptible states
                      i.inds = "numeric", #indexes of the infectious states
                      r.inds = "numeric", #indexes of the recovered states
                      n.age.class = "numeric", #the number of age classes
                      age.class = "numeric",#the actual age classes defined by upper age in class
                      aging.rate = "numeric", #percent that age out of each age class at each time step
                      survival.rate = "numeric", #the percent who survive in this age class at each time step
                      birth.rate = "numeric", #the birth rate for this transition
                      waifw = "matrix", #who acquires infection from whom matrix
                      age.surv.matrix = "matrix", #matrix governing age/survival
                      introduction.rate = "numeric", #an vector of possibly 0 case importation rates in each age class
                      exponent = "numeric",        #exponent on infecteds to correct for discretization
                      frequency.dep = "logical",	#true/false indicating if freq. dependence or density dep. is implemented
                      is.stochastic = "logical",   #true/false indicating if stochasticity desired
                      get.births = "function" #a function that takes state and birth rate; allowing births to be shaped by pop struct
         ))


#Class to represent transitions in the context of vaccination
setClass("ID.transition.SIR.vac",
         slots = list(vac.per="vacc.per.time.step.by.age"), #percent routinely vaccinated in each age class
         contains="ID.transition.SIR")

setClass("ID.transition.SIR.vac.SIA",
         slots = list(sia.vac = "numeric",  #number vaccinated in each age group by SIA on this transition
                      sia.vsucc = "vsucc" #object to determine the probability that an SIA vaccination is successful
         ),
         contains="ID.transition.SIR.vac")

#Class to represent transitions in the context of vaccination
#and maternal antibodies.
setClass("ID.transition.MSIRV",
         slots = list(m.inds = "numeric", #indexes of those maternally protected
                      v.inds = "numeric"), #index of those who have been vaccinated
         contains = "ID.transition.SIR.vac"
)

#Class to represent transitions in the context of vaccination
#and maternal antibodies.
setClass("ID.transition.MSIRV.SIA",
         slots = list(m.inds = "numeric", #indexes of those maternally protected
                      v.inds = "numeric"), #index of those who have been vaccinated
         contains = "ID.transition.SIR.vac.SIA"
)


#### Set cov estimates #### ----------------------------------------------------
setClass("cov.estimates",
         slots = list(sia.cov = "numeric",
                      sia.min.age = "numeric",
                      sia.max.age = "numeric",
                      MR1.cov = "numeric",
                      MR2.cov = "numeric",
                      inaccessible.prop = "numeric"
         ))


#### Demography #### -----------------------------------------------------------

#Class definition for transition matrix. This matrix knows how to shift the
#ID.state.matrix from one time step to another, and contains all the fun stuff
#we need to do the transition (i.e., the matrix)
setClass("nMx",
         slots = list(rates="data.frame", #age specific death rates per 1 from 1950 to 2100 by 5 year increments
                      mid.age="numeric" #mid-age associated with rates,
         ))


#### Experimental Infrastructure (naming experimental objects) #### ------------

#Class to define an experiment.
setClass("experiment",
         slots = list(name = "character", #a human readable name for the experiment
                      R0 = "numeric", #if included the R0 we scale to
                      state.t0 = "ID.state.matrix",#the initial state
                      t.min = "numeric", #the starting time for the experiment in years
                      t.max = "numeric", #the ending time for the experiment in years
                      t0.doy = "numeric", #the starting day of year, should default to 0
                      step.size = "numeric", #the step size in years, 1/no.gens.per.year
                      trans = "ID.transition.SIR", #defines the transitions for this experiment
                      season.obj = "seasonal.force", #seasonal forcing function
                      births.per.1000.each.timestep = "numeric", #the estimated crude birth rate per time step, length is T
                      description = "character" #describe the experiment
         ))

#Class defines an experiment with SIAs or opportunity for triggered SIAs
setClass("experiment.SIAtrigger",
         slots = list(trans = "ID.transition.SIR.vac", #the transitions for this experiment
                      sia.obj = "matrix", #defines the SIAs that will go on during this experiment
                      trigger.list="list" #defined the triggered vaccination events that occur in this experiment
         ),
         contains="experiment")

#Class for an experiment with birth and death rates updated each time step and population rescaling -xxamy
setClass("experiment.updatedemog",
         slots = list(trans = "ID.transition.SIR.vac", #the transitions for this experiment - xxamy added this
                      surv.each.timestep = "matrix", #the estimated age-specific survival rates over the length of the experiment, ncol is T
                      pop.rescale.each.timestep = "numeric", #a vector of 0's if not call for population rescale, if !=0 then rescale by this pop number
                      maternal.obj = "maternal.exp.decay"),  #maternal antibody object
         contains=c("experiment"))

#Class for an experiment with changing vaccination coverage over time
setClass("experiment.updatedemog.vaccinationchange",
         slots = list(time.specific.MR1cov = "vector", #vector of values of the coverage values
                      time.specific.MR2cov = "vector", #vector of values of the coverage values
                      time.specific.SIAcov = "vector", #vector of values of the coverage values
                      time.specific.min.age.MR1 = "vector",
                      time.specific.max.age.MR1 = "vector",
                      time.specific.min.age.MR2 = "vector",
                      time.specific.max.age.MR2 = "vector",
                      time.specific.min.age.SIA = "vector",
                      time.specific.max.age.SIA = "vector",
                      obj.vcdf.MR1 = "vaccine.cdf.byage",
                      obj.vcdf.MR2 = "vaccine.cdf.byage",
                      obj.prob.vsucc = "prob.vsucc.byage",
                      sia.timing.in.year = "vector"), #time of the SIA in years (e.g., 0.5 is July of the year)
         contains=c("experiment.updatedemog"))

#Class for an experiment
setClass("experiment.updatedemog.vaccinationchange.vaccinationlimitations",
         slots = list(MR1MR2correlation = "logical",
                      MR1SIAcorrelation = "logical",
                      MR2SIAcorrelation = "logical",
                      SIAinefficient = "logical",
                      SIAinacc = "logical",
                      prop.inacc = "numeric"),
         contains=c("experiment.updatedemog.vaccinationchange"))


#### Experiment Result Objects (naming experimental results objects) #### ------

# Class to hold SIR simulation results. Hold the results of the simulation as
# SIR compartment values at each time step.
setClass("sim.results.SIR",
         slots = list(s.inds = "numeric", #indexing of susceptibles
                      i.inds = "numeric", #infectious indexes
                      r.inds = "numeric", #recovered indexes
                      age.class = "numeric", #upper end of age classes
                      t = "numeric"#times when we have info in the matrix
         ),
         contains="matrix")

# Holds the results of a simulation with an MSIRV object
setClass("sim.results.MSIRV",
         slots = list(m.inds = "numeric",
                      v.inds = "numeric",
                      routine.intro = "numeric",
                      sia.times = "numeric"),
         contains= "sim.results.SIR")

# Holds the results of a simulation with an experiment.updatedemog object
setClass("sim.results.MSIRV.update.demog",
         slots = list(births.each.timestep = "numeric", #the number of births per time step as output from simulation
                      growth.rate.each.timestep = "numeric"
         ),
         contains="sim.results.MSIRV")

# Hols the results output for a simulation with an experiment.updatedemog.vaccinationchange object
setClass("sim.results.MSIRV.update.demog.vaccine.change",
         slots = list(MR1.fail.each.timestep = "numeric", #the number of births per time step as output from simulation
                      MR2.fail.each.timestep = "numeric",
                      SIA.fail.each.timestep = "numeric"
         ),
         contains="sim.results.MSIRV.update.demog")

# Holds the results of an experiment along with the experiment definition
setClass("experiment.result",
         slots = list(experiment.def = "experiment",
                      result = "sim.results.SIR"))


#### Maternal Immunity #### ----------------------------------------------------

#assumes that protection is the same as exponetial decay
setClass("maternal.exp.decay",
         slots = list(decay.rt = "numeric")#the exponential decay rate
)


#assume that there is some threshold age at which
#you are 100% protected if this age or younger
setClass("maternal.thresh",
         slots = list(thresh="numeric"))

# object for the prob of vaccine success by age
setClass("prob.vsucc.byage",
         slots = list(prob.vsucc = "numeric", #vector of length ages
                      ages = "numeric" #vector of length prob
         )
)

#### Seasonality #### ----------------------------------------------------------

#generic virtual class for seasonal forcing
setClass("seasonal.force")

setClass("seasonal.cosine",
         slots = list(amplitude = "numeric", #amplitude of seasonal forcing
                      offset = "numeric", #offset from t=0 for nadir
                      period = "numeric"), #number of peaks per year
         prototype = list(amplitude=.2, offset = 0, period = 1),
         contains="seasonal.force"
)

setClass("seasonal.periodic.bspline",
         slots = list(parameters = "numeric",
                      nbasis = "numeric",
                      degree = "numeric",
                      period = "numeric"
         ),
         prototype = list(parameters= c(0.5,1,1), nbasis=3, degree=3, period=1),
         contains="seasonal.force",
)

setClass("seasonal.age.specific",
         slots = list(amplitude = "numeric", #amplitude of seasonal forcing
                      offset = "numeric", #offset from t=0 for nadir
                      period = "numeric", #number of peaks per year
                      age.class = "numeric",#age classes in the WAIFW
                      chosen.age = "numeric"),#vector of length 2, indicating start and end of seasonally forced age class
         prototype = list(amplitude=.2, offset = 0, period = 1, age.class=c(1:60,60*12),chosen.age=c(50,60)),
         contains="seasonal.force",
)

#designed to take output from TSIR-type analysis
setClass("seasonal.piecewise.constant",
         slots = list(parameters = "numeric",
                      time.in.year = "numeric"),
         prototype = list(parameters = rep(1,24),time.in.year = (0:24)/24),
         contains="seasonal.force",
)

#designed to take output from TSIR-type analysis
setClass("seasonal.age.specific.piecewise.constant",
         slots = list(parameters = "numeric",
                      time.in.year = "numeric",
                      age.class = "numeric",#age classes in the WAIFW
                      chosen.age = "numeric"),#vector of length 2, indicating start and end of seasonally forced age class
         prototype = list(parameters = rep(1,24),time.in.year = (0:24)/24),
         contains="seasonal.force",
)

#### State Objects (matrix that holds the population state) #### ---------------

#Matrix that holds the population state in terms of infected, susceptible, recovered, etc.
setClass("ID.state.matrix",
         slots = list(n.epi.class = "numeric", #number of different
                      #epi classes
                      epi.class = "numeric", #the epid class of each row,
                      #between 1 and n.epi.class
                      epi.class.label = "character", #label for each epi class
                      n.age.class = "numeric"#the number of age classes
         ),
         contains="matrix")


##################
##### METHODS #### -------------------------------------------------------------
##################

#### Methods to Manipulate or Plot Experiment Output #### ----------------------

#Plot a sim.results.SIR object
#does a 4 panel plot that shows M+S,I,R+V and average age
#
#Parameters -
#     x - the sim.results.MSIRV object
#     y - ignored
#
setGeneric("plot",
           function(sim.res, ...) standardGeneric("plot"))


# plot - name of Generic
# sim.results.MSIRV - name of class (aka signature)

setMethod("plot",
          "sim.results.MSIRV",
          function (x,y,from=0, to=max(x@t), low.age = 0,
                    high.age = max(x@age.class), proportions = FALSE,
                    plot.events = TRUE,...) {
            par(mfcol = c(4,1))
            par(mar=c(3,4,2,2))

            orig.t <- x@t

            x@.Data <- x[,x@t>=from & x@t<=to]
            x@t <- x@t[x@t>=from & x@t<=to]

            age.classes <- c()
            for (i in 1:length(x@age.class)) {
              age.classes <- c(age.classes, rep(x@age.class[i],5))
            }


            #make everything be in the age classes of interest
            x@.Data <- x[age.classes>=low.age & age.classes<=high.age,]
            x@m.inds <- x@m.inds[x@age.class>=low.age & x@age.class<=high.age]
            x@s.inds <- x@s.inds[x@age.class>=low.age & x@age.class<=high.age]
            x@i.inds <- x@i.inds[x@age.class>=low.age & x@age.class<=high.age]
            x@r.inds <- x@r.inds[x@age.class>=low.age & x@age.class<=high.age]
            x@v.inds <- x@v.inds[x@age.class>=low.age & x@age.class<=high.age]
            x@age.class <- x@age.class[x@age.class>=low.age & x@age.class<=high.age]


            tots <- getCompartmentTotals(x)

            if (proportions) {
              tots <- toProportions(tots)
            }

            plt.events <- function() {
              if (!plot.events) return()

              sia.times <- orig.t[which(x@sia.times>0)]
              routine.times <- orig.t[which(x@routine.intro>0)]

              if (length(sia.times)>0) {
                for (i in 1:length(sia.times)) {
                  lines(rep(sia.times[i],2),c(0, 10^10),
                        col=rgb(t(col2rgb("cornflowerblue")/256), alpha=.45),
                        lty=3)

                }
              }


              if (length(routine.times)>0) {
                for (i in 1:length(routine.times)) {
                  lines(rep(routine.times[i],2),c(0, 10^10),
                        col=rgb(t(col2rgb("coral")/256), alpha=.45),
                        lty=3)
                }
              }

            }

            plot(tots@t, tots[tots@s.inds,]+tots[tots@m.inds,], type="b", xlab="t",
                 ylab="M+S", cex=.75, ylim=c(0, max(tots[tots@s.inds,]+tots[tots@m.inds,])))
            lines(tots@t, tots[tots@s.inds,], lty=2, col="green")
            lines(tots@t, tots[tots@m.inds,], lty=3, col="blue")
            plt.events()
            plot(tots@t,tots[tots@i.inds,], type="b", xlab="t", ylab="I", cex=.75)
            plt.events()
            plot(tots@t,tots[tots@r.inds,]+tots[tots@v.inds,], type="b", xlab="t",
                 ylab="R+V", cex=.75,ylim=c(0, max(tots[tots@r.inds,]+tots[tots@v.inds,])))
            lines(tots@t, tots[tots@r.inds,], lty=2, col="green")
            lines(tots@t, tots[tots@v.inds,], lty=3, col="blue")
            plt.events()
            plot(x@t,getAverageInfectionAge(x) , type="b", xlab="t", ylab="Age", cex=.75)
            plt.events()

          })

#Plot a sim.results.SIR object
#does a 4 panel plot that shows S,I,R and average age
#
#Parameters -
#     x - the sim.results.SIR object
#     y - ignored
#
setMethod("plot",
          "sim.results.SIR",
          function (x,y,...) {
            par(mfcol = c(4,1))
            par(mar=c(3,2,2,2))

            tots <- getCompartmentTotals(x)
            plot(tots@t,tots[1,], type="b", xlab="t",
                 ylab="S", cex=.75,...)
            plot(tots@t,tots[2,], type="b", xlab="t",
                 ylab="I", cex=.75, ...)
            plot(tots@t,tots[3,], type="b", xlab="t",
                 ylab="R", cex=.75, ...)
            plot(x@t,getAverageInfectionAge(x) , type="b", xlab="t",
                 ylab="Age", cex=.75, ...)
          })


#Get the total number of susceptible, infected, and
#recovered from and sim.results.SIR
#
#returns -
# a sim results object with no age classes
setGeneric("getCompartmentTotals",
           function(sim.res, ...) standardGeneric("getCompartmentTotals"))

# the name of the generic
# the name of the class (signature)
# the method itself (function)
setMethod("getCompartmentTotals",
          "sim.results.SIR",
          function (sim.res) {
            rc <- new("sim.results.SIR",
                      nrow = 3,
                      ncol = ncol(sim.res),
                      t = sim.res@t,
                      s.inds = 1,
                      i.inds = 2,
                      r.inds = 3,
                      age.class = max(sim.res@age.class))
            rc[1,] <- colSums(sim.res[sim.res@s.inds,])
            rc[2,] <- colSums(sim.res[sim.res@i.inds,])
            rc[3,] <- colSums(sim.res[sim.res@r.inds,])

            rownames(rc) <- c("S","I","R")
            return(rc)
          })
setMethod("getCompartmentTotals",
          "sim.results.MSIRV",
          function (sim.res) {
            rc <- new("sim.results.MSIRV",
                      nrow = 5,
                      ncol = ncol(sim.res),
                      t = sim.res@t,
                      m.inds = 1,
                      s.inds = 2,
                      i.inds = 3,
                      r.inds = 4,
                      v.inds = 5,
                      age.class = max(sim.res@age.class))
            #sia.times = sim.res@sia.times,
            #trigger.times = sim.res@trigger.times)
            rc[rc@m.inds,] <- colSums(sim.res[sim.res@m.inds,])
            rc[rc@s.inds,] <- colSums(sim.res[sim.res@s.inds,])
            rc[rc@i.inds,] <- colSums(sim.res[sim.res@i.inds,])
            rc[rc@r.inds,] <- colSums(sim.res[sim.res@r.inds,])
            rc[rc@v.inds,] <- colSums(sim.res[sim.res@v.inds,])

            rownames(rc) <- c("M","S","I","R","V")
            return(rc)
          })

#Method for getting the average age of infection at each time point.
#takes in a sim.results.SIR object and returns a numeric vector of
#the average age of infection at each time point
setGeneric("getAverageInfectionAge",
           function(sim.res, ...) standardGeneric("getAverageInfectionAge"))

setMethod("getAverageInfectionAge",
          "sim.results.SIR",
          function(sim.res) {
            age.mids <- (sim.res@age.class +
                           c(0,sim.res@age.class[2:length(sim.res@age.class)-1]))/2
            tmp <- sim.res[sim.res@i.inds,]*age.mids

            rc <- colSums(tmp)/
              colSums(sim.res[sim.res@i.inds,])
            return(rc)

          })


#Convert the results to proportions of each age class in each
#compartment, rather than the absolute numbers
setGeneric("toProportions", function(ob, ...) standardGeneric("toProportions"))

setMethod("toProportions",
          "sim.results.MSIRV",
          function (ob) {
            ac <- c()
            for (i in ob@age.class) {
              ac <- c(ac, rep(i,5))
            }



            #print(cbind(round(ob@t,2),t(round(ob[ac==50,],2)),colSums(ob[ac==50,])))
            for (i in ob@age.class) {
              div <- t(matrix(rep(colSums(ob[ac==i,]),5), ncol=5))
              ob[ac==i,] <- ob[ac==i,]/div
            }

            #print(cbind(round(ob@t,2),t(round(ob[ac==50,],2)),colSums(ob[ac==50,])))


            return(ob)
          })

#Function to get the cumulative prob of getting infected over a period of nweeks from a sim.results.MSIRV
# to account for the fact that can get have CRS if get sick during first trimester (1st three months)
#
#Parameters -
#       a sim.result.MSIRV
#       nbiweeks - chosen no biweeks of risk
#
#Returns -
#       an object with age classes set by fertility age classes
setGeneric("getCumRiskTrimester",
           function(sim.res,nbiweeks, ...) standardGeneric("getCumRiskTrimester"))

setMethod("getCumRiskTrimester",
          c("sim.results.MSIRV","numeric"),
          function (sim.res, nbiweeks, ...) {

            n.age.cats <- length(sim.res@age.class)
            ntimesteps <- length(sim.res@t)

            rc <- matrix(NA,n.age.cats,ntimesteps)

            for (t in 2:(ntimesteps-nbiweeks)) {
              cumvals <- rowSums(sim.res[sim.res@i.inds,t:(t+nbiweeks-1)])
              cumprobs <- cumvals/pmax(sim.res[sim.res@s.inds,t-1],1)
              rc[,t] <- cumprobs
              #print(cbind(cumprobs,sim.res[sim.res@s.inds,t-1]))
              #this is the probabilty of getting infected, because to be an I you
              #had to first be an S at time t-1, and then the sum of I's are the number of S's who became infected
            }

            rc[1,] <- 0
            rc[,1] <- rc[,2]
            #print("here")
            rc[,(ntimesteps-nbiweeks):ntimesteps] <- rc[,ntimesteps-nbiweeks]

            return(rc)
          })


#### Methods to get the percent immune by maternal antibodies. #### ------------

#Parameters -
#    age - the age we want the decay rate for
#    mat.obj - the object describing maternally acquired immunity
setGeneric("pmaternal",
           function(age, mat.obj) standardGeneric("pmaternal"))

setMethod("pmaternal",
          c("numeric", "maternal.exp.decay"),
          function (age, mat.obj) {
            return(exp(-age*mat.obj@decay.rt))
          })


setMethod("pmaternal",
          c("numeric", "maternal.thresh"),
          function (age, mat.obj) {
            rc <- sapply(age, function(age) {
              if(age<=mat.obj@thresh) {return(1)}
              return(0)
            })
            return(rc)
          })

#### Seasonality Methods #### --------------------------------------------------

#register a function for getting the seasonal adjustment for
#a time or series of times based on a seasonality object
#
#Parameter:
#     t.in.yr - time of the year in fraction of 1 year
#     season.obj - an object encoding the seasonality
#
#
#Returns:
#   the seasonal multiplier, or series thereof
setGeneric("get.seasonal.mult",
           function(t.in.yr, season.obj) standardGeneric("get.seasonal.mult"))


setMethod("get.seasonal.mult",
          c("numeric","seasonal.cosine"),
          function(t.in.yr, season.obj) {
            mult <- 1+ season.obj@amplitude *
              cos(2*season.obj@period*pi*(t.in.yr+season.obj@offset))
            return(mult)
          })

setMethod("get.seasonal.mult",
          c("numeric","seasonal.periodic.bspline"),
          function(t.in.yr, season.obj) {
            require(pomp)
            y <- periodic.bspline.basis(t.in.yr,nbasis=season.obj@nbasis,
                                        degree=season.obj@degree,period=season.obj@period)
            mult <- y%*%season.obj@parameters
            return(mult)
          })


setMethod("get.seasonal.mult",
          c("numeric","seasonal.piecewise.constant"),
          function(t.in.yr, season.obj) {
            mult <- season.obj@parameters[findInterval(t.in.yr%%1,season.obj@time.in.year)]
            return(mult)
          })


setMethod("get.seasonal.mult",
          c("numeric","seasonal.age.specific"),
          function(t.in.yr, season.obj) {
            mult.vec <- 1+ season.obj@amplitude *
              cos(2*season.obj@period*pi*(t.in.yr+season.obj@offset))

            nage <- length(season.obj@age.class)
            change.age <- which(season.obj@age.class>=season.obj@chosen.age[1] &
                                  season.obj@age.class<=season.obj@chosen.age[2])

            mult <- array(dim=c(nage,nage,length(t.in.yr)))
            for (j in 1:length(t.in.yr)) {
              mult[,,j] <- matrix(1,nage,nage)
              mult[,,j][change.age,change.age] <- mult.vec[j]
            }
            return(mult)
          })


setMethod("get.seasonal.mult",
          c("numeric","seasonal.age.specific.piecewise.constant"),
          function(t.in.yr, season.obj) {
            mult.vec <- season.obj@parameters[findInterval(t.in.yr%%1,season.obj@time.in.year)]

            nage <- length(season.obj@age.class)
            change.age <- which(season.obj@age.class>=season.obj@chosen.age[1] &
                                  season.obj@age.class<=season.obj@chosen.age[2])

            mult <- array(dim=c(nage,nage,length(t.in.yr)))
            for (j in 1:length(t.in.yr)) {
              mult[,,j] <- matrix(1,nage,nage)
              mult[,,j][change.age,change.age] <- mult.vec[j]
            }

            return(mult)
          })

#### Getting State at Next Time Step (next.ID.state methods) ####

#Method that takes an ID.state.matrix and an ID.transition object
#and create the next state of the ID.state.matrix
#
#Parameters -
#   state - current system state (ID.state.matrix)
#   tran - ID transition object that has all of the infor supervising the transition
#
#Returns -
#  an updated ID.state.matrix
setGeneric("next.ID.state",
           function(state, tran, ...) standardGeneric("next.ID.state"))

# default version of next.ID.state
setMethod("next.ID.state",
          c("ID.state.matrix", "ID.transition.SIR"),
          SIRIDTran)

# next ID state for ID.transition.SIR.vac class
setMethod("next.ID.state",
          c("ID.state.matrix", "ID.transition.SIR.vac"),
          function(state, tran) {

            #A bit of a hack, so we do not double deal with vaccination
            if (!is(tran, "ID.transition.MSIRV")) {
              delta.state <-  state[tran@s.inds]*
                (tran@vac.per@pvacc.in.age.class)

              state[tran@s.inds] <- state[tran@s.inds] - delta.state
              state[tran@r.inds] <- state[tran@r.inds] + delta.state

              mid.age.class <- (tran@age.class + c(0, tran@age.class[2:tran@n.age.class-1]))/2
              sia.sv.prob <-tran@sia.vac * pvacsuccess(mid.age.class, tran@sia.vsucc)

              delta.sia <- state[tran@s.inds]*sia.sv.prob
              state[tran@s.inds] <- state[tran@s.inds] - delta.sia
              state[tran@r.inds] <- state[tran@r.inds] + delta.sia
            }

            #rc <- callNextMethod(state, tran)
            rc <- SIRIDTran(state,tran) #avoid callNextMethod overhead
            return(rc)
          })

#next ID state for MSIRV class
setMethod("next.ID.state",
          c("ID.state.matrix", "ID.transition.MSIRV"),
          function (state, tran) {

            tran@age.surv.matrix <- .Call("update_age_surv_MSIRV",
                                          tran@age.surv.matrix,
                                          sz = nrow(tran@age.surv.matrix),
                                          tran@vac.per@pvacc.in.age.class,
                                          tran@v.inds,
                                          tran@m.inds,
                                          tran@s.inds,
                                          tran@i.inds)

            rc <- callNextMethod(state, tran)
            #print(rc[1:5,1])
            return(rc)
          })


#### Run Experiment #### -------------------------------------------------------

# Method to run an experiment that updates using birth rates and age-specific death rates (from the UNPD) - xxawinter
# Obviously considerable redundancy with the "run" method above; simply keeping this separate for now
# so as to not pollute work done with the method above.
#
#Parameters -
#   exper - the experiment object
#   rescaleWAIFW - should we rescale the WAIFW matrix? Don't do it if the experiment is not if a fullus susceptible state
#
#Returns -
#   an experiments result object

setGeneric("run", function(exper,...) standardGeneric("run"))

setMethod("run",
          "experiment.updatedemog",
          function(exper, rescale.WAIFW=T, ...) {

            #print(rescale.WAIFW)
            state <- exper@state.t0

            #rescale the WAIFW if a specific R0 is specified
            if (rescale.WAIFW & length(exper@R0)>0) {
              #print("RESCALING!!!")
              exper@trans@waifw <- scale.WAIFW(exper@R0,
                                               state,exper@trans@waifw,
                                               frequency.dep=exper@trans@frequency.dep,
                                               suscept.state=exper@trans@s.inds[1])
            }

            #get the number of time steps in the experiment
            T <- round((exper@t.max-exper@t.min)/exper@step.size)+1

            if (!is.null(exper@season.obj)) {
              mults <- get.seasonal.mult(exper@t0.doy/365+(1:T-1)*exper@step.size,
                                         exper@season.obj)
            } else {
              mults <- rep(1,T)
            }

            #make a temporary transmission object
            tmp.trans <- exper@trans

            #hold the states as we walk through
            rc <- matrix(ncol = T, nrow = nrow(state))
            rc[,1] <- state

            #need output vector for births overtime
            births.each.timestep <- growth.rate.each.timestep <- rep(NA, T)
            births.each.timestep[1] <- tmp.trans@birth.rate

            #output vector for when SIAs were administered - default in this experiment is 0
            sia.times <- routine.intro <- rep(0,T)

            for (t in 2:T) {

              #scale the waifw by seasonality
              if (!is.array(mults)){
                tmp.trans@waifw <- exper@trans@waifw*mults[t]
              } else {
                tmp.trans@waifw <- exper@trans@waifw*mults[,,t]
              }

              #if exper@pop.rescale.each.timestep[t]==NaN then tmp.trans@pop.rescale==NaN and will not rescale
              #otherwise if any number other than NaN it will rescale
              #or if exper@pop.rescale.each.timestep not completed , then numeric(0) and any index [t] is NA
              if (!is.na(exper@pop.rescale.each.timestep[t])) {state <- exper@pop.rescale.each.timestep[t] *(state/sum(state))}

              #put in the correct birth rate for that time-step, if it varies
              if (length(exper@births.per.1000.each.timestep)>1) {
                # Read 'birth rate'
                tmp.trans@birth.rate = (exper@births.per.1000.each.timestep[t]*sum(state)/1000)
              } else {
                tmp.trans@birth.rate <- exper@trans@birth.rate*sum(rc[,t-1])/sum(exper@state.t0)
              }
              births.each.timestep[t] <-  tmp.trans@birth.rate

              #put in time appropriate survival rate otherwise it uses original surv.matrix set up for trans object and keeps constant over time
              if (!is.na(exper@surv.each.timestep[1,1])) {
                tmp.trans@age.surv.matrix <- ExtractAgeSpecificSurvivalMatrix(tmp.trans=tmp.trans, maternal.obj=exper@maternal.obj,
                                                                              surv.at.timestep.t=exper@surv.each.timestep[,t])
              }

              #put in the correct SIAs
              #tmp.trans@sia.vac <- exper@sia.obj[,((t-1)%%ncol(exper@sia.obj))+1]
              #if (sum(tmp.trans@sia.vac)>0) {sia.times[t] <- 1}

              N0 <- sum(state)
              state <- next.ID.state(state, tmp.trans)
              NT <- sum(state)
              growth.rate.each.timestep[t] <- log(NT/N0) #instaneous biweekly growth rate

              #print(t)
              #print(dim(state))
              rc [,t] <- state

            }


            rc <- new("sim.results.MSIRV.update.demog",
                      data=rc,
                      m.inds = exper@trans@m.inds,
                      s.inds = exper@trans@s.inds,
                      i.inds = exper@trans@i.inds,
                      r.inds = exper@trans@r.inds,
                      v.inds = exper@trans@v.inds,
                      t = exper@t.min+(1:T-1)*exper@step.size,
                      age.class = exper@trans@age.class,
                      births.each.timestep = births.each.timestep,
                      growth.rate.each.timestep = growth.rate.each.timestep,
                      routine.intro = routine.intro,
                      sia.times = sia.times)

            rc <- new("experiment.result",
                      experiment.def = exper,
                      result = rc)


            return(rc)
          }
)

setMethod("run",
          "experiment.updatedemog.vaccinationchange",
          function(exper, rescale.WAIFW=T, ...) {

            #print(rescale.WAIFW)
            state <- exper@state.t0

            #rescale the WAIFW if a specific R0 is specified
            if (rescale.WAIFW & length(exper@R0)>0) {
              #print("RESCALING!!!")
              exper@trans@waifw <- scale.WAIFW(exper@R0,
                                               state,exper@trans@waifw,
                                               frequency.dep=exper@trans@frequency.dep,
                                               suscept.state=exper@trans@s.inds[1])
            }

            #get the number of time steps in the experiment
            T <- round((exper@t.max-exper@t.min)/exper@step.size)+1

            if (!is.null(exper@season.obj)) {
              mults <- get.seasonal.mult(exper@t0.doy/365+(1:T-1)*exper@step.size,
                                         exper@season.obj)
            } else {
              mults <- rep(1,T)
            }

            #make a temporary transmission object
            tmp.trans <- exper@trans

            #hold the states as we walk through
            rc <- matrix(ncol = T, nrow = nrow(state))
            rc[,1] <- state

            #need output vector for births over time
            births.each.timestep <- growth.rate.each.timestep <- rep(NA, T)
            births.each.timestep[1] <- tmp.trans@birth.rate

            #generate the age and time specific vaccination matrix
            routine <- get.routine.time.age.specific(time.step= exper@step.size*12,
                                                     age.classes=exper@trans@age.class,
                                                     time.specific.MR1cov=exper@time.specific.MR1cov,
                                                     age.min.MR1=exper@time.specific.min.age.MR1,
                                                     age.max.MR1=exper@time.specific.max.age.MR1,
                                                     time.specific.MR2cov=exper@time.specific.MR2cov,
                                                     age.min.MR2=exper@time.specific.min.age.MR2,
                                                     age.max.MR2=exper@time.specific.max.age.MR2,
                                                     obj.vcdf.MR1=exper@obj.vcdf.MR1,
                                                     obj.vcdf.MR2=exper@obj.vcdf.MR2,
                                                     obj.prob.vsucc=exper@obj.prob.vsucc,
                                                     MR1MR2correlation=F)
            routine.intro <- rep(0, T)
            if (any(exper@time.specific.MR1cov!=0)) routine.intro[min(which(exper@time.specific.MR1cov>0))*(1/exper@step.size)+1] <- 1
            if (any(exper@time.specific.MR2cov!=0)) routine.intro[min(which(exper@time.specific.MR2cov>0))*(1/exper@step.size)+1] <- 1
            index.routine.vacc <- c(1,rep(1:nrow(routine$age.time.specific.routine), each=(T-1)/exper@t.max))
            SIA <- get.sia.time.age.specific(age.classes=exper@trans@age.class,
                                             time.specific.SIAcov=exper@time.specific.SIAcov,
                                             age.min.sia=exper@time.specific.min.age.SIA,
                                             age.max.sia=exper@time.specific.max.age.SIA,
                                             obj.prob.vsucc=exper@obj.prob.vsucc)
            index.sia.vacc <- rep(NA,T)
            year.sia <- which(exper@time.specific.SIAcov!=0)
            index.sia.vacc[(year.sia-1)*(T-1)/exper@t.max + round((exper@sia.timing.in.year*(T-1)/exper@t.max))] <-  year.sia #minus 1 because adding the sia.timing
            sia.times <- ifelse(!is.na(index.sia.vacc), 1, 0)

            #need output vectors for primary vaccination failure over time
            MR1.fail.each.timestep <- MR2.fail.each.timestep <- SIA.fail.each.timestep <- rep(0, T)
            MR1.fail.each.timestep[1] <- routine$prop.fail.MR1[1]
            MR2.fail.each.timestep[1] <- routine$prop.fail.MR2[1]
            SIA.fail.each.timestep[1] <- 0

            for (t in 2:T) {

              #scale the waifw by seasonality
              if (!is.array(mults)){
                tmp.trans@waifw <- exper@trans@waifw*mults[t]
              } else {
                tmp.trans@waifw <- exper@trans@waifw*mults[,,t]
              }

              #if exper@pop.rescale.each.timestep[t]==NaN then tmp.trans@pop.rescale==NaN and will not rescale
              #otherwise if any number other than NaN it will rescale
              #or if exper@pop.rescale.each.timestep not completed , then numeric(0) and any index [t] is NA
              if (!is.na(exper@pop.rescale.each.timestep[t])) {state <- exper@pop.rescale.each.timestep[t] *(state/sum(state))}

              #put in the correct birth rate for that time-step, if it varies
              if (length(exper@births.per.1000.each.timestep)>1) {
                # Read 'birth rate'
                tmp.trans@birth.rate = (exper@births.per.1000.each.timestep[t]*sum(state)/1000)
              } else {
                tmp.trans@birth.rate <- exper@trans@birth.rate*sum(rc[,t-1])/sum(exper@state.t0)
              }
              births.each.timestep[t] <-  tmp.trans@birth.rate

              #put in time appropriate survival rate otherwise it uses original surv.matrix set up for trans object and keeps constant over time
              if (!is.na(exper@surv.each.timestep[1,1])) {
                tmp.trans@age.surv.matrix <- ExtractAgeSpecificSurvivalMatrix(tmp.trans=tmp.trans, maternal.obj=exper@maternal.obj,
                                                                              surv.at.timestep.t=exper@surv.each.timestep[,t])
              }

              ##put in correct vaccination coverage
              routine.vacc.prob <- routine$age.time.specific.routine
              sia.vacc.prob <- SIA$age.time.specific.SIA
              if (!is.na(index.sia.vacc[t])){ #if SIA
                tmp.trans@vac.per@pvacc.in.age.class <-
                  routine.vacc.prob[index.routine.vacc[t],] +
                  sia.vacc.prob[index.sia.vacc[t],] +
                  (routine.vacc.prob[index.routine.vacc[t],]*sia.vacc.prob[index.sia.vacc[t],])
                #stow output
                MR1.fail.each.timestep[t] <- routine$prop.fail.MR1[index.routine.vacc[t]]
                MR2.fail.each.timestep[t] <- routine$prop.fail.MR2[index.routine.vacc[t]]
                SIA.fail.each.timestep[t] <- SIA$prop.fail.SIA[index.sia.vacc[t]]

              } else { #if no SIA
                tmp.trans@vac.per@pvacc.in.age.class <- routine.vacc.prob[index.routine.vacc[t],]
                #stow output
                MR1.fail.each.timestep[t] <- routine$prop.fail.MR1[index.routine.vacc[t]]
                MR2.fail.each.timestep[t] <- routine$prop.fail.MR2[index.routine.vacc[t]]
              }


              N0 <- sum(state)
              state <- next.ID.state(state, tmp.trans)
              NT <- sum(state)
              growth.rate.each.timestep[t] <- log(NT/N0) #instantaneous biweekly growth rate

              #print(t)
              #print(dim(state))
              rc [,t] <- state

            }


            rc <- new("sim.results.MSIRV.update.demog.vaccine.change",
                      data=rc,
                      m.inds = exper@trans@m.inds,
                      s.inds = exper@trans@s.inds,
                      i.inds = exper@trans@i.inds,
                      r.inds = exper@trans@r.inds,
                      v.inds = exper@trans@v.inds,
                      t = exper@t.min+(1:T-1)* exper@step.size,
                      age.class = exper@trans@age.class,
                      births.each.timestep = births.each.timestep,
                      growth.rate.each.timestep = growth.rate.each.timestep,
                      MR1.fail.each.timestep = MR1.fail.each.timestep,
                      MR2.fail.each.timestep = MR2.fail.each.timestep,
                      SIA.fail.each.timestep = SIA.fail.each.timestep,
                      routine.intro = routine.intro,
                      sia.times = sia.times)

            rc <- new("experiment.result",
                      experiment.def = exper,
                      result = rc)


            return(rc)
          }
)

setMethod("run",
          "experiment.updatedemog.vaccinationchange.vaccinationlimitations",
          function(exper, rescale.WAIFW=T, ...) {

            #print(rescale.WAIFW)
            state <- exper@state.t0

            #rescale the WAIFW if a specific R0 is specified
            if (rescale.WAIFW & length(exper@R0)>0) {
              #print("RESCALING!!!")
              exper@trans@waifw <- scale.WAIFW(exper@R0,
                                               state,exper@trans@waifw,
                                               frequency.dep=exper@trans@frequency.dep,
                                               suscept.state=exper@trans@s.inds[1])
            }

            #get the number of time steps in the experiment
            T <- round((exper@t.max-exper@t.min)/exper@step.size)+1

            if (!is.null(exper@season.obj)) {
              mults <- get.seasonal.mult(exper@t0.doy/365+(1:T-1)*exper@step.size,
                                         exper@season.obj)
            } else {
              mults <- rep(1,T)
            }

            #make a temporary transmission object
            tmp.trans <- exper@trans

            #hold the states as we walk through
            rc <- matrix(ncol = T, nrow = nrow(state))
            rc[,1] <- state

            #need output vector for births over time
            births.each.timestep <- growth.rate.each.timestep <- rep(NA, T)
            births.each.timestep[1] <- tmp.trans@birth.rate

            #there are a couple options that are not yet coded - stop experiment if these are selected
            if ((exper@MR1SIAcorrelation & exper@SIAinacc) | (exper@MR2SIAcorrelation & exper@SIAinacc))  {
              stop("error: SIA limitation can either be correlation with a routine dose OR inaccessible population (with or without inefficiency), not both")
            }
            if ((exper@MR1SIAcorrelation & exper@SIAinefficient) | (exper@MR2SIAcorrelation & exper@SIAinefficient))  {
              stop("error: SIA limitation can either be correlation with a routine dose OR SIA inefficiency population (with or without inaccessible), not both")
            }
            if (!exper@MR1SIAcorrelation & exper@MR2SIAcorrelation){
              stop("error: MR1SIAcorrelation = FALSE and MR2SIAcorrelation= T option isn't yet coded or available")
            }

            #generate the age and time specific vaccination matrices for routine and SIAs
            routine <- get.routine.time.age.specific(time.step= exper@step.size*12,
                                                     age.classes=exper@trans@age.class,
                                                     time.specific.MR1cov=exper@time.specific.MR1cov,
                                                     age.min.MR1=exper@time.specific.min.age.MR1,
                                                     age.max.MR1=exper@time.specific.max.age.MR1,
                                                     time.specific.MR2cov=exper@time.specific.MR2cov,
                                                     age.min.MR2=exper@time.specific.min.age.MR2,
                                                     age.max.MR2=exper@time.specific.max.age.MR2,
                                                     obj.vcdf.MR1=exper@obj.vcdf.MR1,
                                                     obj.vcdf.MR2=exper@obj.vcdf.MR2,
                                                     obj.prob.vsucc=exper@obj.prob.vsucc,
                                                     MR1MR2correlation=exper@MR1MR2correlation)
            routine.intro <- rep(0, T)
            if (any(exper@time.specific.MR1cov!=0)) routine.intro[min(which(exper@time.specific.MR1cov>0))*(1/exper@step.size)+1] <- 1
            if (any(exper@time.specific.MR2cov!=0)) routine.intro[min(which(exper@time.specific.MR2cov>0))*(1/exper@step.size)+1] <- 1
            index.routine.vacc <- c(1,rep(1:nrow(routine$age.time.specific.routine), each=(T-1)/exper@t.max))
            SIA <- get.sia.time.age.specific(age.classes=exper@trans@age.class,
                                             time.specific.SIAcov=exper@time.specific.SIAcov,
                                             age.min.sia=exper@time.specific.min.age.SIA,
                                             age.max.sia=exper@time.specific.max.age.SIA,
                                             obj.prob.vsucc=exper@obj.prob.vsucc)
            index.sia.vacc <- rep(NA,T)
            year.sia <- which(exper@time.specific.SIAcov!=0)
            index.sia.vacc[(year.sia-1)*(T-1)/exper@t.max + round((exper@sia.timing.in.year*(T-1)/exper@t.max))] <-  year.sia #minus 1 because adding the sia.timing
            sia.times <- ifelse(!is.na(index.sia.vacc), 1, 0)

            #need output vectors for primary vaccination failure over time
            MR1.fail.each.timestep <- MR2.fail.each.timestep <- SIA.fail.each.timestep <- rep(0, T)
            MR1.fail.each.timestep[1] <- routine$prop.fail.MR1[1]
            MR2.fail.each.timestep[1] <- routine$prop.fail.MR2[1]
            SIA.fail.each.timestep[1] <- 0

            for (t in 2:T) {

              #scale the waifw by seasonality
              if (!is.array(mults)){
                tmp.trans@waifw <- exper@trans@waifw*mults[t]
              } else {
                tmp.trans@waifw <- exper@trans@waifw*mults[,,t]
              }

              #if exper@pop.rescale.each.timestep[t]==NaN then tmp.trans@pop.rescale==NaN and will not rescale
              #otherwise if any number other than NaN it will rescale
              #or if exper@pop.rescale.each.timestep not completed , then numeric(0) and any index [t] is NA
              if (!is.na(exper@pop.rescale.each.timestep[t])) {state <- exper@pop.rescale.each.timestep[t] *(state/sum(state))}

              #put in the correct birth rate for that time-step, if it varies
              if (length(exper@births.per.1000.each.timestep)>1) {
                # Read 'birth rate'
                tmp.trans@birth.rate = (exper@births.per.1000.each.timestep[t]*sum(state)/1000)
              } else {
                tmp.trans@birth.rate <- exper@trans@birth.rate*sum(rc[,t-1])/sum(exper@state.t0)
              }
              births.each.timestep[t] <-  tmp.trans@birth.rate

              #put in time appropriate survival rate otherwise it uses original surv.matrix set up for trans object and keeps constant over time
              if (!is.na(exper@surv.each.timestep[1,1])) {
                tmp.trans@age.surv.matrix <- ExtractAgeSpecificSurvivalMatrix(tmp.trans=tmp.trans, maternal.obj=exper@maternal.obj,
                                                                              surv.at.timestep.t=exper@surv.each.timestep[,t])
              }

              ##put in correct vaccination coverage
              tmp.trans@vac.per@pvacc.in.age.class <- rep(0, exper@trans@n.age.class) #start with clean slate each time step
              routine.vacc.prob <- routine$age.time.specific.routine
              sia.vacc.prob <- SIA$age.time.specific.SIA

              if (!is.na(index.sia.vacc[t])){ #if SIA

                #IF SIA independent MR1 or MR2
                if (!exper@MR1SIAcorrelation & !exper@MR2SIAcorrelation & !exper@SIAinacc & !exper@SIAinefficient){
                  tmp.trans@vac.per@pvacc.in.age.class <-
                    routine.vacc.prob[index.routine.vacc[t],] +
                    sia.vacc.prob[index.sia.vacc[t],] +
                    (routine.vacc.prob[index.routine.vacc[t],]*sia.vacc.prob[index.sia.vacc[t],])
                }

                #IF SIA and MR1 correlated, SIA & MR2 NOT correlated
                if (exper@MR1SIAcorrelation & !exper@MR2SIAcorrelation & !exper@SIAinacc & !exper@SIAinefficient) {
                  MR1.age.range <- exper@time.specific.min.age.MR1[index.routine.vacc[t]]:exper@time.specific.max.age.MR1[index.routine.vacc[t]]
                  MR2.age.range <- exper@time.specific.min.age.MR2[index.routine.vacc[t]]:exper@time.specific.max.age.MR2[index.routine.vacc[t]]
                  if (exper@time.specific.SIAcov[index.sia.vacc[t]]<exper@time.specific.MR1cov[index.routine.vacc[t]]){
                    #non MR1 or MR2 ages
                    tmp.trans@vac.per@pvacc.in.age.class <- sia.vacc.prob[index.sia.vacc[t],]
                    #MR1 ages
                    tmp.trans@vac.per@pvacc.in.age.class[MR1.age.range] <-
                      routine.vacc.prob[index.routine.vacc[t],MR1.age.range] + routine$one.minus.ve1[index.routine.vacc[t]]*sia.vacc.prob[index.sia.vacc[t],MR1.age.range]
                    #MR2 ages
                    tmp.trans@vac.per@pvacc.in.age.class[MR2.age.range] <-
                      routine.vacc.prob[index.routine.vacc[t],MR2.age.range] +
                      sia.vacc.prob[index.sia.vacc[t],MR2.age.range] +
                      (routine.vacc.prob[index.routine.vacc[t],MR2.age.range]*sia.vacc.prob[index.sia.vacc[t],MR2.age.range])
                  }
                  if (exper@time.specific.SIAcov[index.sia.vacc[t]]>exper@time.specific.MR1cov[index.routine.vacc[t]]){
                    #non MR1 or MR2 ages
                    tmp.trans@vac.per@pvacc.in.age.class <- sia.vacc.prob[index.sia.vacc[t],]
                    #MR1 ages
                    tmp.trans@vac.per@pvacc.in.age.class[MR1.age.range] <-
                      routine.vacc.prob[index.routine.vacc[t],MR1.age.range] + routine$prop.fail.MR1.byage[index.routine.vacc[t],MR1.age.range] +
                      (sia.vacc.prob[index.sia.vacc[t],MR1.age.range] - routine.vacc.prob[index.routine.vacc[t],MR1.age.range])
                    #MR2 ages
                    tmp.trans@vac.per@pvacc.in.age.class[MR2.age.range] <-
                      routine.vacc.prob[index.routine.vacc[t],MR2.age.range] +
                      sia.vacc.prob[index.sia.vacc[t],MR2.age.range] +
                      (routine.vacc.prob[index.routine.vacc[t],MR2.age.range]*sia.vacc.prob[index.sia.vacc[t],MR2.age.range])
                  }
                }

                #IF SIA and MR1 & SIA and MR2 correlated
                if (exper@MR1SIAcorrelation & exper@MR2SIAcorrelation & !exper@SIAinacc & !exper@SIAinefficient) {
                  MR1.age.range <- exper@time.specific.min.age.MR1[index.routine.vacc[t]]:exper@time.specific.max.age.MR1[index.routine.vacc[t]]
                  MR2.age.range <- exper@time.specific.min.age.MR2[index.routine.vacc[t]]:exper@time.specific.max.age.MR2[index.routine.vacc[t]]
                  if (exper@time.specific.SIAcov[index.sia.vacc[t]]<exper@time.specific.MR1cov[index.routine.vacc[t]]){
                    #non MR1 or MR2 ages
                    tmp.trans@vac.per@pvacc.in.age.class <- sia.vacc.prob[index.sia.vacc[t],]
                    #MR1 ages
                    tmp.trans@vac.per@pvacc.in.age.class[MR1.age.range] <-
                      routine.vacc.prob[index.routine.vacc[t],MR1.age.range] + routine$one.minus.ve1[index.routine.vacc[t]]*sia.vacc.prob[index.sia.vacc[t],MR1.age.range]
                    #MR2 ages
                    tmp.trans@vac.per@pvacc.in.age.class[MR2.age.range] <-
                      routine.vacc.prob[index.routine.vacc[t],MR2.age.range] + routine$one.minus.ve1[index.routine.vacc[t]]*sia.vacc.prob[index.sia.vacc[t],MR2.age.range]
                  }
                  if (exper@time.specific.SIAcov[index.sia.vacc[t]]>exper@time.specific.MR1cov[index.routine.vacc[t]]){
                    #non MR1 or MR2 ages
                    tmp.trans@vac.per@pvacc.in.age.class <- sia.vacc.prob[index.sia.vacc[t],]
                    #MR1 ages
                    tmp.trans@vac.per@pvacc.in.age.class[MR1.age.range] <-
                      routine.vacc.prob[index.routine.vacc[t],MR1.age.range] + routine$prop.fail.MR1.byage[index.routine.vacc[t],MR1.age.range] +
                      (sia.vacc.prob[index.sia.vacc[t],MR1.age.range] - routine.vacc.prob[index.routine.vacc[t],MR1.age.range])
                    #MR2 ages
                    tmp.trans@vac.per@pvacc.in.age.class[MR2.age.range] <-
                      routine.vacc.prob[index.routine.vacc[t],MR2.age.range] + routine$prop.fail.MR1.byage[index.routine.vacc[t],MR2.age.range] +
                      (sia.vacc.prob[index.sia.vacc[t],MR2.age.range] - routine.vacc.prob[index.routine.vacc[t],MR2.age.range])
                  }
                }

                #If inaccessible population, but not SIA inefficiency
                if (!exper@MR1SIAcorrelation & !exper@MR2SIAcorrelation & exper@SIAinacc & !exper@SIAinefficient) {
                  tmp.trans@vac.per@pvacc.in.age.class <- 1-(exper@prop.inacc[t]+
                                                               (1-exper@prop.inacc[t])*(1-(routine.vacc.prob[index.routine.vacc[t],] +
                                                                                             sia.vacc.prob[index.sia.vacc[t],] +
                                                                                             (routine.vacc.prob[index.routine.vacc[t],]*sia.vacc.prob[index.sia.vacc[t],]))))
                  #tmp.trans@vac.per@pvacc.in.age.class <- 1-(exper@prop.inacc[t]+(1-exper@prop.inacc[t])*(1-routine.vacc.prob[index.routine.vacc[t],])*(1-sia.vacc.prob[index.sia.vacc[t],]))
                }

                #If inaccessible population AND SIA inefficiency
                if (!exper@MR1SIAcorrelation & !exper@MR2SIAcorrelation & exper@SIAinacc & exper@SIAinefficient) {
                  # with SIA inaccessible & with SIA inefficiency -- VIMC 2017-2021 versions - xxamy
                  z <- routine.vacc.prob[index.routine.vacc[t],] #prob. successful routine given access
                  ro <- (1-exper@prop.inacc[t]) #prob. accessible
                  m <- (1-(sia.vacc.prob[index.sia.vacc[t],]*(1/ro)*(1-0.1)))^(1/(1-0.1)) #prob. not successful campaign vaccination given access, and inefficiency of 0.1
                  m[is.na(m)] <- 0 #if NaN then negative number was raised, which means everyone in accessible population was vaccinated by campaign, therefore prop. not vaccinated =0
                  age.spec.vacc.prob <- 1-((1-ro)+(ro*m)) #ages over 36 months, where only eligible for campaign
                  max.age.routine <- exper@time.specific.max.age.MR2[index.routine.vacc[t]]
                  age.spec.vacc.prob[1:max.age.routine] <- 1-((1-ro)+(ro*m*(1-(z/ro))))[1:max.age.routine] #ages 1:36 months where routine takes place
                  #if(sum(age.spec.vacc.prob)>0) print(age.spec.vacc.prob[1:max.age.routine])
                  tmp.trans@vac.per@pvacc.in.age.class <- age.spec.vacc.prob
                }

                #stow output
                MR1.fail.each.timestep[t] <- routine$prop.fail.MR1[index.routine.vacc[t]]
                MR2.fail.each.timestep[t] <- routine$prop.fail.MR2[index.routine.vacc[t]]
                SIA.fail.each.timestep[t] <- SIA$prop.fail.SIA[index.sia.vacc[t]]

              } else { #if no SIA
                tmp.trans@vac.per@pvacc.in.age.class <- routine.vacc.prob[index.routine.vacc[t],]
                #stow output
                MR1.fail.each.timestep[t] <- routine$prop.fail.MR1[index.routine.vacc[t]]
                MR2.fail.each.timestep[t] <- routine$prop.fail.MR2[index.routine.vacc[t]]
              }

              N0 <- sum(state)
              state <- next.ID.state(state, tmp.trans)
              NT <- sum(state)
              growth.rate.each.timestep[t] <- log(NT/N0) #instantaneous biweekly growth rate

              #print(t)
              #print(dim(state))
              rc [,t] <- state

            }


            rc <- new("sim.results.MSIRV.update.demog.vaccine.change",
                      data=rc,
                      m.inds = exper@trans@m.inds,
                      s.inds = exper@trans@s.inds,
                      i.inds = exper@trans@i.inds,
                      r.inds = exper@trans@r.inds,
                      v.inds = exper@trans@v.inds,
                      t = exper@t.min+(1:T-1)* exper@step.size,
                      age.class = exper@trans@age.class,
                      births.each.timestep = births.each.timestep,
                      growth.rate.each.timestep = growth.rate.each.timestep,
                      MR1.fail.each.timestep = MR1.fail.each.timestep,
                      MR2.fail.each.timestep = MR2.fail.each.timestep,
                      SIA.fail.each.timestep = SIA.fail.each.timestep,
                      routine.intro = routine.intro,
                      sia.times = sia.times)

            rc <- new("experiment.result",
                      experiment.def = exper,
                      result = rc)


            return(rc)
          }
)


#### Vaccine Success #### ------------------------------------------------------

#register the pvaccsuccess generic
#
#Parameters -
#   age - the age in months we are checking
#   vs.obj - the vaccince success object
#
#Returns -
#   the probability of successful vaccination at that
#   age
setGeneric("pvacsuccess",
           function(age, vsucc.obj) standardGeneric("pvacsuccess"))

#default method for vaccine success, just returns 1 (vaccine is always successful)
setMethod("pvacsuccess",
          c("numeric","vsucc.constant"),
          function (age, vsucc.obj) {

            prob.vsucc <- new("prob.vsucc.byage",
                              prob.vsucc=rep(vsucc.obj@success.rate, length(age)),
                              ages=age)

            return(prob.vsucc)
          })

#set the method for the vsucc logistic
setMethod("pvacsuccess",
          c("numeric","vsucc.Logistic"),
          function (age, vsucc.obj) {
            rc <- plogis(vsucc.obj@mo.eff*age,
                         -vsucc.obj@intercept)
            rc[rc>vsucc.obj@full.efficacy] <- vsucc.obj@full.efficacy

            prob.vsucc <- new("prob.vsucc.byage",
                              prob.vsucc=rc,
                              ages=age)
            return(prob.vsucc)
          }
)




