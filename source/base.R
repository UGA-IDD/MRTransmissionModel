## Base Functions for the MR Model
## Authors: Amy Winter & Justin Lessler & Jess Metcalf

#### State Objects (matrix that holds the population state) ####

#Matrix that holds the population state in terms of infected, susceptible, recovered, etc.
setClass("ID.state.matrix",
         representation(n.epi.class = "numeric", #number of different
                        #epi classes
                        epi.class = "numeric", #the epid class of each row,
                        #between 1 and n.epi.class
                        epi.class.label = "character", #label for each epi class
                        n.age.class = "numeric"#the number of age classes
         ),
         contains="matrix")

#Function creates an ID state matrix with the given number of
#epidemiologic categories, age classes and labels
#
#Parameters -
#   n.age.class - the number of age classes
#   n.epi.class - the number of states in the model
#   epi.class.label - the names of the states
#   value - starting values, defauls to all 0
#
#Returns -
#  an ID.state.matrix object
create.ID.state.matrix <- function(n.age.class,
                                   n.epi.class = 3,
                                   epi.class.label = c("S","I","R"),
                                   value = rep(0,n.age.class*n.epi.class)) {
  rc <- new("ID.state.matrix",
            nrow = n.age.class*n.epi.class,
            ncol = 1,
            n.epi.class = n.epi.class,
            #repeat epi classes for each age
            epi.class = rep(1:n.epi.class, n.age.class),
            epi.class.label = epi.class.label,
            n.age.class = n.age.class
  )
  
  rc[,1] <- value
  return(rc)
}

#### Transition Objects (object that holds everything needed to make transitions) ####

#Class definition for transition matrix. This matrix knows how to shift the
#ID.state.matrix from one time step to another, and contains all the fun stuff
#we need to do the transition (i.e., the matrix)
setClass("ID.transition.SIR",
         representation(n.epi.class="numeric", #number of different epi classes
                        epi.class="numeric", #the epi class of each row, between 1 and n.epi.classes
                        epi.class.label = "character", #epi class labels
                        s.inds = "numeric", #indexes of the susceptible states
                        i.inds = "numeric", #indexes of the infectious states
                        r.inds = "numeric", #indexes of the recovered states
                        n.age.class = "numeric", #the number of age classes
                        age.class = "numeric",#the actual age classes defined by upper age in class
                        aging.rate = "numeric", #percent that age out of each age class at each time step
                        survival.rate = "ANY", #the percent who survive in this age class at each time step
                        birth.rate = "ANY", #the birth rate for this transition
                        waifw = "matrix", #who acquires infection from whom matrix
                        age.surv.matrix = "matrix", #matrix governing age/survival
                        introduction.rate = "numeric", #an vector of possibly 0 case importation rates in each age class
                        exponent = "numeric",        #exponent on infecteds to correct for discretization
                        frequency.dep = "logical",	#true/false indicating if freq. dependence or density dep. is implemented
                        is.stochastic = "logical",   #true/false indicating if stochasticity desired
                        get.births = "function" #a function that takes state and birth rate; allowing births to be shaped by pop struct
         ))

#Simple method to create an ID.transition.matrix
#
#Parameters -
#   n.age.class - the number of age classes
#   n.epi.class - the number of states in the model
#   aging.rate  - the rate which people move to the next age class
#   survival.rate - the percent who survive at each time step
#   waifw - whom acquires infection from whom matrix
#   age.class - the actual age classes defined by upper age in class
create.ID.transition.SIR <- function(n.age.class,
                                     aging.rate,
                                     survival.rate,
                                     waifw,
                                     age.class = numeric(),
                                     n.epi.class = 3,
                                     birth.rate = 0,
                                     introduction.rate = 0,
                                     exponent = 1,
                                     frequency.dep=F,
                                     is.stochastic=F,
                                     get.births = function(pop,tran){return(c(tran@birth.rate,rep(0,length(pop)-1)))}) {
  rc <- new("ID.transition.SIR",
            n.epi.class = n.epi.class,
            epi.class = rep(1:n.epi.class, n.age.class),
            epi.class.label = c("S","I","R"),
            n.age.class = n.age.class,
            age.class = age.class,
            aging.rate = aging.rate,
            survival.rate = survival.rate,
            birth.rate = birth.rate,
            waifw = waifw,
            introduction.rate = introduction.rate,
            exponent = exponent,
            frequency.dep = frequency.dep,
            is.stochastic = is.stochastic,
            get.births = get.births)
  
  rc@s.inds <- which(rc@epi.class==1)
  rc@i.inds <- which(rc@epi.class==2)
  rc@r.inds <- which(rc@epi.class==3)
  
  rc@age.surv.matrix <- matrix(0,nrow = n.age.class*n.epi.class,
                               ncol = n.age.class*n.epi.class)
  
  
  #for zeroing out all the appropriate transitions
  template.mtrx <- matrix(c(1,1,0,0,1,1,0,0,1), nrow=3, ncol=3)
  
  #for each age class
  for (i in 1:n.age.class) {
    #first fill in the diagonal matrix
    inds <- (i-1)*n.epi.class+(1:n.epi.class)
    rc@age.surv.matrix[inds,inds] <- (1-aging.rate[i])*survival.rate[i] *
      template.mtrx
    
    #now fill in the off diagnal matrix
    if (i!=n.age.class) {
      #first fill in the diagonal matrix
      inds2 <- (i)*n.epi.class+(1:n.epi.class)
      rc@age.surv.matrix[inds2,inds] <- aging.rate[i]*survival.rate[i] *
        template.mtrx
    }
  }
  
  return(rc)
}


#Class to represent transitions in the context of vaccination
setClass("ID.transition.SIR.vac",
         representation(vac.per="vacc.per.time.step.by.age"), #percent routinely vaccinated in each age class
         contains="ID.transition.SIR")


#Method to create an ID.transtion.SIR.vac object
create.ID.transition.SIR.vac <- function(n.age.class,
                                         aging.rate,
                                         survival.rate,
                                         waifw,
                                         routine.vac = 0, #xxamy
                                         routine.vac.age.index = 12, #xxamy age index based on age classes
                                         time.step = 1/2, #time step in months == generation time
                                         age.class = numeric(),
                                         n.epi.class = 3,
                                         birth.rate = 0,
                                         introduction.rate = 0,
                                         exponent = 1,
                                         frequency.dep=F,
                                         is.stochastic=F,
                                         get.births = function(pop,tran){return(c(tran@birth.rate,rep(0,length(pop)-1)))}) {
  
  tmp <- create.ID.transition.SIR(n.age.class = n.age.class,
                                  aging.rate = aging.rate,
                                  survival.rate = survival.rate,
                                  waifw = waifw,
                                  age.class = age.class,
                                  n.epi.class = n.epi.class,
                                  birth.rate = birth.rate,
                                  introduction.rate = introduction.rate,
                                  exponent = exponent,
                                  frequency.dep=frequency.dep,
                                  is.stochastic=is.stochastic,
                                  get.births = get.births)
  
  #get the %vaccinated in each age group per time step
  #time step must be in years!
  #vac.per <- create.vacc.per.time.step.by.age(vpdf, age.class, time.step) #xxamy
  vac.per <- new("vacc.per.time.step.by.age", #xxamy
                 pvacc.in.age.class = numeric(length=n.age.class))
  vac.per@pvacc.in.age.class[routine.vac.age.index] <- routine.vac #xxamy
  
  
  rc <- new("ID.transition.SIR.vac",
            n.epi.class = tmp@n.epi.class,
            epi.class = tmp@epi.class,
            s.inds = tmp@s.inds,
            i.inds = tmp@i.inds,
            r.inds = tmp@r.inds,
            epi.class.label = c("S","I","R"),
            n.age.class = tmp@n.age.class,
            aging.rate = tmp@aging.rate,
            survival.rate = tmp@survival.rate,
            birth.rate = tmp@birth.rate,
            waifw = tmp@waifw,
            age.surv.matrix = tmp@age.surv.matrix,
            vac.per = vac.per,
            introduction.rate = introduction.rate,
            exponent = exponent,
            age.class = tmp@age.class,
            frequency.dep=frequency.dep,
            is.stochastic=is.stochastic,
            get.births = get.births
  )
  
  return(rc)
}

setClass("ID.transition.SIR.vac.SIA",
         representation(sia.vac = "numeric",  #number vaccinated in each age group by SIA on this transition
                        sia.vsucc = "vsucc" #object to determine the probability that an SIA vaccination is successful
         ),
         contains="ID.transition.SIR.vac")


#Method to create an ID.transtion.SIR.vac object
create.ID.transition.SIR.vac.SIA <- function(n.age.class,
                                             aging.rate,
                                             survival.rate,
                                             waifw,
                                             routine.vac = 0, #xxamy
                                             routine.vac.age.index = 12, #xxamy age index based on age classes
                                             time.step = 1/2, #time step in months == generation time
                                             age.class = numeric(),
                                             n.epi.class = 3,
                                             birth.rate = 0,
                                             sia.vac = 0,
                                             sia.vsucc = new("vsucc.constant", success.rate=1),
                                             introduction.rate = 0,
                                             exponent = 1,
                                             frequency.dep=F,
                                             is.stochastic=F,
                                             get.births = function(pop,tran){return(c(tran@birth.rate,rep(0,length(pop)-1)))}) {
  
  tmp <- create.ID.transition.SIR.vac(n.age.class = n.age.class,
                                      aging.rate = aging.rate,
                                      survival.rate = survival.rate,
                                      waifw = waifw,
                                      age.class = age.class,
                                      n.epi.class = n.epi.class,
                                      birth.rate = birth.rate,
                                      introduction.rate = introduction.rate,
                                      exponent = exponent,
                                      frequency.dep=frequency.dep,
                                      is.stochastic=is.stochastic,
                                      get.births = get.births)
  
  #get the %vaccinated in each age group per time step
  #time step must be in years!
  #vac.per <- create.vacc.per.time.step.by.age(vpdf, age.class, time.step) #xxamy
  vac.per <- new("vacc.per.time.step.by.age", #xxamy
                 pvacc.in.age.class = numeric(length=n.age.class))
  vac.per[routine.vac.age.index] <- routine.vac #xxamy
  
  rc <- new("ID.transition.SIR.vac.SIA",
            n.epi.class = tmp@n.epi.class,
            epi.class = tmp@epi.class,
            s.inds = tmp@s.inds,
            i.inds = tmp@i.inds,
            r.inds = tmp@r.inds,
            epi.class.label = c("S","I","R"),
            n.age.class = tmp@n.age.class,
            aging.rate = tmp@aging.rate,
            survival.rate = tmp@survival.rate,
            birth.rate = tmp@birth.rate,
            waifw = tmp@waifw,
            age.surv.matrix = tmp@age.surv.matrix,
            vac.per = vac.per,
            sia.vac = sia.vac,
            sia.vsucc = sia.vsucc,
            introduction.rate = introduction.rate,
            exponent = exponent,
            age.class = tmp@age.class,
            frequency.dep=frequency.dep,
            is.stochastic=is.stochastic,
            get.births = get.births
  )
  
  return(rc)
}

#Class to represent transitions in the context of vaccination
#and maternal antibodies.
setClass("ID.transition.MSIRV",
         representation(m.inds = "numeric", #indexes of those maternally protected
                        v.inds = "numeric"), #index of those who have been vaccinated
         contains = "ID.transition.SIR.vac"
)

#method to create an ID.transition.MSIRV object
create.ID.transition.MSIRV <- function(n.age.class,
                                       aging.rate,
                                       survival.rate,
                                       waifw,
                                       routine.vac = 0, #xxamy
                                       routine.vac.age.index = 12, #xxamy age index based on age classes
                                       maternal.obj, #maternal antibody object
                                       time.step = 1/2, #time step in months == generation time
                                       age.class = numeric(),
                                       birth.rate = 0,
                                       introduction.rate = 0,
                                       exponent = exponent,
                                       frequency.dep = F,
                                       is.stochastic = F,
                                       get.births = function(pop,tran){return(c(tran@birth.rate,rep(0,length(pop)-1)))}) {
  
  
  #First do all of the standard logic
  rc <- new("ID.transition.MSIRV",
            n.age.class = n.age.class,
            n.epi.class = 5,
            epi.class.label = c("M","S","I","R","V"),
            epi.class = rep(1:5, n.age.class),
            age.class = age.class,
            aging.rate = aging.rate,
            survival.rate = survival.rate,
            birth.rate = birth.rate,
            waifw = waifw,
            introduction.rate = introduction.rate,
            exponent = exponent,
            frequency.dep=frequency.dep,
            is.stochastic=is.stochastic,
            get.births = get.births)
  
  rc@m.inds = which(rc@epi.class==1)
  rc@s.inds = which(rc@epi.class==2)
  rc@i.inds = which(rc@epi.class==3)
  rc@r.inds = which(rc@epi.class==4)
  rc@v.inds = which(rc@epi.class==5)
  
  rc@age.surv.matrix <- matrix(0,nrow = n.age.class*rc@n.epi.class,
                               ncol = n.age.class*rc@n.epi.class)
  
  
  #for zeroing out all the appropriate transitions
  template.mtrx <- matrix(c(1,1,0,0,1,0,1,1,0,1,0,0,1,1,0,0,0,0,1,0,0,0,0,0,1), nrow=5, ncol=5)
  
  #useful low age
  low.age <-c(0, age.class[2:length(age.class)-1])
  
  #for each age class
  for (i in 1:n.age.class) {
    #first fill in the diagonal matrix
    tmp.template <- template.mtrx
    tmp.template[2,1] <- 0
    inds <- (i-1)*rc@n.epi.class+(1:rc@n.epi.class)
    rc@age.surv.matrix[inds,inds] <- (1-aging.rate[i])*survival.rate[i] *
      tmp.template
    
    #now fill in the off diagnal matrix
    if (i!=n.age.class) {
      #create a modified template matrix to move
      #people from the M class to the S class
      tmp.template <- template.mtrx
      
      #calculate the probability of losing maternal protection during
      #an age class
      p.mat.loss <- (pmaternal(low.age[i], maternal.obj)-
                       pmaternal(age.class[i], maternal.obj))/pmaternal(low.age[i],maternal.obj)
      
      #debug - Justin?
      p.mat.loss[pmaternal(low.age[i],maternal.obj)==0] <- 1
      #print(c(age.class[i],p.mat.loss))
      
      tmp.template[1,1] <- 1-p.mat.loss
      tmp.template[2,1] <- 1-tmp.template[1,1]
      
      #cat("age.class.i",age.class[i],"\n")
      #print(tmp.template)
      
      
      inds2 <- (i)*rc@n.epi.class+(1:rc@n.epi.class)
      rc@age.surv.matrix[inds2,inds] <- aging.rate[i]*survival.rate[i] *
        tmp.template
    }
  }
  
  #Next do all of the vaccination logic
  
  #get the %vaccinated in each age group per time step time step must be in years!
  #vac.per <- create.vacc.per.time.step.by.age(vpdf, age.class, time.step) #xxamy
  vac.per <- new("vacc.per.time.step.by.age", #xxamy
                 pvacc.in.age.class = numeric(length=n.age.class))
  vac.per@pvacc.in.age.class[routine.vac.age.index] <- routine.vac #xxamy
  rc@vac.per <- vac.per
  
  #print(range(vac.per@pvacc.in.age.class))
  
  return(rc)
}


#Class to represent transitions in the context of vaccination
#and maternal antibodies.
setClass("ID.transition.MSIRV.SIA",
         representation(m.inds = "numeric", #indexes of those maternally protected
                        v.inds = "numeric"), #index of those who have been vaccinated
         contains = "ID.transition.SIR.vac.SIA"
)

#method to create an ID.transition.MSIRV object
create.ID.transition.MSIRV.SIA <- function(n.age.class,
                                           aging.rate,
                                           survival.rate,
                                           waifw,
                                           routine.vac = 0, #xxamy
                                           routine.vac.age.index = 12, #xxamy age index based on age classes
                                           maternal.obj, #maternal antibody object
                                           time.step = 1/2, #time step in months == generation time
                                           age.class = numeric(),
                                           birth.rate = 0,
                                           sia.vac = 0,
                                           sia.vsucc =  new("vsucc.constant", success.rate=1),
                                           introduction.rate = 0,
                                           exponent = exponent,
                                           frequency.dep = F,
                                           is.stochastic = F,
                                           get.births = function(pop,tran){return(c(tran@birth.rate,rep(0,length(pop)-1)))}) {
  
  
  #First do all of the standard logic
  rc <- new("ID.transition.MSIRV.SIA",
            n.age.class = n.age.class,
            n.epi.class = 5,
            epi.class.label = c("M","S","I","R","V"),
            epi.class = rep(1:5, n.age.class),
            age.class = age.class,
            aging.rate = aging.rate,
            survival.rate = survival.rate,
            birth.rate = birth.rate,
            waifw = waifw,
            sia.vsucc = sia.vsucc,
            introduction.rate = introduction.rate,
            exponent = exponent,
            frequency.dep=frequency.dep,
            is.stochastic=is.stochastic,
            get.births = get.births)
  
  rc@m.inds = which(rc@epi.class==1)
  rc@s.inds = which(rc@epi.class==2)
  rc@i.inds = which(rc@epi.class==3)
  rc@r.inds = which(rc@epi.class==4)
  rc@v.inds = which(rc@epi.class==5)
  
  rc@age.surv.matrix <- matrix(0,nrow = n.age.class*rc@n.epi.class,
                               ncol = n.age.class*rc@n.epi.class)
  
  
  #for zeroing out all the appropriate transitions
  template.mtrx <- matrix(c(1,1,0,0,1,0,1,1,0,1,0,0,1,1,0,0,0,0,1,0,0,0,0,0,1), nrow=5, ncol=5)
  
  #useful low age
  low.age <-c(0, age.class[2:length(age.class)-1])
  
  #for each age class
  for (i in 1:n.age.class) {
    #first fill in the diagonal matrix
    tmp.template <- template.mtrx
    tmp.template[2,1] <- 0
    inds <- (i-1)*rc@n.epi.class+(1:rc@n.epi.class)
    rc@age.surv.matrix[inds,inds] <- (1-aging.rate[i])*survival.rate[i] *
      tmp.template
    
    #now fill in the off diagnal matrix
    if (i!=n.age.class) {
      #create a modified template matrix to move
      #people from the M class to the S class
      tmp.template <- template.mtrx
      
      #calculate the probability of losing maternal protection during
      #an age class
      p.mat.loss <- (pmaternal(low.age[i], maternal.obj)-
                       pmaternal(age.class[i], maternal.obj))/pmaternal(low.age[i],maternal.obj)
      
      #debug - Justin?
      p.mat.loss[pmaternal(low.age[i],maternal.obj)==0] <- 1
      #print(c(age.class[i],p.mat.loss))
      
      tmp.template[1,1] <- 1-p.mat.loss
      tmp.template[2,1] <- 1-tmp.template[1,1]
      
      #cat("age.class.i",age.class[i],"\n")
      #print(tmp.template)
      
      
      inds2 <- (i)*rc@n.epi.class+(1:rc@n.epi.class)
      rc@age.surv.matrix[inds2,inds] <- aging.rate[i]*survival.rate[i] *
        tmp.template
    }
  }
  
  #Next do all of the vaccination logic
  
  #get the %vaccinated in each age group per time step time step must be in years!
  #vac.per <- create.vacc.per.time.step.by.age(vpdf, age.class, time.step) #xxamy
  vac.per <- new("vacc.per.time.step.by.age", #xxamy
                 pvacc.in.age.class = numeric(length=n.age.class))
  vac.per@pvacc.in.age.class[routine.vac.age.index] <- routine.vac #xxamy
  rc@vac.per <- vac.per
  rc@sia.vac <- sia.vac
  
  #print(range(vac.per@pvacc.in.age.class))
  
  return(rc)
}


#### Experimental Infrastructure (naming experimental objects) ####

#Class to define an experiment.
setClass("experiment",
         representation(name = "character", #a human readable name for the experiment
                        R0 = "numeric", #if included the R0 we scale to
                        state.t0 = "ID.state.matrix",#the initial state
                        t.min = "numeric", #the starting time for the experiment in years
                        t.max = "numeric", #the ending time for the experiment in years
                        t0.doy = "numeric", #the starting day of year, should default to 0
                        step.size = "numeric", #the step size in years, 1/no.gens.per.year
                        trans = "ID.transition.SIR", #defines the transitions for this experiment
                        season.obj = "seasonal.force", #seasonal forcing function
                        births.per.1000.each.timestep = "ANY", #the estimated crude birth rate per time step, length is T
                        description = "character" #describe the experiment
         ))


#Class defines an experiment with SIAs or opportunity for triggered SIAs
setClass("experiment.SIAtrigger",
         representation(trans = "ID.transition.SIR.vac", #the transitions for this experiment
                        sia.obj = "matrix", #defines the SIAs that will go on during this experiment
                        trigger.list="list" #defined the triggered vaccination events that occur in this experiment
         ),
         contains="experiment")

#Class for an experiment with birth and death rates updated each time step and population rescaling -xxamy
setClass("experiment.updatedemog",
         representation(trans = "ID.transition.SIR.vac", #the transitions for this experiment - xxamy added this
                        surv.each.timestep = "matrix", #the estimated age-specific survival rates over the length of the experiment, ncol is T
                        pop.rescale.each.timestep = "ANY", #a vector of NA's if not call for population rescale, if !=0 then rescale by this pop number
                        maternal.obj = "maternal.exp.decay"),  #maternal antibody object
         contains=c("experiment")) 

#Class for an experiment with changing vaccination coverage over time
setClass("experiment.updatedemog.vaccinationchange",
         representation(time.specific.MR1cov = "ANY", #vector of values of the coverage values
                        time.specific.MR2cov = "ANY", #vector of values of the coverage values
                        time.specific.SIAcov = "ANY", #vector of values of the coverage values
                        time.specific.min.age.MR1 = "ANY",
                        time.specific.max.age.MR1 = "ANY",
                        time.specific.min.age.MR2 = "ANY",
                        time.specific.max.age.MR2 = "ANY",
                        time.specific.min.age.SIA = "ANY",
                        time.specific.max.age.SIA = "ANY",
                        obj.vcdf.MR1 = "vaccine.cdf.byage",
                        obj.vcdf.MR2 = "vaccine.cdf.byage",
                        obj.prob.vsucc = "prob.vsucc.byage",
                        sia.timing.in.year = "ANY"), #time of the SIA in years (e.g., 0.5 is July of the year)
         contains=c("experiment.updatedemog"))

#Class for an experiment 
setClass("experiment.updatedemog.vaccinationchange.vaccinationlimitations",
         representation(MR1MR2correlation = "logical", 
                        MR1SIAcorrelation = "logical", 
                        MR2SIAcorrelation = "logical", 
                        SIAinefficient = "logical",
                        SIAinacc = "logical",
                        prop.inacc = "ANY"),
         contains=c("experiment.updatedemog.vaccinationchange"))

#### Running the Experiment (run methods) ####

setGeneric("run", function(exper,...) standardGeneric("run"))

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
                  sia.vacc.prob[index.sia.vacc[t],] - 
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
                    sia.vacc.prob[index.sia.vacc[t],] -
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


# Last and final transition function to be called
SIRIDTran <- function (state, tran) {
  
  #define denominator of phi matrix
  if (tran@frequency.dep) {denom<-sum(state)} else {denom<-1}
  
  ## ## #define the phi matrix
  ## phi <- (tran@waifw%*%(state[tran@i.inds,]^tran@exponent))/denom
  
  ## phi <- 1 - exp(-phi)
  
  ## phi <- matrix(phi, nrow=state@n.age.class ,
  ##               ncol=state@n.age.class)
  
  ## phi <- t(phi)
  
  
  ## #get the transition matrix
  ## tran.matrix <- tran@age.surv.matrix
  
  ## #make susceptible part of matrix
  ## tran.matrix[tran@s.inds, tran@s.inds] <-
  ##     tran.matrix[tran@s.inds, tran@s.inds] * (1-phi)
  
  ## #make susceptible part of matrix
  ## tran.matrix[tran@i.inds, tran@s.inds] <-
  ##     tran.matrix[tran@i.inds, tran@s.inds] * (phi)
  
  ## #no one stays infecrted...might get more sophisticated later
  ## tran.matrix[tran@i.inds, tran@i.inds] <- 0
  
  tran.matrix <- .Call("calc_phi_and_update_tran",
                       tran@waifw,
                       state[,1],
                       tran@s.inds,
                       tran@i.inds,
                       tran@exponent,
                       denom,
                       tran@age.surv.matrix)
  
  
  if (!tran@is.stochastic) {
    
    birthst <- tran@get.births(state[,1],tran)
    #print(tran@get.births)
    #print(tran@birth.rate)
    #print(birthst[1])#xxjnow
    
    state[,] <- tran.matrix%*%state
    
    #add in the births to 1,1 for now, assuming that is correct.
    #might go back on this later.
    #state[1,1] <-  state[1,1] + tran@birth.rate #A BIT OF A HACK
    state[,1] <-  state[,1] + birthst #xxj - current births based on function
    #crude way of handling introductions
    state[tran@i.inds,] <- state[tran@i.inds,] +
      tran@introduction.rate
    
    #print(range(state[tran@i.inds,]))
    
  } else {
    birthst <- tran@get.births(state[,1],tran)
    #mortality probablity in each category of the n.age.class * no classes
    mort <- 1-rep(tran@survival.rate, each=state@n.epi.class)
    
    if (is.loaded("do_ID_transition_SIR_stochastic_moves")) {
      #C implementation of this.
      #state[,1] <- .C("do_ID_transition_SIR_stochastic_moves",
      #                as.integer(state[,1]),
      #                as.integer(length(state[,1])),
      #                as.double(tran.matrix),
      #                newstate=as.integer(state[,1]))$newstate
      
      state[,1] <- .Call("do_ID_transition_SIR_stochastic_moves_cl",
                         as.integer(state[,1]),
                         tran.matrix)
    } else {
      #loop over and distribute via a multinomial
      newstate <- rep(0,1+length(tran.matrix[,1]))
      for (k in 1:length(tran.matrix[1,])) {
        newstate <- newstate+
          rmultinom(1,state[k,1],c(tran.matrix[,k],mort[k]))
      }
      
      state[,1] <- newstate[1:length(tran.matrix[,1])]
      
    }
    
    
    #add in the births to 1,1 for now, assuming that is correct.
    #might go back on this later.
    #state[1,1] <-  state[1,1] + rpois(1,tran@birth.rate) #A BIT OF A HACK
    state[,1] <-  state[,1] +rpois(length(birthst),birthst)#xxj - current births based on function
    
    
    #crude way of handling introductions
    state[tran@i.inds,1] <- state[tran@i.inds,1] +
      rpois(length(tran@i.inds),tran@introduction.rate)
    
    
  }
  
  
  
  
  return(state)
}


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

#### Experiment Result Objects (naming experimental results objects) ####

# Class to hold SIR simulation results. Hold the results of the simulation as
# SIR compartment values at each time step.
setClass("sim.results.SIR",
         representation(s.inds = "numeric", #indexing of susceptibles
                        i.inds = "numeric", #infectious indexes
                        r.inds = "numeric", #recovered indexes
                        age.class = "numeric", #upper end of age classes
                        t = "numeric"#times when we have info in the matrix
         ),
         contains="matrix")

# Holds the results of a simulation with an MSIRV object
setClass("sim.results.MSIRV",
         representation(m.inds = "numeric",
                        v.inds = "numeric",
                        routine.intro = "numeric",
                        sia.times = "numeric"),
         contains= "sim.results.SIR")

# Holds the results of a simulation with an experiment.updatedemog object
setClass("sim.results.MSIRV.update.demog",
         representation(births.each.timestep = "ANY", #the number of births per time step as output from simulation
                        growth.rate.each.timestep = "ANY"
         ),
         contains="sim.results.MSIRV")

# Hols the results output for a simulation with an experiment.updatedemog.vaccinationchange object
setClass("sim.results.MSIRV.update.demog.vaccine.change",
         representation(MR1.fail.each.timestep = "ANY", #the number of births per time step as output from simulation
                        MR2.fail.each.timestep = "ANY",
                        SIA.fail.each.timestep = "ANY"
         ),
         contains="sim.results.MSIRV.update.demog")

# Holds the results of an experiment along with the experiment definition
setClass("experiment.result",
         representation(experiment.def = "experiment",
                        result = "sim.results.SIR"))



#### Methods to Manipulate or Plot Experiment Output ####

#Plot a sim.results.SIR object
#does a 4 panel plot that shows M+S,I,R+V and average age
#
#Parameters -
#     x - the sim.results.MSIRV object
#     y - ignored
#
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


