#### Experimental Infrastructure (naming experimental objects) ####


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

