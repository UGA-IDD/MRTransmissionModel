#' Object for the probability of vaccine success by age
#'
#' @slot prob.vsucc numeric.
#' @slot ages numeric.
#'
#' @export
#' @docType class
#' @rdname prob.vsucc.byage-class
setClass("prob.vsucc.byage",
         slots = list(prob.vsucc = "numeric", #vector of length ages
                      ages = "numeric" #vector of length prob
         )
)


#' Class for maternal immunity
#'
#' Assumes that protection is the same as exponential decay
#'
#' @slot decay.rt numeric.
#'
#' @export
#' @docType class
#' @rdname maternal.exp.decay-class
setClass("maternal.exp.decay",
         slots = list(decay.rt = "numeric")#the exponential decay rate
)


#' Class for maternal immunity threshold
#'
#' Assume that there is some threshold age at which
#' you are 100% protected if this age or younger
#'
#' @slot thresh numeric.
#'
#' @export
#' @docType class
#' @rdname maternal.thresh-class
setClass("maternal.thresh",
         slots = list(thresh="numeric"))



#' Class for seasonal forcing
#'
#' @export
#' @docType class
#' @rdname seasonal.force-class
setClass("seasonal.force")



#' Class for seasonal cosine
#'
#' @slot amplitude numeric.
#' @slot offset numeric.
#' @slot period numeric.
#'
#' @export
#' @docType class
#' @rdname seasonal.cosine-class
setClass("seasonal.cosine",
         slots = list(amplitude = "numeric", #amplitude of seasonal forcing
                      offset = "numeric", #offset from t=0 for nadir
                      period = "numeric"), #number of peaks per year
         prototype = list(amplitude=.2, offset = 0, period = 1),
         contains="seasonal.force"
)

#' Class for seasonal periodic bspline
#'
#' @slot parameters numeric.
#' @slot nbasis numeric.
#' @slot degree numeric.
#' @slot period numeric.
#'
#' @export
#' @docType class
#' @rdname seasonal.periodic.bspline-class
setClass("seasonal.periodic.bspline",
         slots = list(parameters = "numeric",
                      nbasis = "numeric",
                      degree = "numeric",
                      period = "numeric"
         ),
         prototype = list(parameters= c(0.5,1,1), nbasis=3, degree=3, period=1),
         contains="seasonal.force",
)

#' Class for age specific seasonal trend
#'
#' @slot amplitude numeric.
#' @slot offset numeric.
#' @slot period numeric.
#' @slot age.class numeric.
#' @slot chosen.age numeric.
#'
#' @export
#' @docType class
#' @rdname seasonal.age.specific-class
setClass("seasonal.age.specific",
         slots = list(amplitude = "numeric", #amplitude of seasonal forcing
                      offset = "numeric", #offset from t=0 for nadir
                      period = "numeric", #number of peaks per year
                      age.class = "numeric",#age classes in the WAIFW
                      chosen.age = "numeric"),#vector of length 2, indicating start and end of seasonally forced age class
         prototype = list(amplitude=.2, offset = 0, period = 1, age.class=c(1:60,60*12),chosen.age=c(50,60)),
         contains="seasonal.force",
)


#' Seasonal class
#'
#' Designed to take output from TSIR-type analysis
#'
#' @slot parameters numeric.
#' @slot time.in.year numeric.
#'
#' @export
#' @docType class
#' @rdname seasonal.piecewise.constant-class
setClass("seasonal.piecewise.constant",
         slots = list(parameters = "numeric",
                      time.in.year = "numeric"),
         prototype = list(parameters = rep(1,24),time.in.year = (0:24)/24),
         contains="seasonal.force",
)


#' Seasonal class
#'
#' Designed to take output from TSIR-type analysis
#'
#' @slot parameters numeric.
#' @slot time.in.year numeric.
#' @slot age.class numeric.
#' @slot chosen.age numeric.
#'
#' @export
#' @docType class
#' @rdname seasonal.age.specific.piecewise.constant-class
setClass("seasonal.age.specific.piecewise.constant",
         slots = list(parameters = "numeric",
                      time.in.year = "numeric",
                      age.class = "numeric",#age classes in the WAIFW
                      chosen.age = "numeric"),#vector of length 2, indicating start and end of seasonally forced age class
         prototype = list(parameters = rep(1,24),time.in.year = (0:24)/24),
         contains="seasonal.force",
)


#' Generate object of class ID.state.matrix
#'
#' @slot n.epi.class numeric.
#' @slot epi.class numeric.
#' @slot epi.class.label character.
#' @slot n.age.class numeric.
#'
#' @docType class
#' @rdname ID.state.matrix-class

setClass("ID.state.matrix",
         slots = list(n.epi.class = "numeric", #number of different
                      #epi classes
                      epi.class = "numeric", #the epid class of each row,
                      #between 1 and n.epi.class
                      epi.class.label = "character", #label for each epi class
                      n.age.class = "numeric"#the number of age classes
         ),
         contains="matrix")

#' Generate object of class vacc.per.time.step.by.age
#'
#' Class that holds the information on how many people should be
#' vaccinated during each time step in each age class
#' in the TSIR framework.
#'
#' @slot pvacc.in.age.class numeric.
#'
#' @export
#' @docType class
#' @rdname vacc.per.time.step.by.age-class

setClass("vacc.per.time.step.by.age",
         slots = list(pvacc.in.age.class = "numeric")
         # vector of how many are vaccinated in the age class per time step
)

#' Generates an object of class ID.transition.SIR
#'
#' Class definition for transition matrix. This matrix knows how to shift the
#' ID.state.matrix from one time step to another, and contains all the fun stuff
#' we need to do the transition (i.e., the matrix)
#'
#'
#' @slot n.epi.class numeric.
#' @slot epi.class numeric.
#' @slot epi.class.label character.
#' @slot s.inds numeric.
#' @slot i.inds numeric.
#' @slot r.inds numeric.
#' @slot n.age.class numeric.
#' @slot age.class numeric.
#' @slot aging.rate numeric.
#' @slot survival.rate numeric.
#' @slot birth.rate numeric.
#' @slot waifw matrix.
#' @slot age.surv.matrix matrix.
#' @slot introduction.rate numeric.
#' @slot exponent numeric.
#' @slot frequency.dep logical.
#' @slot is.stochastic logical.
#' @slot get.births function.
#'
#' @export
#'
#' @docType class
#' @rdname ID.transition.SIR-class
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
                        survival.rate = "ANY", #the percent who survive in this age class at each time step -- xxamy changed to "ANY"
                        birth.rate = "ANY", #the birth rate for this transition -- xxamy changed to "ANY"
                        waifw = "matrix", #who acquires infection from whom matrix
                        age.surv.matrix = "matrix", #matrix governing age/survival
                        introduction.rate = "numeric", #an vector of possibly 0 case importation rates in each age class
                        exponent = "numeric",        #exponent on infecteds to correct for discretization
                        frequency.dep = "logical",	#true/false indicating if freq. dependence or density dep. is implemented
                        is.stochastic = "logical",   #true/false indicating if stochasticity desired
                        get.births = "function" #a function that takes state and birth rate; allowing births to be shaped by pop struct
         ))

#' Generates object of class ID.transition.SIR.vac
#'
#' Class to represent transitions in the context of vaccination
#'
#' @slot vac.per vacc.per.time.step.by.age.
#'
#' @export
#'
#' @docType class
#' @rdname ID.transition.SIR.vac-class

setClass("ID.transition.SIR.vac",
         slots = list(vac.per="vacc.per.time.step.by.age"), #percent routinely vaccinated in each age class
         contains="ID.transition.SIR")


#' Make the generic vaccine success object
#'
#' @export
#'
#' @docType class
#' @rdname vsucc-class
setClass("vsucc")

#' Make a constants vaccine success object
#'
#' @slot success.rate numeric.
#'
#' @export
#' @docType class
#' @rdname vsucc.constant-class
setClass("vsucc.constant",
         slots = list(success.rate="numeric"),
         contains="vsucc")


#' Make an object to hold a logistic function for vaccine success
#'
#' @slot intercept numeric
#' @slot mo.eff numeric
#' @slot full.efficacy numeric
#'
#' @export
#' @docType class
#' @rdname vsucc.Logistic-class
setClass("vsucc.Logistic",
         slots = list(intercept = "numeric",
                      mo.eff = "numeric",
                      full.efficacy = "numeric" #the max reachable efficacy
         ),
         contains = "vsucc")



#' Generates object of class ID.transition.SIR.vac.SIA
#'
#' @slot sia.vac numeric.
#' @slot sia.vsucc vsucc.
#'
#' @export
#'
#' @docType class
#' @rdname ID.transition.SIR.vac.SIA-class
setClass("ID.transition.SIR.vac.SIA",
         slots = list(sia.vac = "numeric",  #number vaccinated in each age group by SIA on this transition
                      sia.vsucc = "vsucc" #object to determine the probability that an SIA vaccination is successful
         ),
         contains="ID.transition.SIR.vac")


#' Generates an object of class ID.transitions.MSIRV.SIA
#'
#' Class to represent transitions in the context of vaccination
#' and maternal antibodies.
#'
#' @slot m.inds numeric
#' @slot v.inds numeric
#'
#' @export
#' @docType class
#' @rdname ID.transisition.MSIRV.SIA-class
setClass("ID.transition.MSIRV.SIA",
         slots = list(m.inds = "numeric", #indexes of those maternally protected
                      v.inds = "numeric"), #index of those who have been vaccinated
         contains = "ID.transition.SIR.vac.SIA"
)

#' Class to represent transitions in the context of vaccination and maternal antibodies with multiple locations
#'
#' @slot n.subpops numeric
#' @slot subpop.class.label numeric
#' @slot coupling matrix
#'
#' @export
#' @docType class
#' @rdname ID.transition.MSIRV.space-class
setClass("ID.transition.MSIRV.space",
         representation(n.subpops = "numeric", #total number of locations
                        subpop.class.label = "numeric", #location index
                        coupling = "matrix"), #matrix of connection between each sub-population
         contains = "ID.transition.MSIRV.SIA"
)


#' Generate object of class ID.transition.MSIRV
#'
#' Class to represent transitions in the context of vaccination
#' and maternal antibodies.
#'
#' @slot m.inds numeric.
#' @slot v.inds numeric.
#'
#' @export
#'
#' @docType class
#' @rdname ID.transition.MSIRV-class
#'

setClass("ID.transition.MSIRV",
         slots = list(m.inds = "numeric", #indexes of those maternally protected
                      v.inds = "numeric"), #index of those who have been vaccinated
         contains = "ID.transition.SIR.vac"
)

#' Generate object of class nMx
#'
#' Class definition for transition matrix. This matrix knows how to shift the
#' ID.state.matrix from one time step to another, and contains all the fun stuff
#' we need to do the transition (i.e., the matrix)
#'
#' @slot rate.years numeric.
#' @slot rates data.frame.
#' @slot mid.age numeric.
#'
#' @export
#' @docType class
#' @rdname nMx-class
setClass("nMx",
         slots = list(rate.years="numeric", #year associated with each time-specific rates --xxamy added this new object rate.years
                        rates="data.frame", #age specific death rates per 1 (rows) by time (columns)
                        mid.age="numeric" #mid-age associated with age-specific rates
         ))

#' Generate and object of class vaccine.cdf.byage
#'
#' Object for the CDF of vaccination by age (non scaled)
#'
#' @slot cdf numeric.
#' @slot ages numeric.
#'
#' @export
#' @docType class
#' @rdname vaccine.cdf.byage-class
setClass("vaccine.cdf.byage",
         slots = list(cdf = "numeric", #vector of length ages
                      ages = "numeric" #vector of length cdf
         )
)


#' Generate object of class cov.estimates
#'
#' @slot sia.cov numeric.
#' @slot sia.min.age numeric.
#' @slot sia.max.age numeric.
#' @slot MR1.cov numeric.
#' @slot MR2.cov numeric.
#' @slot inaccessible.prop numeric.
#'
#' @export
#' @docType class
#' @rdname cov.estimates-class
setClass("cov.estimates",
         slots = list(sia.cov = "numeric",
                      sia.min.age = "numeric",
                      sia.max.age = "numeric",
                      MR1.cov = "numeric",
                      MR2.cov = "numeric",
                      inaccessible.prop = "numeric"
         ))


#' Class to define an experiment.
#'
#' @slot name character.
#' @slot R0 numeric.
#' @slot state.t0 ID.state.matrix.
#' @slot t.min numeric.
#' @slot t.max numeric.
#' @slot t0.doy numeric.
#' @slot step.size numeric.
#' @slot trans ID.transition.SIR.
#' @slot season.obj seasonal.force.
#' @slot births.per.1000.each.timestep numeric.
#' @slot description character.
#'
#' @export
#' @docType class
#' @rdname experiment-class
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
                      births.per.1000.each.timestep = "ANY", #the estimated crude birth rate per time step, length is T
                      description = "character" #describe the experiment
         ))


#' Class defines an experiment with SIAs or opportunity for triggered SIAs
#' @slot trans ID.transition.SIR.vac.
#' @slot sia.obj matrix.
#' @slot trigger.list list.
#'
#' @export
#' @docType class
#' @rdname experiment.SIAtrigger-class
setClass("experiment.SIAtrigger",
         slots = list(trans = "ID.transition.SIR.vac", #the transitions for this experiment
                      sia.obj = "matrix", #defines the SIAs that will go on during this experiment
                      trigger.list="list" #defined the triggered vaccination events that occur in this experiment
         ),
         contains="experiment")


#' Class for an experiment with birth and death rates updated each time step and population rescaling
#'
#' @slot trans ID.transition.SIR.vac.
#' @slot surv.each.timestep matrix.
#' @slot pop.rescale.each.timestep numeric.
#' @slot maternal.obj maternal.exp.decay.
#'
#' @export
#' @docType class
#' @rdname experiment.updatedemog-class
setClass("experiment.updatedemog",
         slots = list(trans = "ID.transition.SIR.vac", #the transitions for this experiment - xxamy added this
                      surv.each.timestep = "matrix", #the estimated age-specific survival rates over the length of the experiment, ncol is T
                      pop.rescale.each.timestep = "ANY", #a vector of 0's if not call for population rescale, if !=0 then rescale by this pop number
                      maternal.obj = "maternal.exp.decay"),  #maternal antibody object
         contains=c("experiment"))



#' Class for a spatial experiment updating demography over time
#'
#' @slot trans ID.transition.MSIRV.space.
#'
#' @export
#' @docType class
#' @rdname experiment.updatedemog.spatial-class
setClass("experiment.updatedemog.spatial",
         slots = list(trans = "ID.transition.MSIRV.space"),
         contains=c("experiment.updatedemog"))



#' Class for an experiment with changing vaccination coverage over time
#'
#' @slot time.specific.MR1cov vector.
#' @slot time.specific.MR2cov vector.
#' @slot time.specific.SIAcov vector.
#' @slot time.specific.min.age.MR1 vector.
#' @slot time.specific.max.age.MR1 vector.
#' @slot time.specific.min.age.MR2 vector.
#' @slot time.specific.max.age.MR2 vector.
#' @slot time.specific.min.age.SIA vector.
#' @slot time.specific.max.age.SIA vector.
#' @slot obj.vcdf.MR1 vaccine.cdf.byage.
#' @slot obj.vcdf.MR2 vaccine.cdf.byage.
#' @slot obj.prob.vsucc prob.vsucc.byage.
#' @slot sia.timing.in.year vector.
#'
#' @export
#' @docType class
#' @rdname experiment.updatedemog.vaccinationchange-class
setClass("experiment.updatedemog.vaccinationchange",
         slots = list(time.specific.MR1cov = "ANY", #vector of values of the coverage values --xxamy changed to "ANY"
                        time.specific.MR2cov = "ANY", #vector of values of the coverage values --xxamy changed to "ANY"
                        time.specific.SIAcov = "ANY", #vector of values of the coverage values --xxamy changed to "ANY"
                        time.specific.min.age.MR1 = "ANY", #--xxamy changed to "ANY"
                        time.specific.max.age.MR1 = "ANY", #--xxamy changed to "ANY"
                        time.specific.min.age.MR2 = "ANY", #--xxamy changed to "ANY"
                        time.specific.max.age.MR2 = "ANY", #--xxamy changed to "ANY"
                        time.specific.min.age.SIA = "ANY", #--xxamy changed to "ANY"
                        time.specific.max.age.SIA = "ANY", #--xxamy changed to "ANY"
                        obj.vcdf.MR1 = "vaccine.cdf.byage",
                        obj.vcdf.MR2 = "vaccine.cdf.byage",
                        obj.prob.vsucc = "prob.vsucc.byage",
                        sia.timing.in.year = "ANY"), #time of the SIA in years (e.g., 0.5 is July of the year) --xxamy changed to "ANY"
         contains=c("experiment.updatedemog"))


#' Class for a spatial experiment updating demography and with changing vaccination coverage over time
#'
#' @slot trans ID.transition.MSIRV.space.
#'
#' @export
#' @docType class
#' @rdname experiment.updatedemog.vaccinationchange.spatial-class
setClass("experiment.updatedemog.vaccinationchange.spatial",
         slots = list(trans = "ID.transition.MSIRV.space"),
         contains=c("experiment.updatedemog.vaccinationchange"))


#' Class for an experiment
#'
#' @slot MR1MR2correlation logical.
#' @slot MR1SIAcorrelation logical.
#' @slot MR2SIAcorrelation logical.
#' @slot SIAinefficient logical.
#' @slot SIAinacc logical.
#' @slot prop.inacc numeric.
#'
#' @export
#' @docType class
#' @rdname experiment.updatedemog.vaccinationchange.vaccinationlimitations-class
setClass("experiment.updatedemog.vaccinationchange.vaccinationlimitations",
         slots = list(MR1MR2correlation = "logical",
                      MR1SIAcorrelation = "logical",
                      MR2SIAcorrelation = "logical",
                      SIAinefficient = "logical",
                      SIAinacc = "logical",
                      prop.inacc = "ANY"),
         contains=c("experiment.updatedemog.vaccinationchange"))


#' Class for a spatial vaccination experiment updating demography and with changing vaccination coverage over time, with vaccination limitations
#'
#' @slot trans ID.transition.MSIRV.space.
#'
#' @export
#'
#' @docType class
#' @rdname experiment.updatedemog.vaccinationchange.vaccinationlimitations.spatial-class
setClass("experiment.updatedemog.vaccinationchange.vaccinationlimitations.spatial",
         representation(trans = "ID.transition.MSIRV.space"),
         contains=c("experiment.updatedemog.vaccinationchange.vaccinationlimitations"))



#' Generate of object of class experiment result (sim.result)
#'
#' Class to hold SIR simulation results. Hold the results of the simulation as
#' SIR compartment values at each time step.
#'
#' @slot .Data matrix.
#' @slot s.inds numeric.
#' @slot i.inds numeric.
#' @slot r.inds numeric.
#' @slot age.class numeric.
#' @slot t numeric.
#'
#' @export
#' @docType class
#' @rdname sim.results.SIR-class
setClass("sim.results.SIR",
         slots = list(.Data = "matrix",
                      s.inds = "numeric", #indexing of susceptibles
                      i.inds = "numeric", #infectious indexes
                      r.inds = "numeric", #recovered indexes
                      age.class = "numeric", #upper end of age classes
                      t = "numeric"#times when we have info in the matrix
         ),
         contains="matrix")

#' Holds the results of a simulation with an MSIRV object
#'
#' @slot m.inds numeric.
#' @slot v.inds numeric.
#' @slot routine.intro numeric.
#' @slot sia.times numeric.
#'
#' @export
#' @docType class
#' @rdname sim.results.MSIRV-class
setClass("sim.results.MSIRV",
         slots = list(m.inds = "numeric",
                      v.inds = "numeric",
                      routine.intro = "numeric",
                      sia.times = "numeric"),
         contains= "sim.results.SIR")


#' Holds the results of a simulation with a MSIRV.space object
#'
#' @slot n.subpops numeric.
#' @slot subpop.class.label numeric.
#'
#' @export
#' @docType class
#' @rdname sim.results.MSIRV.space-class
setClass("sim.results.MSIRV.space",
         representation(n.subpops = "numeric",
                        subpop.class.label = "numeric"),
         contains= "sim.results.MSIRV")


#' Holds the results of a simulation with an experiment.updatedemog object
#'
#' @slot births.each.timestep ANY.
#' @slot growth.rate.each.timestep ANY.
#'
#' @export
#' @docType class
#' @rdname sim.results.MSIRV.update.demog.space-class
setClass("sim.results.MSIRV.update.demog.space",
         representation(births.each.timestep = "ANY", #the number of births per time step as output from simulation
                        growth.rate.each.timestep = "ANY"
         ),
         contains="sim.results.MSIRV.space")


#' Holds the results output for a simulation with an experiment.updatedemog.vaccinationchange object
#'
#' @slot MR1.fail.each.timestep ANY.
#' @slot MR2.fail.each.timestep ANY.
#' @slot SIA.fail.each.timestep ANY.
#'
#' @export
#' @docType class
#' @rdname sim.results.MSIRV.update.demog.vaccine.change.space-class
setClass("sim.results.MSIRV.update.demog.vaccine.change.space",
         representation(MR1.fail.each.timestep = "ANY", #the number of births per time step as output from simulation
                        MR2.fail.each.timestep = "ANY",
                        SIA.fail.each.timestep = "ANY"
         ),
         contains="sim.results.MSIRV.update.demog.space")



#' Holds the results of a simulation with an experiment.updatedemog object
#'
#' @slot births.each.timestep numeric.
#' @slot growth.rate.each.timestep numeric.
#'
#' @export
#' @docType class
#' @rdname sim.results.MSIRV.update.demog-class
setClass("sim.results.MSIRV.update.demog",
         slots = list(births.each.timestep = "ANY", #the number of births per time step as output from simulation
                      growth.rate.each.timestep = "ANY"
         ),
         contains="sim.results.MSIRV")

#' Holds the results output for a simulation with an experiment.updatedemog.vaccinationchange object
#'
#' @slot MR1.fail.each.timestep numeric.
#' @slot MR2.fail.each.timestep numeric.
#' @slot SIA.fail.each.timestep numeric.
#'
#' @export
#' @docType class
#' @rdname sim.results.MSIRV.update.demog.vaccine.change-class
setClass("sim.results.MSIRV.update.demog.vaccine.change",
         slots = list(MR1.fail.each.timestep = "ANY", #the number of births per time step as output from simulation
                      MR2.fail.each.timestep = "ANY",
                      SIA.fail.each.timestep = "ANY"
         ),
         contains="sim.results.MSIRV.update.demog")

#' Holds the results of an experiment along with the experiment definition
#'
#' @slot experiment.def experiment.
#' @slot result sim.results.SIR.
#'
#' @export
#' @docType class
#' @rdname experiment.result-class
setClass("experiment.result",
         slots = list(experiment.def = "experiment",
                      result = "sim.results.SIR"))



#' Class definition for space.nMx object - aka Age Specific Death Rates (ASDR) - by space
#'
#' @slot rate.years numeric; year associated with each time-specific rate
#' @slot rates data.frame; age and space specific death rates per 1 (rows) by year (columns)
#' @slot mid.age numeric; mid-age associated with each age-specific rates
#' @slot n.subpops numeric; number of subpopulations
#'
#' @export
#' @docType class
#' @rdname space.nMx-class
setClass("space.nMx",
         representation(rate.years="numeric", #year associated with each time-specific rate
                        rates="data.frame", #age and space specific death rates per 1 (rows) by year (columns)
                        mid.age="numeric", #mid-age associated with each age-specific rates
                        n.subpops="numeric" #number of subpopulations
         ))

