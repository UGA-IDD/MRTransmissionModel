# S4 classes
# look at an existing S4 package
install.packages("distrMod")
library("distrMod")
library("methods")

?promptClass

# document the S4 Classes
promptClass("vacc.per.time.step.by.age")
promptClass("vaccine.cdf.byage")
promptClass("vsucc")
promptClass("vsucc.constant")
promptClass("vsucc.Logistic")
promptClass("ID.transition.SIR")
promptClass("ID.transition.SIR.vac")
promptClass("ID.transition.SIR.vac.SIA")
promptClass("ID.transition.MSIRV")
promptClass("ID.transition.MSIRV.SIA")
promptClass("cov.estimates")
promptClass("nMx")
promptClass("experiment")
promptClass("experiment.SIAtrigger")
promptClass("experiment.updatedemog")
promptClass("experiment.updatedemog.vaccinationchange")
promptClass("experiment.updatedemog.vaccinationchange.vaccinationlimitations")
promptClass("sim.results.SIR")
promptClass("sim.results.MSIRV")
promptClass("sim.results.MSIRV.update.demog")
promptClass("sim.results.MSIRV.update.demog.vaccine.change")
promptClass("experiment.result")
promptClass("maternal.exp.decay")
promptClass("maternal.thresh")
promptClass("prob.vsucc.byage")
promptClass("seasonal.force")
promptClass("seasonal.cosine")
promptClass("seasonal.periodic.bspline")
promptClass("seasonal.age.specific")
promptClass("seasonal.piecewise.constant")
promptClass("seasonal.age.specific.piecewise.constant")
promptClass("ID.state.matrix")

# class documentation already includes the methods for each class


##################
### START HERE ###
##################

# figure out how to document S4 methods

# document the S4 Methods
promptMethods("plot")
promptMethods()


args(getGeneric("plot"))





























