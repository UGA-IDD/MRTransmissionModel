#### Experiment Result Objects (naming experimental results objects) ###

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
