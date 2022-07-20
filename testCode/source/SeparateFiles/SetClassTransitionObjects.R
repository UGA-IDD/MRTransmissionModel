#### Transition Objects (object that holds everything needed to make transitions) ####

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




