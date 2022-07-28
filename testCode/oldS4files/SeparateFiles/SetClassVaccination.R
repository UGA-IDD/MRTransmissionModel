#### Vaccination ####

#Class that holds the information on how many people should be
#vaccinated during each time step in each age class
#in the TSIR framework.
setClass("vacc.per.time.step.by.age",
         slots = list(pvacc.in.age.class = "numeric") #vector of how many are
         #vaccinated in the age
         #class per time step
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
