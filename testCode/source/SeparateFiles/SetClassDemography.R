#### Demography ####

#Class definition for transition matrix. This matrix knows how to shift the
#ID.state.matrix from one time step to another, and contains all the fun stuff
#we need to do the transition (i.e., the matrix)
setClass("nMx",
         slots = list(rates="data.frame", # age specific death rates per 1 from 1950 to 2100 by 5 year increments
                      mid.age="numeric" # mid-age associated with rates,
         ))
