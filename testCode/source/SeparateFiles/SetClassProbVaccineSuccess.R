## object for the prob of vaccine success by age
setClass("prob.vsucc.byage",
         slots = list(prob.vsucc = "numeric", #vector of length ages
                        ages = "numeric" #vector of length prob
         )
)
