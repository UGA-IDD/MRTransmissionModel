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
