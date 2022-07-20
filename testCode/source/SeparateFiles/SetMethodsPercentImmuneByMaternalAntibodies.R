# Methods to get the percent immune by maternal antibodies.
#
#Parameters -
#    age - the age we want the decay rate for
#    mat.obj - the object describing maternally acquired immunity
setGeneric("pmaternal",
           function(age, mat.obj) standardGeneric("pmaternal"))

setMethod("pmaternal",
          c("numeric", "maternal.exp.decay"),
          function (age, mat.obj) {
            return(exp(-age*mat.obj@decay.rt))
          })


setMethod("pmaternal",
          c("numeric", "maternal.thresh"),
          function (age, mat.obj) {
            rc <- sapply(age, function(age) {
              if(age<=mat.obj@thresh) {return(1)}
              return(0)
            })
            return(rc)
          })
