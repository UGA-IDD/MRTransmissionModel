#' Methods to get the percent immune by maternal antibodies.
#'
#'
#' @param age numeric; the age we want the decay rate for
#' @param mat.obj the object describing maternally acquired immunity
#'
#' @include setClasses.R
#' @return percent immune by maternal antibodies
#'
#' @export
#' @docType methods
#' @rdname pmaternal-methods


setGeneric("pmaternal",
           function(age, mat.obj) standardGeneric("pmaternal"))


#' @rdname pmaternal-methods
#' @aliases pmaternal,numeric,maternal.exp.decay-method
setMethod("pmaternal",
          c("numeric", "maternal.exp.decay"),
          function (age, mat.obj) {
            return(exp(-age*mat.obj@decay.rt))
          })

#' @rdname pmaternal-methods
#' @aliases pmaternal,numeric,maternal.thresh-method
setMethod("pmaternal",
          c("numeric", "maternal.thresh"),
          function (age, mat.obj) {
            rc <- sapply(age, function(age) {
              if(age<=mat.obj@thresh) {return(1)}
              return(0)
            })
            return(rc)
          })
