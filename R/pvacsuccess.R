#' Methods for probability of successful vaccination
#'
#' @param age the age in months we are checking
#' @param vsucc.obj the vaccine success object
#'
#' @importFrom stats plogis
#'
#' @include setClasses.R
#' @return the probability of successful vaccination at that age
#'
#' @export
#' @docType methods
#' @rdname pvacsuccess-methods

setGeneric("pvacsuccess",
           function(age, vsucc.obj) standardGeneric("pvacsuccess"))

#' @rdname pvacsuccess-methods
#' @aliases pvacsuccess,numeric,vsucc.constant-method
setMethod("pvacsuccess",
          c("numeric","vsucc.constant"),
          function (age, vsucc.obj) {

            prob.vsucc <- new("prob.vsucc.byage",
                              prob.vsucc=rep(vsucc.obj@success.rate, length(age)),
                              ages=age)

            return(prob.vsucc)
          })

#' @rdname pvacsuccess-methods
#' @aliases pvacsuccess,numeric,vsucc.Logistic-method
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
