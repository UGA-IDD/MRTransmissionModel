#' Get the average age of infection at each time point
#'
#' @param sim.res object of type sim.results.SIR
#' @param ... ignored
#'
#' @return numeric vector of the average age of infection at each time point
#'
#' @include setClasses.R
#' @export
#' @docType methods
#' @rdname getAverageInfectionAge-methods
setGeneric("getAverageInfectionAge",
           function(sim.res, ...) standardGeneric("getAverageInfectionAge"))

#' @rdname getAverageInfectionAge-methods
#' @aliases getAverageInfectionAge,sim.results.SIR-method
setMethod("getAverageInfectionAge",
          "sim.results.SIR",
          function(sim.res) {
            age.mids <- (sim.res@age.class +
                           c(0,sim.res@age.class[2:length(sim.res@age.class)-1]))/2
            tmp <- sim.res[sim.res@i.inds,]*age.mids

            rc <- colSums(tmp)/
              colSums(sim.res[sim.res@i.inds,])
            return(rc)

          })

#' @rdname getAverageInfectionAge-methods
#' @aliases getAverageInfectionAge,sim.results.MSIRV-method
setMethod("getAverageInfectionAge",
          "sim.results.MSIRV",
          function(sim.res) {
            age.mids <- (sim.res@age.class +
                           c(0,sim.res@age.class[2:length(sim.res@age.class)-1]))/2
            tmp <- sim.res[sim.res@i.inds,]*age.mids

            rc <- colSums(tmp)/
              colSums(sim.res[sim.res@i.inds,])
            return(rc)

          })

