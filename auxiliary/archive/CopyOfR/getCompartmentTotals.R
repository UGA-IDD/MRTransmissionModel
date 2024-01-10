#' Get the total number of susceptible, infected, and recovered from sim.results.SIR
#'
#' @param sim.res sim result object
#' @param ... ignored
#'
#' @return a sim result object with no age classes
#'
#' @include setClasses.R
#'
#' @export
#' @docType methods
#' @rdname getCompartmentTotals-methods


setGeneric("getCompartmentTotals",
           function(sim.res, ...) standardGeneric("getCompartmentTotals"))


#' @rdname getCompartmentTotals-methods
#' @aliases getCompartmentTotals,sim.results.SIR-method
setMethod("getCompartmentTotals",
          "sim.results.SIR",
          function (sim.res) {
            rc <- new("sim.results.SIR",
                      nrow = 3,
                      ncol = ncol(sim.res),
                      t = sim.res@t,
                      s.inds = 1,
                      i.inds = 2,
                      r.inds = 3,
                      age.class = max(sim.res@age.class))
            rc[1,] <- colSums(sim.res[sim.res@s.inds,])
            rc[2,] <- colSums(sim.res[sim.res@i.inds,])
            rc[3,] <- colSums(sim.res[sim.res@r.inds,])

            rownames(rc) <- c("S","I","R")
            return(rc)
          })

#' @rdname getCompartmentTotals-methods
#' @aliases getCompartmentTotals,sim.results.MSIRV-method
setMethod("getCompartmentTotals",
          "sim.results.MSIRV",
          function (sim.res) {
            rc <- new("sim.results.MSIRV",
                      nrow = 5,
                      ncol = ncol(sim.res),
                      t = sim.res@t,
                      m.inds = 1,
                      s.inds = 2,
                      i.inds = 3,
                      r.inds = 4,
                      v.inds = 5,
                      age.class = max(sim.res@age.class))
            #sia.times = sim.res@sia.times,
            #trigger.times = sim.res@trigger.times)
            rc[rc@m.inds,] <- colSums(sim.res[sim.res@m.inds,])
            rc[rc@s.inds,] <- colSums(sim.res[sim.res@s.inds,])
            rc[rc@i.inds,] <- colSums(sim.res[sim.res@i.inds,])
            rc[rc@r.inds,] <- colSums(sim.res[sim.res@r.inds,])
            rc[rc@v.inds,] <- colSums(sim.res[sim.res@v.inds,])

            rownames(rc) <- c("M","S","I","R","V")
            return(rc)
          })
