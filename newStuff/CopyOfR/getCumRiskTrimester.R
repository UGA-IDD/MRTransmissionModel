#' Get cumulative probability of infection during pregnancy
#'
#' Calculate the cumulative probability of getting infected
#' over a period of n weeks from a sim.results.MSIRV object
#' to account for the fact that baby can get CRS if
#' mother gets sick during first trimester (1st three months)
#'
#' @param sim.res a sim.result.MSIRV object
#' @param nbiweeks chosen number of biweeks of risk
#' @param ... ignored
#'
#' @include setClasses.R
#' @return an object with age classes set by fertility age classes
#'
#' @export
#' @docType methods
#' @rdname getCumRiskTrimester-methods

setGeneric("getCumRiskTrimester",
           function(sim.res,nbiweeks, ...) standardGeneric("getCumRiskTrimester"))


#' @rdname getCumRiskTrimester-methods
#' @aliases getCumRiskTrimester,sim.results.MSIRV,numeric-method
setMethod("getCumRiskTrimester",
          c("sim.results.MSIRV","numeric"),
          function (sim.res, nbiweeks, ...) {

            n.age.cats <- length(sim.res@age.class)
            ntimesteps <- length(sim.res@t)

            rc <- matrix(NA,n.age.cats,ntimesteps)

            for (t in 2:(ntimesteps-nbiweeks)) {
              cumvals <- rowSums(sim.res[sim.res@i.inds,t:(t+nbiweeks-1)])
              cumprobs <- cumvals/pmax(sim.res[sim.res@s.inds,t-1],1)
              rc[,t] <- cumprobs
              #print(cbind(cumprobs,sim.res[sim.res@s.inds,t-1]))
              #this is the probabilty of getting infected, because to be an I you
              #had to first be an S at time t-1, and then the sum of I's are the number of S's who became infected
            }

            rc[1,] <- 0
            rc[,1] <- rc[,2]
            #print("here")
            rc[,(ntimesteps-nbiweeks):ntimesteps] <- rc[,ntimesteps-nbiweeks]

            return(rc)
          })
