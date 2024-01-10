#' Plot function for class sim.results.SIR
#'
#' @param x object of class sim.results.SIR
#' @param ... ignored
#'
#' @return plot of sim results
#' @method plot sim.results.SIR
#' @export
#'
plot.sim.results.SIR <- function (x,...) {

  par(mfcol = c(4,1))
  par(mar=c(3,2,2,2))

  tots <- getCompartmentTotals(x)
  plot(tots@t,tots[1,], type="b", xlab="t",
       ylab="S", cex=.75,...)
  plot(tots@t,tots[2,], type="b", xlab="t",
       ylab="I", cex=.75, ...)
  plot(tots@t,tots[3,], type="b", xlab="t",
       ylab="R", cex=.75, ...)
  plot(x@t,getAverageInfectionAge(x) , type="b", xlab="t",
       ylab="Age", cex=.75, ...)
}


