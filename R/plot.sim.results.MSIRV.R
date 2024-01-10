#' Plot method for class "sim.results.MSIRV"
#'
#' @param x object of class "sim.results.MSIRV"
#' @param main character, plot title
#' @param ... ignored
#'
#' @importFrom graphics par lines
#' @importFrom grDevices rgb col2rgb
#'
#' @return plot object
#' @method plot sim.results.MSIRV
#' @export

plot.sim.results.MSIRV <- function (x, main,...) {

  # plotting parameters
  from=0
  to=max(x@t)
  low.age = 0
  high.age = max(x@age.class)
  proportions = FALSE
  plot.events = TRUE

  par(mfcol = c(4,1))
  par(mar=c(3,4,2,2))

  orig.t <- x@t

  x@.Data <- x[,x@t>=from & x@t<=to]
  x@t <- x@t[x@t>=from & x@t<=to]

  age.classes <- c()
  for (i in 1:length(x@age.class)) {
    age.classes <- c(age.classes, rep(x@age.class[i],5))
  }

  #make everything be in the age classes of interest
  x@.Data <- x[age.classes>=low.age & age.classes<=high.age,]
  x@m.inds <- x@m.inds[x@age.class>=low.age & x@age.class<=high.age]
  x@s.inds <- x@s.inds[x@age.class>=low.age & x@age.class<=high.age]
  x@i.inds <- x@i.inds[x@age.class>=low.age & x@age.class<=high.age]
  x@r.inds <- x@r.inds[x@age.class>=low.age & x@age.class<=high.age]
  x@v.inds <- x@v.inds[x@age.class>=low.age & x@age.class<=high.age]
  x@age.class <- x@age.class[x@age.class>=low.age & x@age.class<=high.age]


  tots <- getCompartmentTotals(x)

  if (proportions) {
    tots <- toProportions(tots)
  }

  plt.events <- function() {
    if (!plot.events) return()

    sia.times <- orig.t[which(x@sia.times>0)]
    routine.times <- orig.t[which(x@routine.intro>0)]

    if (length(sia.times)>0) {
      for (i in 1:length(sia.times)) {
        lines(rep(sia.times[i],2),c(0, 10^10),
              col=rgb(t(col2rgb("cornflowerblue")/256), alpha=.45),
              lty=3)

      }
    }


    if (length(routine.times)>0) {
      for (i in 1:length(routine.times)) {
        lines(rep(routine.times[i],2),c(0, 10^10),
              col=rgb(t(col2rgb("coral")/256), alpha=.45),
              lty=3)
      }
    }

  }

  plot(tots@t, tots[tots@s.inds,]+tots[tots@m.inds,], type="b", xlab="t",
       ylab="M+S", cex=.75, ylim=c(0, max(tots[tots@s.inds,]+tots[tots@m.inds,])), main=main)
  lines(tots@t, tots[tots@s.inds,], lty=2, col="green")
  lines(tots@t, tots[tots@m.inds,], lty=3, col="blue")
  plt.events()
  plot(tots@t,tots[tots@i.inds,], type="b", xlab="t", ylab="I", cex=.75)
  plt.events()
  plot(tots@t,tots[tots@r.inds,]+tots[tots@v.inds,], type="b", xlab="t",
       ylab="R+V", cex=.75,ylim=c(0, max(tots[tots@r.inds,]+tots[tots@v.inds,])))
  lines(tots@t, tots[tots@r.inds,], lty=2, col="green")
  lines(tots@t, tots[tots@v.inds,], lty=3, col="blue")
  plt.events()
  plot(x@t,getAverageInfectionAge(x) , type="b", xlab="t", ylab="Age", cex=.75)
  plt.events()

}

