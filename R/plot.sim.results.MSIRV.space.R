#' Plot a sim.results.MSIRV.space object
#'
#' @param x object of class sim.results.SIR
#' @param ... ignored
#'
#' @importFrom graphics par lines mtext
#' @importFrom grDevices rgb col2rgb
#'
#' @return plot of sim results
#' @method plot sim.results.MSIRV.space
#' @export

plot.sim.results.MSIRV.space <- function(x, ...){


  # plotting parameters
  from <- 0
  to <- max(x@t)
  low.age <- 0
  high.age <- max(x@age.class)
  proportions <- FALSE
  plot.events <- TRUE

  result.all <- x
  orig.t <- x@t

  for (s in 1:result.all@n.subpops){

    #limit the data to only the subpopulation
    x <- result.all
    here <- which(x@subpop.class.label==s,arr.ind=TRUE)
    x@.Data <- x[here,x@t>=from & x@t<=to]

    #limit the time to the time of interest
    x@t <- x@t[x@t>=from & x@t<=to]

    #limit to the age classes of interest and the subpopulation of interest
    #the age classes will be the same regardless of the subpopulation, so I didn't specify the indeces for age.class vector
    age.classes <- c()
    for (i in 1:length(x@age.class)) {
      age.classes <- c(age.classes, rep(x@age.class[i],5))
    }
    #similarly the msirv indeces will be the same once we have shrunken the .Data object to the correct rows of the subpopulation
    last.index <- length(x@age.class)
    first.index <- last.index-(length(x@age.class))+1
    #make everything be in the age classes of interest
    x@.Data <- x[age.classes>=low.age & age.classes<=high.age,]
    x@m.inds <- x@m.inds[(first.index:last.index)[x@age.class>=low.age & x@age.class<=high.age]]
    x@s.inds <- x@s.inds[(first.index:last.index)[x@age.class>=low.age & x@age.class<=high.age]]
    x@i.inds <- x@i.inds[(first.index:last.index)[x@age.class>=low.age & x@age.class<=high.age]]
    x@r.inds <- x@r.inds[(first.index:last.index)[x@age.class>=low.age & x@age.class<=high.age]]
    x@v.inds <- x@v.inds[(first.index:last.index)[x@age.class>=low.age & x@age.class<=high.age]]
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

    par(mfcol = c(4,1), mar=c(3,4,2,2), oma=c(0,0,2,0))
    plot(tots@t, tots[tots@s.inds,]+tots[tots@m.inds,], type="b", xlab="t",
         ylab="M+S", cex=.75, ylim=c(0, max(tots[tots@s.inds,]+tots[tots@m.inds,])))
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
    mtext(paste("Subpopulation",s), line=0, side=3, outer=TRUE, cex=2)
  }

}
