#### Methods to Manipulate or Plot Experiment Output ####

#Plot a sim.results.SIR object
#does a 4 panel plot that shows M+S,I,R+V and average age
#
#Parameters -
#     x - the sim.results.MSIRV object
#     y - ignored
#
setMethod("plot",
          "sim.results.MSIRV",
          function (x,y,from=0, to=max(x@t), low.age = 0,
                    high.age = max(x@age.class), proportions = FALSE,
                    plot.events = TRUE,...) {
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

          })

#Plot a sim.results.SIR object
#does a 4 panel plot that shows S,I,R and average age
#
#Parameters -
#     x - the sim.results.SIR object
#     y - ignored
#
setMethod("plot",
          "sim.results.SIR",
          function (x,y,...) {
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
          })


#Get the total number of susceptible, infected, and
#recovered from and sim.results.SIR
#
#returns -
# a sim results object with no age classes
setGeneric("getCompartmentTotals",
           function(sim.res, ...) standardGeneric("getCompartmentTotals"))
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

#Method for getting the average age of infection at each time point.
#takes in a sim.results.SIR object and returns a numeric vector of
#the average age of infection at each time point
setGeneric("getAverageInfectionAge",
           function(sim.res, ...) standardGeneric("getAverageInfectionAge"))
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


#Convert the results to proportions of each age class in each
#compartment, rather than the absolute numbers
setGeneric("toProportions", function(ob, ...) standardGeneric("toProportions"))
setMethod("toProportions",
          "sim.results.MSIRV",
          function (ob) {
            ac <- c()
            for (i in ob@age.class) {
              ac <- c(ac, rep(i,5))
            }



            #print(cbind(round(ob@t,2),t(round(ob[ac==50,],2)),colSums(ob[ac==50,])))
            for (i in ob@age.class) {
              div <- t(matrix(rep(colSums(ob[ac==i,]),5), ncol=5))
              ob[ac==i,] <- ob[ac==i,]/div
            }

            #print(cbind(round(ob@t,2),t(round(ob[ac==50,],2)),colSums(ob[ac==50,])))


            return(ob)
          })

#Function to get the cumulative prob of getting infected over a period of nweeks from a sim.results.MSIRV
# to account for the fact that can get have CRS if get sick during first trimester (1st three months)
#
#Parameters -
#       a sim.result.MSIRV
#       nbiweeks - chosen no biweeks of risk
#
#Returns -
#       an object with age classes set by fertility age classes
setGeneric("getCumRiskTrimester",
           function(sim.res,nbiweeks, ...) standardGeneric("getCumRiskTrimester"))

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
