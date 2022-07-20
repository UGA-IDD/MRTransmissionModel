#register a function for getting the seasonal adjustment for
#a time or series of times based on a seasonality object
#
#Parameter:
#     t.in.yr - time of the year in fraction of 1 year
#     season.obj - an object encoding the seasonality
#
#
#Returns:
#   the seasonal multiplier, or series thereof
setGeneric("get.seasonal.mult",
           function(t.in.yr, season.obj) standardGeneric("get.seasonal.mult"))


setMethod("get.seasonal.mult",
          c("numeric","seasonal.cosine"),
          function(t.in.yr, season.obj) {
            mult <- 1+ season.obj@amplitude *
              cos(2*season.obj@period*pi*(t.in.yr+season.obj@offset))
            return(mult)
          })

setMethod("get.seasonal.mult",
          c("numeric","seasonal.periodic.bspline"),
          function(t.in.yr, season.obj) {
            require(pomp)
            y <- periodic.bspline.basis(t.in.yr,nbasis=season.obj@nbasis,
                                        degree=season.obj@degree,period=season.obj@period)
            mult <- y%*%season.obj@parameters
            return(mult)
          })


setMethod("get.seasonal.mult",
          c("numeric","seasonal.piecewise.constant"),
          function(t.in.yr, season.obj) {
            mult <- season.obj@parameters[findInterval(t.in.yr%%1,season.obj@time.in.year)]
            return(mult)
          })


setMethod("get.seasonal.mult",
          c("numeric","seasonal.age.specific"),
          function(t.in.yr, season.obj) {
            mult.vec <- 1+ season.obj@amplitude *
              cos(2*season.obj@period*pi*(t.in.yr+season.obj@offset))

            nage <- length(season.obj@age.class)
            change.age <- which(season.obj@age.class>=season.obj@chosen.age[1] &
                                  season.obj@age.class<=season.obj@chosen.age[2])

            mult <- array(dim=c(nage,nage,length(t.in.yr)))
            for (j in 1:length(t.in.yr)) {
              mult[,,j] <- matrix(1,nage,nage)
              mult[,,j][change.age,change.age] <- mult.vec[j]
            }
            return(mult)
          })


setMethod("get.seasonal.mult",
          c("numeric","seasonal.age.specific.piecewise.constant"),
          function(t.in.yr, season.obj) {
            mult.vec <- season.obj@parameters[findInterval(t.in.yr%%1,season.obj@time.in.year)]

            nage <- length(season.obj@age.class)
            change.age <- which(season.obj@age.class>=season.obj@chosen.age[1] &
                                  season.obj@age.class<=season.obj@chosen.age[2])

            mult <- array(dim=c(nage,nage,length(t.in.yr)))
            for (j in 1:length(t.in.yr)) {
              mult[,,j] <- matrix(1,nage,nage)
              mult[,,j][change.age,change.age] <- mult.vec[j]
            }

            return(mult)
          })
