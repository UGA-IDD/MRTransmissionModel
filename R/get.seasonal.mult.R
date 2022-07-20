#' Get seasonal adjustment for times based on a seasonality object
#'
#'
#' @param t.in.yr time of the year in fraction of 1 year
#' @param season.obj an object encoding the seasonality
#'
#' @import pomp
#' @include setClasses.R
#' @return the seasonal multiplier, or series thereof
#'
#' @export
#' @docType methods
#' @rdname get.seasonal.mult-methods

setGeneric("get.seasonal.mult",
           function(t.in.yr, season.obj) standardGeneric("get.seasonal.mult"))

#' @rdname get.seasonal.mult-methods
#' @aliases get.seasonal.mult,numeric,seasonal.cosine-method
setMethod("get.seasonal.mult",
          c("numeric","seasonal.cosine"),
          function(t.in.yr, season.obj) {
            mult <- 1+ season.obj@amplitude *
              cos(2*season.obj@period*pi*(t.in.yr+season.obj@offset))
            return(mult)
          })

#' @rdname get.seasonal.mult-methods
#' @aliases get.seasonal.mult,numeric,seasonal.periodic.bspline-method
setMethod("get.seasonal.mult",
          c("numeric","seasonal.periodic.bspline"),
          function(t.in.yr, season.obj) {

            y <- pomp::periodic.bspline.basis(t.in.yr,nbasis=season.obj@nbasis,
                                        degree=season.obj@degree,period=season.obj@period)
            mult <- y%*%season.obj@parameters
            return(mult)
          })

#' @rdname get.seasonal.mult-methods
#' @aliases get.seasonal.mult,numeric,seasonal.piecewise.constant-method
setMethod("get.seasonal.mult",
          c("numeric","seasonal.piecewise.constant"),
          function(t.in.yr, season.obj) {
            mult <- season.obj@parameters[findInterval(t.in.yr%%1,season.obj@time.in.year)]
            return(mult)
          })

#' @rdname get.seasonal.mult-methods
#' @aliases get.seasonal.mult,numeric,seasonal.age.specific-method
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

#' @rdname get.seasonal.mult-methods
#' @aliases get.seasonal.mult,numeric,seasonal.age.specific.piecewise.constant-method
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
