#' Methods to run experiment
#'
#'
#' @param exper the experiment object
#' @param ... additional arguments
#'
#' @include setClasses.R
#' @return an experiment result object
#'
#' @export
#' @docType methods
#' @rdname run-methods
#'

setGeneric("run", function(exper,...) standardGeneric("run"))

#' @param rescale.WAIFW logical; should we rescale the WAIFW matrix? Set to FALSE if experiment is not a fully susceptible state
#' @rdname run-methods
#' @aliases run,experiment.updatedmog,ANY-method
setMethod("run",
          "experiment.updatedemog",
          function(exper, rescale.WAIFW=TRUE, ...) {

            #print(rescale.WAIFW)
            state <- exper@state.t0

            #rescale the WAIFW if a specific R0 is specified
            if (rescale.WAIFW & length(exper@R0)>0) {
              #print("RESCALING!!!")
              exper@trans@waifw <- scaleWAIFW(exper@R0,
                                               state,exper@trans@waifw,
                                               frequency.dep=exper@trans@frequency.dep,
                                               suscept.state=exper@trans@s.inds[1])
            }

            #get the number of time steps in the experiment
            numTimeSteps <- round((exper@t.max-exper@t.min)/exper@step.size)+1

            if (!is.null(exper@season.obj)) {
              mults <- get.seasonal.mult(exper@t0.doy/365+(1:numTimeSteps-1)*exper@step.size,
                                         exper@season.obj)
            } else {
              mults <- rep(1,numTimeSteps)
            }

            #make a temporary transmission object
            tmp.trans <- exper@trans

            #hold the states as we walk through
            rc <- matrix(ncol = numTimeSteps, nrow = nrow(state))
            rc[,1] <- state

            #need output vector for births overtime
            births.each.timestep <- growth.rate.each.timestep <- rep(NA, numTimeSteps)
            births.each.timestep[1] <- tmp.trans@birth.rate

            #output vector for when SIAs were administered - default in this experiment is 0
            sia.times <- routine.intro <- rep(0,numTimeSteps)

            for (t in 2:numTimeSteps) {

              #scale the waifw by seasonality
              if (!is.array(mults)){
                tmp.trans@waifw <- exper@trans@waifw*mults[t]
              } else {
                tmp.trans@waifw <- exper@trans@waifw*mults[,,t]
              }

              #if exper@pop.rescale.each.timestep[t]==NaN then tmp.trans@pop.rescale==NaN and will not rescale
              #otherwise if any number other than NaN it will rescale
              #or if exper@pop.rescale.each.timestep not completed , then numeric(0) and any index [t] is NA
              if (!is.na(exper@pop.rescale.each.timestep[t])) {state <- exper@pop.rescale.each.timestep[t] *(state/sum(state))}

              #put in the correct birth rate for that time-step, if it varies
              if (length(exper@births.per.1000.each.timestep)>1) {
                # Read 'birth rate'
                tmp.trans@birth.rate = (exper@births.per.1000.each.timestep[t]*sum(state)/1000)
              } else {
                tmp.trans@birth.rate <- exper@trans@birth.rate*sum(rc[,t-1])/sum(exper@state.t0)
              }
              births.each.timestep[t] <-  tmp.trans@birth.rate

              #put in time appropriate survival rate otherwise it uses original surv.matrix set up for trans object and keeps constant over time
              if (!is.na(exper@surv.each.timestep[1,1])) {
                tmp.trans@age.surv.matrix <- ExtractAgeSpecificSurvivalMatrix(tmp.trans=tmp.trans, maternal.obj=exper@maternal.obj,
                                                                              surv.at.timestep.t=exper@surv.each.timestep[,t])
              }

              #put in the correct SIAs
              #tmp.trans@sia.vac <- exper@sia.obj[,((t-1)%%ncol(exper@sia.obj))+1]
              #if (sum(tmp.trans@sia.vac)>0) {sia.times[t] <- 1}

              N0 <- sum(state)
              state <- next.ID.state(state, tmp.trans)
              NT <- sum(state)
              growth.rate.each.timestep[t] <- log(NT/N0) #instaneous biweekly growth rate

              #print(t)
              #print(dim(state))
              rc [,t] <- state

            }


            rc <- new("sim.results.MSIRV.update.demog",
                      data=rc,
                      m.inds = exper@trans@m.inds,
                      s.inds = exper@trans@s.inds,
                      i.inds = exper@trans@i.inds,
                      r.inds = exper@trans@r.inds,
                      v.inds = exper@trans@v.inds,
                      t = exper@t.min+(1:numTimeSteps-1)*exper@step.size,
                      age.class = exper@trans@age.class,
                      births.each.timestep = births.each.timestep,
                      growth.rate.each.timestep = growth.rate.each.timestep,
                      routine.intro = routine.intro,
                      sia.times = sia.times)

            rc <- new("experiment.result",
                      experiment.def = exper,
                      result = rc)


            return(rc)
          }
)

#' @param rescale.WAIFW logical; should we rescale WAIFW matrix
#' @rdname run-methods
#' @aliases run,experiment.updatedmog.vaccinationchange,ANY-method
setMethod("run",
          "experiment.updatedemog.vaccinationchange",
          function(exper, rescale.WAIFW=TRUE, ...) {

            #print(rescale.WAIFW)
            state <- exper@state.t0

            #rescale the WAIFW if a specific R0 is specified
            if (rescale.WAIFW & length(exper@R0)>0) {
              #print("RESCALING!!!")
              exper@trans@waifw <- scaleWAIFW(exper@R0,
                                               state,exper@trans@waifw,
                                               frequency.dep=exper@trans@frequency.dep,
                                               suscept.state=exper@trans@s.inds[1])
            }

            #get the number of time steps in the experiment
            numTimeSteps <- round((exper@t.max-exper@t.min)/exper@step.size)+1

            if (!is.null(exper@season.obj)) {
              mults <- get.seasonal.mult(exper@t0.doy/365+(1:numTimeSteps-1)*exper@step.size,
                                         exper@season.obj)
            } else {
              mults <- rep(1,numTimeSteps)
            }

            #make a temporary transmission object
            tmp.trans <- exper@trans

            #hold the states as we walk through
            rc <- matrix(ncol = numTimeSteps, nrow = nrow(state))
            rc[,1] <- state

            #need output vector for births over time
            births.each.timestep <- growth.rate.each.timestep <- rep(NA, numTimeSteps)
            births.each.timestep[1] <- tmp.trans@birth.rate

            #generate the age and time specific vaccination matrix
            routine <- get.routine.time.age.specific(time.step= exper@step.size*12,
                                                     age.classes=exper@trans@age.class,
                                                     time.specific.MR1cov=exper@time.specific.MR1cov,
                                                     age.min.MR1=exper@time.specific.min.age.MR1,
                                                     age.max.MR1=exper@time.specific.max.age.MR1,
                                                     time.specific.MR2cov=exper@time.specific.MR2cov,
                                                     age.min.MR2=exper@time.specific.min.age.MR2,
                                                     age.max.MR2=exper@time.specific.max.age.MR2,
                                                     obj.vcdf.MR1=exper@obj.vcdf.MR1,
                                                     obj.vcdf.MR2=exper@obj.vcdf.MR2,
                                                     obj.prob.vsucc=exper@obj.prob.vsucc,
                                                     MR1MR2correlation=FALSE)
            routine.intro <- rep(0, numTimeSteps)
            if (any(exper@time.specific.MR1cov!=0)) routine.intro[min(which(exper@time.specific.MR1cov>0))*(1/exper@step.size)+1] <- 1
            if (any(exper@time.specific.MR2cov!=0)) routine.intro[min(which(exper@time.specific.MR2cov>0))*(1/exper@step.size)+1] <- 1
            index.routine.vacc <- c(1,rep(1:nrow(routine$age.time.specific.routine), each=(numTimeSteps-1)/exper@t.max))
            if (length(index.routine.vacc)<numTimeSteps) index.routine.vacc[(length(index.routine.vacc)+1):numTimeSteps] <- index.routine.vacc[length(index.routine.vacc)]

            SIA <- get.sia.time.age.specific(age.classes=exper@trans@age.class,
                                             time.specific.SIAcov=exper@time.specific.SIAcov,
                                             age.min.sia=exper@time.specific.min.age.SIA,
                                             age.max.sia=exper@time.specific.max.age.SIA,
                                             obj.prob.vsucc=exper@obj.prob.vsucc)
            index.sia.vacc <- rep(NA,numTimeSteps)
            year.sia <- which(exper@time.specific.SIAcov!=0)
            index.sia.vacc[(year.sia-1)*(numTimeSteps-1)/exper@t.max + round((exper@sia.timing.in.year*(numTimeSteps-1)/exper@t.max))] <-  year.sia #minus 1 because adding the sia.timing
            sia.times <- ifelse(!is.na(index.sia.vacc), 1, 0)

            #need output vectors for primary vaccination failure over time
            MR1.fail.each.timestep <- MR2.fail.each.timestep <- SIA.fail.each.timestep <- rep(0, numTimeSteps)
            MR1.fail.each.timestep[1] <- routine$prop.fail.MR1[1]
            MR2.fail.each.timestep[1] <- routine$prop.fail.MR2[1]
            SIA.fail.each.timestep[1] <- 0

            for (t in 2:numTimeSteps) {

              #introduction rates
              if (length(exper@intro.rate)>1) tmp.trans@introduction.rate <- rep(exper@intro.rate[t],exper@trans@n.age.class)

              #scale the waifw by seasonality
              if (!is.array(mults)){
                tmp.trans@waifw <- exper@trans@waifw*mults[t]
              } else {
                tmp.trans@waifw <- exper@trans@waifw*mults[,,t]
              }

              #if exper@pop.rescale.each.timestep[t]==NaN then tmp.trans@pop.rescale==NaN and will not rescale
              #otherwise if any number other than NaN it will rescale
              #or if exper@pop.rescale.each.timestep not completed , then numeric(0) and any index [t] is NA
              if (!is.na(exper@pop.rescale.each.timestep[t])) {state <- exper@pop.rescale.each.timestep[t] *(state/sum(state))}

              #put in the correct birth rate for that time-step, if it varies
              if (length(exper@births.per.1000.each.timestep)>1) {
                # Read 'birth rate'
                tmp.trans@birth.rate = (exper@births.per.1000.each.timestep[t]*sum(state)/1000)
              } else {
                tmp.trans@birth.rate <- exper@trans@birth.rate*sum(rc[,t-1])/sum(exper@state.t0)
              }
              births.each.timestep[t] <-  tmp.trans@birth.rate

              #put in time appropriate survival rate otherwise it uses original surv.matrix set up for trans object and keeps constant over time
              if (!is.na(exper@surv.each.timestep[1,1])) {
                tmp.trans@age.surv.matrix <- ExtractAgeSpecificSurvivalMatrix(tmp.trans=tmp.trans, maternal.obj=exper@maternal.obj,
                                                                              surv.at.timestep.t=exper@surv.each.timestep[,t])
              }

              ##put in correct vaccination coverage
              routine.vacc.prob <- routine$age.time.specific.routine
              sia.vacc.prob <- SIA$age.time.specific.SIA
              if (!is.na(index.sia.vacc[t])){ #if SIA
                tmp.trans@vac.per@pvacc.in.age.class <-
                  routine.vacc.prob[index.routine.vacc[t],] +
                  sia.vacc.prob[index.sia.vacc[t],] - #xxamy changed from addition to subtraction
                  (routine.vacc.prob[index.routine.vacc[t],]*sia.vacc.prob[index.sia.vacc[t],])
                #stow output
                MR1.fail.each.timestep[t] <- routine$prop.fail.MR1[index.routine.vacc[t]]
                MR2.fail.each.timestep[t] <- routine$prop.fail.MR2[index.routine.vacc[t]]
                SIA.fail.each.timestep[t] <- SIA$prop.fail.SIA[index.sia.vacc[t]]

              } else { #if no SIA
                tmp.trans@vac.per@pvacc.in.age.class <- routine.vacc.prob[index.routine.vacc[t],]
                #stow output
                MR1.fail.each.timestep[t] <- routine$prop.fail.MR1[index.routine.vacc[t]]
                MR2.fail.each.timestep[t] <- routine$prop.fail.MR2[index.routine.vacc[t]]
              }


              N0 <- sum(state)
              state <- next.ID.state(state, tmp.trans)
              NT <- sum(state)
              growth.rate.each.timestep[t] <- log(NT/N0) #instantaneous biweekly growth rate

              #print(t)
              #print(dim(state))
              rc [,t] <- state

            }


            rc <- new("sim.results.MSIRV.update.demog.vaccine.change",
                      data=rc,
                      m.inds = exper@trans@m.inds,
                      s.inds = exper@trans@s.inds,
                      i.inds = exper@trans@i.inds,
                      r.inds = exper@trans@r.inds,
                      v.inds = exper@trans@v.inds,
                      t = exper@t.min+(1:numTimeSteps-1)* exper@step.size,
                      age.class = exper@trans@age.class,
                      births.each.timestep = births.each.timestep,
                      growth.rate.each.timestep = growth.rate.each.timestep,
                      MR1.fail.each.timestep = MR1.fail.each.timestep,
                      MR2.fail.each.timestep = MR2.fail.each.timestep,
                      SIA.fail.each.timestep = SIA.fail.each.timestep,
                      routine.intro = routine.intro,
                      sia.times = sia.times)

            rc <- new("experiment.result",
                      experiment.def = exper,
                      result = rc)


            return(rc)
          }
)

#' @rdname run-methods
#' @aliases run,experiment.updatedemog.vaccinationchange.vaccinationcorrelation,ANY-method
setMethod("run",
          "experiment.updatedemog.vaccinationchange.vaccinationcorrelation",
          function(exper, rescale.WAIFW=TRUE, ...) {

            state <- exper@state.t0

            if (rescale.WAIFW & length(exper@R0)>0) {
              exper@trans@waifw <- scaleWAIFW(exper@R0,
                                               state, exper@trans@waifw,
                                               frequency.dep=exper@trans@frequency.dep,
                                               suscept.state=exper@trans@s.inds[1])
            }

            numTimeSteps <- round((exper@t.max-exper@t.min)/exper@step.size)+1

            if (!is.null(exper@season.obj)) {
              mults <- get.seasonal.mult(exper@t0.doy/365+(1:numTimeSteps-1)*exper@step.size,
                                         exper@season.obj)
            } else {
              mults <- rep(1, numTimeSteps)
            }

            tmp.trans <- exper@trans

            rc <- matrix(ncol=numTimeSteps, nrow=nrow(state))
            rc[,1] <- state

            births.each.timestep <- growth.rate.each.timestep <- rep(NA, numTimeSteps)
            births.each.timestep[1] <- tmp.trans@birth.rate

            routine <- get.routine.time.age.specific(
              time.step              = exper@step.size*12,
              age.classes            = exper@trans@age.class,
              time.specific.MR1cov  = exper@time.specific.MR1cov,
              age.min.MR1            = exper@time.specific.min.age.MR1,
              age.max.MR1            = exper@time.specific.max.age.MR1,
              time.specific.MR2cov  = exper@time.specific.MR2cov,
              age.min.MR2            = exper@time.specific.min.age.MR2,
              age.max.MR2            = exper@time.specific.max.age.MR2,
              obj.vcdf.MR1           = exper@obj.vcdf.MR1,
              obj.vcdf.MR2           = exper@obj.vcdf.MR2,
              obj.prob.vsucc         = exper@obj.prob.vsucc,
              MR1MR2correlation      = exper@MR1MR2correlation)

            routine.intro <- rep(0, numTimeSteps)
            if (any(exper@time.specific.MR1cov!=0)) routine.intro[min(which(exper@time.specific.MR1cov>0))*(1/exper@step.size)+1] <- 1
            if (any(exper@time.specific.MR2cov!=0)) routine.intro[min(which(exper@time.specific.MR2cov>0))*(1/exper@step.size)+1] <- 1
            index.routine.vacc <- c(1, rep(1:nrow(routine$age.time.specific.routine), each=(numTimeSteps-1)/exper@t.max))
            if (length(index.routine.vacc)<numTimeSteps) index.routine.vacc[(length(index.routine.vacc)+1):numTimeSteps] <- index.routine.vacc[length(index.routine.vacc)]

            SIA <- get.sia.time.age.specific(
              age.classes          = exper@trans@age.class,
              time.specific.SIAcov = exper@time.specific.SIAcov,
              age.min.sia          = exper@time.specific.min.age.SIA,
              age.max.sia          = exper@time.specific.max.age.SIA,
              obj.prob.vsucc       = exper@obj.prob.vsucc)

            index.sia.vacc <- rep(NA, numTimeSteps)
            year.sia <- which(exper@time.specific.SIAcov!=0)
            index.sia.vacc[(year.sia-1)*(numTimeSteps-1)/exper@t.max + round((exper@sia.timing.in.year*(numTimeSteps-1)/exper@t.max))] <- year.sia
            sia.times <- ifelse(!is.na(index.sia.vacc), 1, 0)

            MR1.fail.each.timestep <- MR2.fail.each.timestep <- SIA.fail.each.timestep <- rep(0, numTimeSteps)
            MR1.fail.each.timestep[1] <- routine$prop.fail.MR1[1]
            MR2.fail.each.timestep[1] <- routine$prop.fail.MR2[1]
            SIA.fail.each.timestep[1] <- 0

            for (t in 2:numTimeSteps) {

              if (length(exper@intro.rate)>1) tmp.trans@introduction.rate <- rep(exper@intro.rate[index.routine.vacc[t]], exper@trans@n.age.class)

              if (!is.array(mults)){
                tmp.trans@waifw <- exper@trans@waifw*mults[t]
              } else {
                tmp.trans@waifw <- exper@trans@waifw*mults[,,t]
              }

              if (!is.na(exper@pop.rescale.each.timestep[t])) {state <- exper@pop.rescale.each.timestep[t]*(state/sum(state))}

              if (length(exper@births.per.1000.each.timestep)>1) {
                tmp.trans@birth.rate <- exper@births.per.1000.each.timestep[t]*sum(state)/1000
              } else {
                tmp.trans@birth.rate <- exper@trans@birth.rate*sum(rc[,t-1])/sum(exper@state.t0)
              }
              births.each.timestep[t] <- tmp.trans@birth.rate

              if (!is.na(exper@surv.each.timestep[1,1])) {
                tmp.trans@age.surv.matrix <- ExtractAgeSpecificSurvivalMatrix(
                  tmp.trans=tmp.trans, maternal.obj=exper@maternal.obj,
                  surv.at.timestep.t=exper@surv.each.timestep[,t])
              }

              if (!is.na(index.sia.vacc[t])) {
                r <- routine$age.time.specific.routine[index.routine.vacc[t],]
                s <- SIA$age.time.specific.SIA[index.sia.vacc[t],]
                rho.vec <- rep(0, exper@trans@n.age.class)
                MR1.age.range <- exper@time.specific.min.age.MR1[index.routine.vacc[t]]:exper@time.specific.max.age.MR1[index.routine.vacc[t]]
                MR2.age.range <- exper@time.specific.min.age.MR2[index.routine.vacc[t]]:exper@time.specific.max.age.MR2[index.routine.vacc[t]]
                rho.vec[MR1.age.range] <- exper@MR1SIAcorrelation
                rho.vec[MR2.age.range] <- exper@MR2SIAcorrelation
                tmp.trans@vac.per@pvacc.in.age.class <-
                  r + s - r*s - rho.vec * sqrt(pmax(0, r*(1-r)*s*(1-s)))
                MR1.fail.each.timestep[t] <- routine$prop.fail.MR1[index.routine.vacc[t]]
                MR2.fail.each.timestep[t] <- routine$prop.fail.MR2[index.routine.vacc[t]]
                SIA.fail.each.timestep[t] <- SIA$prop.fail.SIA[index.sia.vacc[t]]
              } else {
                tmp.trans@vac.per@pvacc.in.age.class <- routine$age.time.specific.routine[index.routine.vacc[t],]
                MR1.fail.each.timestep[t] <- routine$prop.fail.MR1[index.routine.vacc[t]]
                MR2.fail.each.timestep[t] <- routine$prop.fail.MR2[index.routine.vacc[t]]
              }

              N0 <- sum(state)
              state <- next.ID.state(state, tmp.trans)
              NT <- sum(state)
              growth.rate.each.timestep[t] <- log(NT/N0)

              rc[,t] <- state
            }

            rc <- new("sim.results.MSIRV.update.demog.vaccine.change",
                      data=rc,
                      m.inds = exper@trans@m.inds,
                      s.inds = exper@trans@s.inds,
                      i.inds = exper@trans@i.inds,
                      r.inds = exper@trans@r.inds,
                      v.inds = exper@trans@v.inds,
                      t = exper@t.min+(1:numTimeSteps-1)*exper@step.size,
                      age.class = exper@trans@age.class,
                      births.each.timestep = births.each.timestep,
                      growth.rate.each.timestep = growth.rate.each.timestep,
                      MR1.fail.each.timestep = MR1.fail.each.timestep,
                      MR2.fail.each.timestep = MR2.fail.each.timestep,
                      SIA.fail.each.timestep = SIA.fail.each.timestep,
                      routine.intro = routine.intro,
                      sia.times = sia.times)

            rc <- new("experiment.result",
                      experiment.def = exper,
                      result = rc)

            return(rc)
          }
)

#' @rdname run-methods
#' @aliases run,experiment.updatedemog.vaccinationchange.vaccinationcorrelation.outbreakresponse,ANY-method
setMethod(
  "run",
  "experiment.updatedemog.vaccinationchange.vaccinationcorrelation.outbreakresponse",
  function(exper, rescale.WAIFW = TRUE) {

    state <- exper@state.t0

    if (rescale.WAIFW && length(exper@R0) > 0) {
      exper@trans@waifw <- scaleWAIFW(exper@R0,
                                       state, exper@trans@waifw,
                                       frequency.dep = exper@trans@frequency.dep,
                                       suscept.state = exper@trans@s.inds[1])
    }

    numTimeSteps <- round((exper@t.max - exper@t.min) / exper@step.size) + 1

    if (!is.null(exper@season.obj)) {
      mults <- get.seasonal.mult(
        exper@t0.doy / 365 + (1:numTimeSteps - 1) * exper@step.size,
        exper@season.obj)
    } else {
      mults <- rep(1, numTimeSteps)
    }

    tmp.trans <- exper@trans

    rc    <- matrix(ncol = numTimeSteps, nrow = nrow(state))
    rc[, 1] <- state

    births.each.timestep <- growth.rate.each.timestep <- rep(NA, numTimeSteps)
    births.each.timestep[1] <- tmp.trans@birth.rate

    # --- Routine vaccination schedule (unchanged from parent class) --

    routine <- get.routine.time.age.specific(
      time.step             = exper@step.size * 12,
      age.classes           = exper@trans@age.class,
      time.specific.MR1cov  = exper@time.specific.MR1cov,
      age.min.MR1           = exper@time.specific.min.age.MR1,
      age.max.MR1           = exper@time.specific.max.age.MR1,
      time.specific.MR2cov  = exper@time.specific.MR2cov,
      age.min.MR2           = exper@time.specific.min.age.MR2,
      age.max.MR2           = exper@time.specific.max.age.MR2,
      obj.vcdf.MR1          = exper@obj.vcdf.MR1,
      obj.vcdf.MR2          = exper@obj.vcdf.MR2,
      obj.prob.vsucc        = exper@obj.prob.vsucc,
      MR1MR2correlation     = exper@MR1MR2correlation)

    routine.intro <- rep(0, numTimeSteps)
    if (any(exper@time.specific.MR1cov != 0))
      routine.intro[min(which(exper@time.specific.MR1cov > 0)) * (1 / exper@step.size) + 1] <- 1
    if (any(exper@time.specific.MR2cov != 0))
      routine.intro[min(which(exper@time.specific.MR2cov > 0)) * (1 / exper@step.size) + 1] <- 1

    index.routine.vacc <- c(1, rep(1:nrow(routine$age.time.specific.routine),
                                   each = (numTimeSteps - 1) / exper@t.max))
    if (length(index.routine.vacc) < numTimeSteps)
      index.routine.vacc[(length(index.routine.vacc) + 1):numTimeSteps] <-
        index.routine.vacc[length(index.routine.vacc)]

    # --- Scheduled SIA campaigns (unchanged from parent class) --

    SIA <- get.sia.time.age.specific(
      age.classes          = exper@trans@age.class,
      time.specific.SIAcov = exper@time.specific.SIAcov,
      age.min.sia          = exper@time.specific.min.age.SIA,
      age.max.sia          = exper@time.specific.max.age.SIA,
      obj.prob.vsucc       = exper@obj.prob.vsucc)

    index.sia.vacc <- rep(NA, numTimeSteps)
    year.sia <- which(exper@time.specific.SIAcov != 0)
    index.sia.vacc[(year.sia - 1) * (numTimeSteps - 1) / exper@t.max +
                   round(exper@sia.timing.in.year * (numTimeSteps - 1) / exper@t.max)] <- year.sia
    sia.times <- ifelse(!is.na(index.sia.vacc), 1, 0)

    MR1.fail.each.timestep <- MR2.fail.each.timestep <- SIA.fail.each.timestep <-
      rep(0, numTimeSteps)
    MR1.fail.each.timestep[1] <- routine$prop.fail.MR1[1]
    MR2.fail.each.timestep[1] <- routine$prop.fail.MR2[1]

    # --- OBR pre-computation ---

    or.times            <- rep(0, numTimeSteps)
    last.or.t           <- -Inf
    pending.trigger.t   <- NA_integer_
    pending.case.by.age <- numeric(length(exper@trans@age.class))
    age.classes         <- exper@trans@age.class

    # Row indices of I compartment for the trigger age group
    if (is.na(exper@or.trigger.age.lower) || is.na(exper@or.trigger.age.upper)) {
      trigger.age.pos <- seq_along(age.classes)
    } else {
      trigger.age.pos <- which(age.classes > exper@or.trigger.age.lower &
                               age.classes <= exper@or.trigger.age.upper)
    }
    trigger.i.inds <- exper@trans@i.inds[trigger.age.pos]

    # Reporting rate: trigger age group (for metric) and all ages (for response accumulation)
    or.rep.rate <- if (length(exper@or.reporting.rate) == 1) {
      rep(exper@or.reporting.rate, length(trigger.age.pos))
    } else {
      exper@or.reporting.rate[trigger.age.pos]
    }
    or.rep.rate.all <- if (length(exper@or.reporting.rate) == 1) {
      rep(exper@or.reporting.rate, length(age.classes))
    } else {
      exper@or.reporting.rate
    }

    # Steps per year (for calendar month mapping in observational model and accumulation)
    steps.per.year.or <- round(1 / exper@step.size)

    # Earliest time step at which OBR trigger can be evaluated
    or.min.t <- max(exper@or.total.delay + exper@or.trigger.window,
                    exper@or.start.timestep)

    # "confirmed" trigger mode: pre-allocate per-timestep confirmed case vector
    if (exper@or.trigger.mode == "confirmed") {
      confirmed.trigger <- numeric(numTimeSteps)
    }

    # Observational model: pre-allocate 6 age x timestep output matrices.
    # Computed every step when Se, Sp, and non_meas are available.
    has.obs.model <- nrow(exper@or.non.meas.cases.by.age.month) > 0 &&
                     !is.na(exper@or.Se) && !is.na(exper@or.Sp)
    n.age.all     <- length(age.classes)
    obs.TP          <- matrix(0, nrow = n.age.all, ncol = numTimeSteps)
    obs.FN_test     <- matrix(0, nrow = n.age.all, ncol = numTimeSteps)
    obs.FP_test     <- matrix(0, nrow = n.age.all, ncol = numTimeSteps)
    obs.TN          <- matrix(0, nrow = n.age.all, ncol = numTimeSteps)
    obs.TP_clinical <- matrix(0, nrow = n.age.all, ncol = numTimeSteps)
    obs.FP_clinical <- matrix(0, nrow = n.age.all, ncol = numTimeSteps)

    for (t in 2:numTimeSteps) {

      if (length(exper@intro.rate) > 1)
        tmp.trans@introduction.rate <- rep(exper@intro.rate[index.routine.vacc[t]],
                                           exper@trans@n.age.class)

      if (!is.array(mults)) {
        tmp.trans@waifw <- exper@trans@waifw * mults[t]
      } else {
        tmp.trans@waifw <- exper@trans@waifw * mults[,, t]
      }

      if (!is.na(exper@pop.rescale.each.timestep[t]))
        state <- exper@pop.rescale.each.timestep[t] * (state / sum(state))

      if (length(exper@births.per.1000.each.timestep) > 1) {
        tmp.trans@birth.rate <- exper@births.per.1000.each.timestep[t] * sum(state) / 1000
      } else {
        tmp.trans@birth.rate <- exper@trans@birth.rate * sum(rc[, t - 1]) / sum(exper@state.t0)
      }
      births.each.timestep[t] <- tmp.trans@birth.rate

      if (!is.na(exper@surv.each.timestep[1, 1])) {
        tmp.trans@age.surv.matrix <- ExtractAgeSpecificSurvivalMatrix(
          tmp.trans          = tmp.trans,
          maternal.obj       = exper@maternal.obj,
          surv.at.timestep.t = exper@surv.each.timestep[, t])
      }

      # --- Routine + scheduled SIA vaccination (identical to parent class) --

      if (!is.na(index.sia.vacc[t])) {
        r       <- routine$age.time.specific.routine[index.routine.vacc[t], ]
        s       <- SIA$age.time.specific.SIA[index.sia.vacc[t], ]
        rho.vec <- rep(0, exper@trans@n.age.class)
        MR1.age.range <- exper@time.specific.min.age.MR1[index.routine.vacc[t]]:
                         exper@time.specific.max.age.MR1[index.routine.vacc[t]]
        MR2.age.range <- exper@time.specific.min.age.MR2[index.routine.vacc[t]]:
                         exper@time.specific.max.age.MR2[index.routine.vacc[t]]
        rho.vec[MR1.age.range] <- exper@MR1SIAcorrelation
        rho.vec[MR2.age.range] <- exper@MR2SIAcorrelation
        tmp.trans@vac.per@pvacc.in.age.class <-
          r + s - r * s - rho.vec * sqrt(pmax(0, r * (1 - r) * s * (1 - s)))
        MR1.fail.each.timestep[t] <- routine$prop.fail.MR1[index.routine.vacc[t]]
        MR2.fail.each.timestep[t] <- routine$prop.fail.MR2[index.routine.vacc[t]]
        SIA.fail.each.timestep[t] <- SIA$prop.fail.SIA[index.sia.vacc[t]]
      } else {
        tmp.trans@vac.per@pvacc.in.age.class <-
          routine$age.time.specific.routine[index.routine.vacc[t], ]
        MR1.fail.each.timestep[t] <- routine$prop.fail.MR1[index.routine.vacc[t]]
        MR2.fail.each.timestep[t] <- routine$prop.fail.MR2[index.routine.vacc[t]]
      }

      # --- OBR: trigger check, response accumulation, and campaign delivery ---
      #
      # Step 1 (pre-state-update): check for new trigger; fire campaign if response delay elapsed.
      # Step 2 (post-state-update): accumulate case age distribution for pending response window.
      #
      # Trigger modes (or.trigger.mode):
      #   "I_scaled"  : metric = sum(I * reporting.rate) over surveillance window
      #   "confirmed" : metric = sum(confirmed.trigger[]) over surveillance window;
      #                 confirmed.trigger[] filled post-update (see below)
      #
      # Campaign age targeting:
      #   or.vacc.age.upper not NA -> use it directly
      #   or.vacc.age.upper NA     -> use CDF of pending.case.by.age at or.vacc.agedist.percentile
      #   or.response.case.type "suspected": CDF from I*r + non-measles background
      #   or.response.case.type "confirmed": CDF from estimated confirmed cases

      # 1a. Check for new trigger (only when no campaign is already pending)
      if (is.na(pending.trigger.t) &&
          t >= or.min.t &&
          (t - last.or.t) >= exper@or.min.interval) {

        t.end   <- t - exper@or.total.delay
        t.start <- t.end - exper@or.trigger.window + 1

        if (exper@or.trigger.mode == "confirmed") {
          metric <- sum(confirmed.trigger[t.start:t.end])
        } else {  # "I_scaled"
          metric <- sum(rc[trigger.i.inds, t.start:t.end] * or.rep.rate)
        }

        if (metric >= exper@or.n.confirmations.target) {
          pending.trigger.t <- t
        }
      }

      # 1b. Fire campaign when response delay has elapsed
      if (!is.na(pending.trigger.t) &&
          t == pending.trigger.t + exper@or.response.delay) {

        # Determine upper vaccination age bound
        if (!is.na(exper@or.vacc.age.upper)) {
          vacc.upper <- exper@or.vacc.age.upper
        } else if (sum(pending.case.by.age) > 0 &&
                   !is.na(exper@or.vacc.agedist.percentile)) {
          cdf        <- cumsum(pending.case.by.age) / sum(pending.case.by.age)
          pct.idx    <- which(cdf >= exper@or.vacc.agedist.percentile)[1]
          vacc.upper <- age.classes[pct.idx]
        } else {
          vacc.upper <- max(age.classes)
        }

        vacc.lower          <- exper@or.vacc.age.lower
        or.vacc.age.pos.dyn <- which(age.classes > vacc.lower &
                                     age.classes <= vacc.upper)
        or.vacc.prob.dyn    <- rep(0, exper@trans@n.age.class)
        or.vacc.prob.dyn[or.vacc.age.pos.dyn] <-
          exper@or.vacc.coverage *
          exper@obj.prob.vsucc@prob.vsucc[or.vacc.age.pos.dyn]

        p <- tmp.trans@vac.per@pvacc.in.age.class
        tmp.trans@vac.per@pvacc.in.age.class <- p + or.vacc.prob.dyn - p * or.vacc.prob.dyn
        or.times[t]         <- 1
        last.or.t           <- t
        pending.trigger.t   <- NA_integer_
        pending.case.by.age <- numeric(length(age.classes))
      }

      N0    <- sum(state)
      state <- next.ID.state(state, tmp.trans)
      NT    <- sum(state)
      growth.rate.each.timestep[t] <- log(NT / N0)

      rc[, t] <- state

      # Calendar month for this timestep
      month.t <- floor(((t - 1) %% steps.per.year.or) / steps.per.year.or * 12) + 1

      # Always compute true.reported.all (used by obs model and "suspected" accumulation)
      true.reported.all <- rc[exper@trans@i.inds, t] * or.rep.rate.all

      # Full observational model: 6-category age x timestep breakdown.
      # Runs every step when Se, Sp, and non_meas are provided.
      if (has.obs.model) {

        m.all          <- exper@or.non.meas.cases.by.age.month[, month.t]
        total.susp.all <- true.reported.all + m.all
        total.susp.sum <- sum(total.susp.all)

        if (total.susp.sum > 0) {
          pos.rate <- (sum(true.reported.all) / total.susp.sum) * exper@or.Se +
                      (sum(m.all)             / total.susp.sum) * (1 - exper@or.Sp)
          n.tests  <- if (pos.rate > 0)
            min(exper@or.n.confirmations.target / pos.rate, total.susp.sum)
          else 0

          prop.age         <- total.susp.all / total.susp.sum
          n.tested         <- n.tests * prop.age
          prop.true        <- ifelse(total.susp.all > 0,
                                     true.reported.all / total.susp.all, 0)
          n.true.tested    <- n.tested * prop.true
          n.nonmeas.tested <- n.tested * (1 - prop.true)

          obs.TP[, t]          <- n.true.tested    * exper@or.Se
          obs.FN_test[, t]     <- n.true.tested    * (1 - exper@or.Se)
          obs.FP_test[, t]     <- n.nonmeas.tested * (1 - exper@or.Sp)
          obs.TN[, t]          <- n.nonmeas.tested * exper@or.Sp
          obs.TP_clinical[, t] <- true.reported.all - n.true.tested
          obs.FP_clinical[, t] <- m.all             - n.nonmeas.tested
        }

      } else {
        m.all <- numeric(n.age.all)
      }

      # "confirmed" trigger mode: derive scalar from obs model (trigger age group only)
      if (exper@or.trigger.mode == "confirmed") {
        confirmed.trigger[t] <- sum(
          obs.TP[trigger.age.pos, t]          +
          obs.FP_test[trigger.age.pos, t]     +
          obs.TP_clinical[trigger.age.pos, t] +
          obs.FP_clinical[trigger.age.pos, t]
        )
      }

      # Accumulate case age distribution for the pending response window.
      # Uses already-computed obs model matrices — no duplicate calculation.
      if (!is.na(pending.trigger.t) &&
          t >= pending.trigger.t &&
          t < pending.trigger.t + exper@or.response.delay) {

        if (exper@or.trigger.mode == "I_scaled") {
          pending.case.by.age <- pending.case.by.age + true.reported.all + m.all
        } else {  # "confirmed"
          pending.case.by.age <- pending.case.by.age +
            obs.TP[, t]          +
            obs.FP_test[, t]     +
            obs.TP_clinical[, t] +
            obs.FP_clinical[, t]
        }
      }
    }

    result.rc <- new(
      "sim.results.MSIRV.update.demog.vaccine.change.outbreakresponse",
      data                      = rc,
      m.inds                    = exper@trans@m.inds,
      s.inds                    = exper@trans@s.inds,
      i.inds                    = exper@trans@i.inds,
      r.inds                    = exper@trans@r.inds,
      v.inds                    = exper@trans@v.inds,
      t                         = exper@t.min + (1:numTimeSteps - 1) * exper@step.size,
      age.class                 = exper@trans@age.class,
      births.each.timestep      = births.each.timestep,
      growth.rate.each.timestep = growth.rate.each.timestep,
      MR1.fail.each.timestep    = MR1.fail.each.timestep,
      MR2.fail.each.timestep    = MR2.fail.each.timestep,
      SIA.fail.each.timestep    = SIA.fail.each.timestep,
      routine.intro             = routine.intro,
      sia.times                 = sia.times,
      or.times        = or.times,
      obs.TP          = `rownames<-`(obs.TP,          age.classes),
      obs.FN_test     = `rownames<-`(obs.FN_test,     age.classes),
      obs.FP_test     = `rownames<-`(obs.FP_test,     age.classes),
      obs.TN          = `rownames<-`(obs.TN,          age.classes),
      obs.TP_clinical = `rownames<-`(obs.TP_clinical, age.classes),
      obs.FP_clinical = `rownames<-`(obs.FP_clinical, age.classes))

    rc <- new("experiment.result",
              experiment.def = exper,
              result         = result.rc)

    return(rc)
  }
)

#' @rdname run-methods
#' @aliases run,experiment.updatedemog.vaccinationchange.vaccinationlimitations,ANY-method
setMethod("run",
          "experiment.updatedemog.vaccinationchange.vaccinationlimitations",
          function(exper, rescale.WAIFW=TRUE, ...) {

            #print(rescale.WAIFW)
            state <- exper@state.t0

            #rescale the WAIFW if a specific R0 is specified
            if (rescale.WAIFW & length(exper@R0)>0) {
              #print("RESCALING!!!")
              exper@trans@waifw <- scaleWAIFW(exper@R0,
                                               state,exper@trans@waifw,
                                               frequency.dep=exper@trans@frequency.dep,
                                               suscept.state=exper@trans@s.inds[1])
            }

            #get the number of time steps in the experiment
            numTimeSteps <- round((exper@t.max-exper@t.min)/exper@step.size)+1

            if (!is.null(exper@season.obj)) {
              mults <- get.seasonal.mult(exper@t0.doy/365+(1:numTimeSteps-1)*exper@step.size,
                                         exper@season.obj)
            } else {
              mults <- rep(1,numTimeSteps)
            }

            #make a temporary transmission object
            tmp.trans <- exper@trans

            #hold the states as we walk through
            rc <- matrix(ncol = numTimeSteps, nrow = nrow(state))
            rc[,1] <- state

            #need output vector for births over time
            births.each.timestep <- growth.rate.each.timestep <- rep(NA, numTimeSteps)
            births.each.timestep[1] <- tmp.trans@birth.rate

            #there are a couple options that are not yet coded - stop experiment if these are selected
            if ((exper@MR1SIAcorrelation & exper@SIAinacc) | (exper@MR2SIAcorrelation & exper@SIAinacc))  {
              stop("error: SIA limitation can either be correlation with a routine dose OR inaccessible population (with or without inefficiency), not both")
            }
            if ((exper@MR1SIAcorrelation & exper@SIAinefficient) | (exper@MR2SIAcorrelation & exper@SIAinefficient))  {
              stop("error: SIA limitation can either be correlation with a routine dose OR SIA inefficiency population (with or without inaccessible), not both")
            }
            if (!exper@MR1SIAcorrelation & exper@MR2SIAcorrelation){
              stop("error: MR1SIAcorrelation = FALSE and MR2SIAcorrelation= T option isn't yet coded or available")
            }

            #generate the age and time specific vaccination matrices for routine and SIAs
            routine <- get.routine.time.age.specific(time.step= exper@step.size*12,
                                                     age.classes=exper@trans@age.class,
                                                     time.specific.MR1cov=exper@time.specific.MR1cov,
                                                     age.min.MR1=exper@time.specific.min.age.MR1,
                                                     age.max.MR1=exper@time.specific.max.age.MR1,
                                                     time.specific.MR2cov=exper@time.specific.MR2cov,
                                                     age.min.MR2=exper@time.specific.min.age.MR2,
                                                     age.max.MR2=exper@time.specific.max.age.MR2,
                                                     obj.vcdf.MR1=exper@obj.vcdf.MR1,
                                                     obj.vcdf.MR2=exper@obj.vcdf.MR2,
                                                     obj.prob.vsucc=exper@obj.prob.vsucc,
                                                     MR1MR2correlation=exper@MR1MR2correlation)
            routine.intro <- rep(0, numTimeSteps)
            if (any(exper@time.specific.MR1cov!=0)) routine.intro[min(which(exper@time.specific.MR1cov>0))*(1/exper@step.size)+1] <- 1
            if (any(exper@time.specific.MR2cov!=0)) routine.intro[min(which(exper@time.specific.MR2cov>0))*(1/exper@step.size)+1] <- 1
            index.routine.vacc <- c(1,rep(1:nrow(routine$age.time.specific.routine), each=(numTimeSteps-1)/exper@t.max))
            if (length(index.routine.vacc)<numTimeSteps) index.routine.vacc[(length(index.routine.vacc)+1):numTimeSteps] <- index.routine.vacc[length(index.routine.vacc)]

            SIA <- get.sia.time.age.specific(age.classes=exper@trans@age.class,
                                             time.specific.SIAcov=exper@time.specific.SIAcov,
                                             age.min.sia=exper@time.specific.min.age.SIA,
                                             age.max.sia=exper@time.specific.max.age.SIA,
                                             obj.prob.vsucc=exper@obj.prob.vsucc)
            index.sia.vacc <- rep(NA,numTimeSteps)
            year.sia <- which(exper@time.specific.SIAcov!=0)
            index.sia.vacc[(year.sia-1)*(numTimeSteps-1)/exper@t.max + round((exper@sia.timing.in.year*(numTimeSteps-1)/exper@t.max))] <-  year.sia #minus 1 because adding the sia.timing
            sia.times <- ifelse(!is.na(index.sia.vacc), 1, 0)

            #need output vectors for primary vaccination failure over time
            MR1.fail.each.timestep <- MR2.fail.each.timestep <- SIA.fail.each.timestep <- rep(0, numTimeSteps)
            MR1.fail.each.timestep[1] <- routine$prop.fail.MR1[1]
            MR2.fail.each.timestep[1] <- routine$prop.fail.MR2[1]
            SIA.fail.each.timestep[1] <- 0

            for (t in 2:numTimeSteps) {

              #introduction rates
              if (length(exper@intro.rate)>1) tmp.trans@introduction.rate <- rep(exper@intro.rate[index.routine.vacc[t]],exper@trans@n.age.class)

              #scale the waifw by seasonality
              if (!is.array(mults)){
                tmp.trans@waifw <- exper@trans@waifw*mults[t]
              } else {
                tmp.trans@waifw <- exper@trans@waifw*mults[,,t]
              }

              #if exper@pop.rescale.each.timestep[t]==NaN then tmp.trans@pop.rescale==NaN and will not rescale
              #otherwise if any number other than NaN it will rescale
              #or if exper@pop.rescale.each.timestep not completed , then numeric(0) and any index [t] is NA
              if (!is.na(exper@pop.rescale.each.timestep[t])) {state <- exper@pop.rescale.each.timestep[t] *(state/sum(state))}

              #put in the correct birth rate for that time-step, if it varies
              if (length(exper@births.per.1000.each.timestep)>1) {
                # Read 'birth rate'
                tmp.trans@birth.rate = (exper@births.per.1000.each.timestep[t]*sum(state)/1000)
              } else {
                tmp.trans@birth.rate <- exper@trans@birth.rate*sum(rc[,t-1])/sum(exper@state.t0)
              }
              births.each.timestep[t] <-  tmp.trans@birth.rate

              #put in time appropriate survival rate otherwise it uses original surv.matrix set up for trans object and keeps constant over time
              if (!is.na(exper@surv.each.timestep[1,1])) {
                tmp.trans@age.surv.matrix <- ExtractAgeSpecificSurvivalMatrix(tmp.trans=tmp.trans, maternal.obj=exper@maternal.obj,
                                                                              surv.at.timestep.t=exper@surv.each.timestep[,t])
              }

              ##put in correct vaccination coverage
              tmp.trans@vac.per@pvacc.in.age.class <- rep(0, exper@trans@n.age.class) #start with clean slate each time step
              routine.vacc.prob <- routine$age.time.specific.routine
              sia.vacc.prob <- SIA$age.time.specific.SIA

              if (!is.na(index.sia.vacc[t])){ #if SIA

                #IF SIA independent MR1 or MR2
                if (!exper@MR1SIAcorrelation & !exper@MR2SIAcorrelation & !exper@SIAinacc & !exper@SIAinefficient){
                  tmp.trans@vac.per@pvacc.in.age.class <-
                    routine.vacc.prob[index.routine.vacc[t],] +
                    sia.vacc.prob[index.sia.vacc[t],] - #xxamy changed from addition to subtraction
                    (routine.vacc.prob[index.routine.vacc[t],]*sia.vacc.prob[index.sia.vacc[t],])
                }

                #IF SIA and MR1 correlated, SIA & MR2 NOT correlated
                if (exper@MR1SIAcorrelation & !exper@MR2SIAcorrelation & !exper@SIAinacc & !exper@SIAinefficient) {
                  MR1.age.range <- exper@time.specific.min.age.MR1[index.routine.vacc[t]]:exper@time.specific.max.age.MR1[index.routine.vacc[t]]
                  MR2.age.range <- exper@time.specific.min.age.MR2[index.routine.vacc[t]]:exper@time.specific.max.age.MR2[index.routine.vacc[t]]
                  if (exper@time.specific.SIAcov[index.sia.vacc[t]]<exper@time.specific.MR1cov[index.routine.vacc[t]]){
                    #non MR1 or MR2 ages
                    tmp.trans@vac.per@pvacc.in.age.class <- sia.vacc.prob[index.sia.vacc[t],]
                    #MR1 ages
                    tmp.trans@vac.per@pvacc.in.age.class[MR1.age.range] <-
                      routine.vacc.prob[index.routine.vacc[t],MR1.age.range] + routine$one.minus.ve1[index.routine.vacc[t]]*sia.vacc.prob[index.sia.vacc[t],MR1.age.range]
                    #MR2 ages
                    tmp.trans@vac.per@pvacc.in.age.class[MR2.age.range] <-
                      routine.vacc.prob[index.routine.vacc[t],MR2.age.range] +
                      sia.vacc.prob[index.sia.vacc[t],MR2.age.range] +
                      (routine.vacc.prob[index.routine.vacc[t],MR2.age.range]*sia.vacc.prob[index.sia.vacc[t],MR2.age.range])
                  }
                  if (exper@time.specific.SIAcov[index.sia.vacc[t]]>exper@time.specific.MR1cov[index.routine.vacc[t]]){
                    #non MR1 or MR2 ages
                    tmp.trans@vac.per@pvacc.in.age.class <- sia.vacc.prob[index.sia.vacc[t],]
                    #MR1 ages
                    tmp.trans@vac.per@pvacc.in.age.class[MR1.age.range] <-
                      routine.vacc.prob[index.routine.vacc[t],MR1.age.range] + routine$prop.fail.MR1.byage[index.routine.vacc[t],MR1.age.range] +
                      (sia.vacc.prob[index.sia.vacc[t],MR1.age.range] - routine.vacc.prob[index.routine.vacc[t],MR1.age.range])
                    #MR2 ages
                    tmp.trans@vac.per@pvacc.in.age.class[MR2.age.range] <-
                      routine.vacc.prob[index.routine.vacc[t],MR2.age.range] +
                      sia.vacc.prob[index.sia.vacc[t],MR2.age.range] +
                      (routine.vacc.prob[index.routine.vacc[t],MR2.age.range]*sia.vacc.prob[index.sia.vacc[t],MR2.age.range])
                  }
                }

                #IF SIA and MR1 & SIA and MR2 correlated
                if (exper@MR1SIAcorrelation & exper@MR2SIAcorrelation & !exper@SIAinacc & !exper@SIAinefficient) {
                  MR1.age.range <- exper@time.specific.min.age.MR1[index.routine.vacc[t]]:exper@time.specific.max.age.MR1[index.routine.vacc[t]]
                  MR2.age.range <- exper@time.specific.min.age.MR2[index.routine.vacc[t]]:exper@time.specific.max.age.MR2[index.routine.vacc[t]]
                  if (exper@time.specific.SIAcov[index.sia.vacc[t]]<exper@time.specific.MR1cov[index.routine.vacc[t]]){
                    #non MR1 or MR2 ages
                    tmp.trans@vac.per@pvacc.in.age.class <- sia.vacc.prob[index.sia.vacc[t],]
                    #MR1 ages
                    tmp.trans@vac.per@pvacc.in.age.class[MR1.age.range] <-
                      routine.vacc.prob[index.routine.vacc[t],MR1.age.range] + routine$one.minus.ve1[index.routine.vacc[t]]*sia.vacc.prob[index.sia.vacc[t],MR1.age.range]
                    #MR2 ages
                    tmp.trans@vac.per@pvacc.in.age.class[MR2.age.range] <-
                      routine.vacc.prob[index.routine.vacc[t],MR2.age.range] + routine$one.minus.ve1[index.routine.vacc[t]]*sia.vacc.prob[index.sia.vacc[t],MR2.age.range]
                  }
                  if (exper@time.specific.SIAcov[index.sia.vacc[t]]>exper@time.specific.MR1cov[index.routine.vacc[t]]){
                    #non MR1 or MR2 ages
                    tmp.trans@vac.per@pvacc.in.age.class <- sia.vacc.prob[index.sia.vacc[t],]
                    #MR1 ages
                    tmp.trans@vac.per@pvacc.in.age.class[MR1.age.range] <-
                      routine.vacc.prob[index.routine.vacc[t],MR1.age.range] + routine$prop.fail.MR1.byage[index.routine.vacc[t],MR1.age.range] +
                      (sia.vacc.prob[index.sia.vacc[t],MR1.age.range] - routine.vacc.prob[index.routine.vacc[t],MR1.age.range])
                    #MR2 ages
                    tmp.trans@vac.per@pvacc.in.age.class[MR2.age.range] <-
                      routine.vacc.prob[index.routine.vacc[t],MR2.age.range] + routine$prop.fail.MR1.byage[index.routine.vacc[t],MR2.age.range] +
                      (sia.vacc.prob[index.sia.vacc[t],MR2.age.range] - routine.vacc.prob[index.routine.vacc[t],MR2.age.range])
                  }
                }

                #If inaccessible population, but not SIA inefficiency
                if (!exper@MR1SIAcorrelation & !exper@MR2SIAcorrelation & exper@SIAinacc & !exper@SIAinefficient) {
                  #assumes coverage and accessibility independent processes - if coverage is 90% and prop.inacc = 0.1, effective population-level coverage of 0.81
                  tmp.trans@vac.per@pvacc.in.age.class <- 1-(exper@prop.inacc[t]+
                                                               (1-exper@prop.inacc[t])*(1-(routine.vacc.prob[index.routine.vacc[t],] +
                                                                                             sia.vacc.prob[index.sia.vacc[t],] -
                                                                                             (routine.vacc.prob[index.routine.vacc[t],]*sia.vacc.prob[index.sia.vacc[t],]))))
                  #tmp.trans@vac.per@pvacc.in.age.class <- 1-(exper@prop.inacc[t]+(1-exper@prop.inacc[t])*(1-routine.vacc.prob[index.routine.vacc[t],])*(1-sia.vacc.prob[index.sia.vacc[t],]))
                }

                #If inaccessible population AND SIA inefficiency
                if (!exper@MR1SIAcorrelation & !exper@MR2SIAcorrelation & exper@SIAinacc & exper@SIAinefficient) {
                  # with SIA inaccessible & with SIA inefficiency -- VIMC 2017-2021 versions - xxamy
                  z <- routine.vacc.prob[index.routine.vacc[t],] #prob. successful routine given access
                  ro <- (1-exper@prop.inacc[t]) #prob. accessible
                  m <- (1-(sia.vacc.prob[index.sia.vacc[t],]*(1/ro)*(1-0.1)))^(1/(1-0.1)) #prob. not successful campaign vaccination given access, and inefficiency of 0.1
                  m[is.na(m)] <- 0 #if NaN then negative number was raised, which means everyone in accessible population was vaccinated by campaign, therefore prop. not vaccinated =0
                  age.spec.vacc.prob <- 1-((1-ro)+(ro*m)) #ages over 36 months, where only eligible for campaign
                  max.age.routine <- exper@time.specific.max.age.MR2[index.routine.vacc[t]]
                  age.spec.vacc.prob[1:max.age.routine] <- 1-((1-ro)+(ro*m*(1-(z/ro))))[1:max.age.routine] #ages 1:36 months where routine takes place
                  #if(sum(age.spec.vacc.prob)>0) print(age.spec.vacc.prob[1:max.age.routine])
                  tmp.trans@vac.per@pvacc.in.age.class <- age.spec.vacc.prob
                }

                #stow output
                MR1.fail.each.timestep[t] <- routine$prop.fail.MR1[index.routine.vacc[t]]
                MR2.fail.each.timestep[t] <- routine$prop.fail.MR2[index.routine.vacc[t]]
                SIA.fail.each.timestep[t] <- SIA$prop.fail.SIA[index.sia.vacc[t]]

              } else { #if no SIA
                tmp.trans@vac.per@pvacc.in.age.class <- routine.vacc.prob[index.routine.vacc[t],]
                #stow output
                MR1.fail.each.timestep[t] <- routine$prop.fail.MR1[index.routine.vacc[t]]
                MR2.fail.each.timestep[t] <- routine$prop.fail.MR2[index.routine.vacc[t]]
              }

              N0 <- sum(state)
              state <- next.ID.state(state, tmp.trans)
              NT <- sum(state)
              growth.rate.each.timestep[t] <- log(NT/N0) #instantaneous biweekly growth rate

              #print(t)
              #print(dim(state))
              rc [,t] <- state

            }


            rc <- new("sim.results.MSIRV.update.demog.vaccine.change",
                      data=rc,
                      m.inds = exper@trans@m.inds,
                      s.inds = exper@trans@s.inds,
                      i.inds = exper@trans@i.inds,
                      r.inds = exper@trans@r.inds,
                      v.inds = exper@trans@v.inds,
                      t = exper@t.min+(1:numTimeSteps-1)* exper@step.size,
                      age.class = exper@trans@age.class,
                      births.each.timestep = births.each.timestep,
                      growth.rate.each.timestep = growth.rate.each.timestep,
                      MR1.fail.each.timestep = MR1.fail.each.timestep,
                      MR2.fail.each.timestep = MR2.fail.each.timestep,
                      SIA.fail.each.timestep = SIA.fail.each.timestep,
                      routine.intro = routine.intro,
                      sia.times = sia.times)

            rc <- new("experiment.result",
                      experiment.def = exper,
                      result = rc)


            return(rc)
          }
)


# Run method for experiment.updatedemog.vaccinationchange.spatial objects

#' @rdname run-methods
#' @aliases run,experiment.updatedemog.vaccinationchange.spatial,ANY-method
setMethod("run",
          "experiment.updatedemog.vaccinationchange.spatial",
          function(exper, rescale.WAIFW=T, ...) {

            #print(rescale.WAIFW)
            state <- exper@state.t0

            #number of sub-populations
            n.subpops <- exper@trans@n.subpops

            #rescale the WAIFW if a specific R0 is specified
            if (rescale.WAIFW & length(exper@R0)>0) {
              #print("RESCALING!!!")
              exper@trans@waifw <- scaleWAIFW.space(exper@R0,
                                                     state,exper@trans@waifw,
                                                     frequency.dep=exper@trans@frequency.dep,
                                                     suscept.state=exper@trans@s.inds[1])
            }

            #get the number of time steps in the experiment
            T <- round((exper@t.max-exper@t.min)/exper@step.size)+1

            #get seasonal mult
            if (!is.null(exper@season.obj)) {
              mults <- get.seasonal.mult(exper@t0.doy/365+(1:T-1)*exper@step.size,
                                         exper@season.obj)
            } else {
              mults <- rep(1,T)
            }

            #make a temporary transmission object
            tmp.trans <- exper@trans

            #hold the states as we walk through
            rc <- matrix(ncol = T, nrow = nrow(state))
            rc[,1] <- state

            #need output vector for births over time
            births.each.timestep <- growth.rate.each.timestep <- N0 <- matrix(NA, n.subpops, T)
            births.each.timestep[,1] <- tmp.trans@birth.rate

            #generate the age and space (rows) and year (columns) specific vaccination matrix
            routine <- space.wrapper.get.routine.time.age.specific(time.step= exper@step.size*12,
                                                                   age.classes=exper@trans@age.class,
                                                                   space.time.specific.MR1cov=exper@time.specific.MR1cov,
                                                                   age.min.MR1=exper@time.specific.min.age.MR1,
                                                                   age.max.MR1=exper@time.specific.max.age.MR1,
                                                                   space.time.specific.MR2cov=exper@time.specific.MR2cov,
                                                                   age.min.MR2=exper@time.specific.min.age.MR2,
                                                                   age.max.MR2=exper@time.specific.max.age.MR2,
                                                                   obj.vcdf.MR1=exper@obj.vcdf.MR1,
                                                                   obj.vcdf.MR2=exper@obj.vcdf.MR2,
                                                                   obj.prob.vsucc=exper@obj.prob.vsucc,
                                                                   MR1MR2correlation=F)

            #getting year of routine introductions
            routine.intro <- rep(0, T)
            if (any(colSums(exper@time.specific.MR1cov)!=0)) routine.intro[min(which(colSums(exper@time.specific.MR1cov)>0))*(1/exper@step.size)+1] <- 1
            if (any(colSums(exper@time.specific.MR2cov)!=0)) routine.intro[min(which(colSums(exper@time.specific.MR2cov)>0))*(1/exper@step.size)+1] <- 1

            #specifying that each year (column) of routine should be repeated 24 times
            index.routine.vacc <- c(1,rep(1:ncol(routine$age.time.specific.routine), each=(T-1)/exper@t.max))
            if (length(index.routine.vacc)<T) index.routine.vacc[(length(index.routine.vacc)+1):T] <- index.routine.vacc[length(index.routine.vacc)]

            #generate the age and space (rows) and year (columns) specific vaccination matrix - same size as `routine`
            SIA <- space.wrapper.get.sia.time.age.specific(age.classes=exper@trans@age.class,
                                                           space.time.specific.SIAcov=exper@time.specific.SIAcov,
                                                           age.min.sia=exper@time.specific.min.age.SIA,
                                                           age.max.sia=exper@time.specific.max.age.SIA,
                                                           obj.prob.vsucc=exper@obj.prob.vsucc)

            #getting sia.times as a vector of 0 / 1 to represent when the SIA is to take place
            index.sia.vacc <- rep(NA,T)
            year.sia <- which(colSums(exper@time.specific.SIAcov)!=0)
            index.sia.vacc[(year.sia-1)*(T-1)/exper@t.max + round((exper@sia.timing.in.year*(T-1)/exper@t.max))] <-  year.sia #minus 1 because adding the sia.timing
            sia.times <- ifelse(!is.na(index.sia.vacc), 1, 0)

            #need output vectors for primary vaccination failure over time
            MR1.fail.each.timestep <- MR2.fail.each.timestep <- SIA.fail.each.timestep <- matrix(0, n.subpops, T)
            MR1.fail.each.timestep[,1] <- routine$prop.fail.MR1[,1]
            MR2.fail.each.timestep[,1] <- routine$prop.fail.MR2[,1]
            SIA.fail.each.timestep[,1] <- 0

            for (t in 2:T) {

              print(t)

              #scale the waifw by seasonality
              if (!is.array(mults)){
                tmp.trans@waifw <- exper@trans@waifw*mults[t]
              } else {
                tmp.trans@waifw <- exper@trans@waifw*mults[,,t]
              }

              #if exper@pop.rescale.each.timestep[t]==NaN then tmp.trans@pop.rescale==NaN and will not rescale
              #otherwise if any number other than NaN it will rescale
              #or if exper@pop.rescale.each.timestep not completed , then numeric(0) and any index [t] is NA
              for (s in 1:n.subpops) {
                if (!is.na(exper@pop.rescale.each.timestep[s,t])) {
                  last.index.tmp <- (s*exper@trans@n.epi.class*exper@trans@n.age.class)
                  first.index.tmp <- last.index.tmp-(exper@trans@n.epi.class*exper@trans@n.age.class)+1
                  state[first.index.tmp:last.index.tmp] <- exper@pop.rescale.each.timestep[s,t]*
                    (state[first.index.tmp:last.index.tmp]/sum(state[first.index.tmp:last.index.tmp]))
                }
              }

              #put in the correct birth rate for that time-step, if it varies
              for (s in 1:n.subpops) {
                if (ncol(exper@births.per.1000.each.timestep)>1) {
                  last.index.tmp <- (s*exper@trans@n.epi.class*exper@trans@n.age.class)
                  first.index.tmp <- last.index.tmp-(exper@trans@n.epi.class*exper@trans@n.age.class)+1
                  tmp.trans@birth.rate[s] = (exper@births.per.1000.each.timestep[s,t]*
                                               sum(state[first.index.tmp:last.index.tmp])/1000)
                } else {
                  tmp.trans@birth.rate[s] <- exper@trans@birth.rate[s]
                }
              }
              births.each.timestep[,t] <-  tmp.trans@birth.rate

              #put in time appropriate survival rate otherwise it uses original surv.matrix set up for trans object and keeps constant over time
              if (!is.na(exper@surv.each.timestep[1,1])) {
                for (s in 1:n.subpops) {
                  last.index.tmp1 <- (s*exper@trans@n.age.class)
                  first.index.tmp1 <- last.index.tmp1-(exper@trans@n.age.class)+1
                  tmp.age.surv.matrix <- ExtractAgeSpecificSurvivalMatrix(tmp.trans=tmp.trans, maternal.obj=exper@maternal.obj,
                                                                          surv.at.timestep.t=exper@surv.each.timestep[first.index.tmp1:last.index.tmp1,t])

                  last.index.tmp2 <- (s*exper@trans@n.epi.class*exper@trans@n.age.class)
                  first.index.tmp2 <- last.index.tmp2-(exper@trans@n.epi.class*exper@trans@n.age.class)+1
                  tmp.trans@age.surv.matrix[first.index.tmp2:last.index.tmp2,] <- tmp.age.surv.matrix
                }
              }

              ##put in correct vaccination coverage
              for (s in 1:n.subpops) {
                last.index.tmp <- (s*exper@trans@n.age.class)
                first.index.tmp <- last.index.tmp-(exper@trans@n.age.class)+1

                routine.vacc.prob <- routine$age.time.specific.routine[first.index.tmp:last.index.tmp,]
                sia.vacc.prob <- SIA$age.time.specific.SIA[first.index.tmp:last.index.tmp,]

                if (!is.na(index.sia.vacc[t])){ #if SIA
                  tmp.trans@vac.per@pvacc.in.age.class[first.index.tmp:last.index.tmp] <-
                    routine.vacc.prob[,index.routine.vacc[t]] +
                    sia.vacc.prob[,index.sia.vacc[t]] -
                    (routine.vacc.prob[,index.routine.vacc[t]]*sia.vacc.prob[,index.sia.vacc[t]])
                  #stow output
                  MR1.fail.each.timestep[t] <- routine$prop.fail.MR1[s,index.routine.vacc[t]]
                  MR2.fail.each.timestep[t] <- routine$prop.fail.MR2[s,index.routine.vacc[t]]
                  SIA.fail.each.timestep[t] <- SIA$prop.fail.SIA[s,index.sia.vacc[t]]

                } else { #if no SIA
                  tmp.trans@vac.per@pvacc.in.age.class[first.index.tmp:last.index.tmp] <-
                    routine.vacc.prob[,index.routine.vacc[t]]
                  #stow output
                  MR1.fail.each.timestep[t] <- routine$prop.fail.MR1[s,index.routine.vacc[t]]
                  MR2.fail.each.timestep[t] <- routine$prop.fail.MR2[s,index.routine.vacc[t]]
                }
              }

              #stow the previous population size
              for (s in 1:n.subpops) {
                last.index.tmp <- (s*exper@trans@n.epi.class*exper@trans@n.age.class)
                first.index.tmp <- last.index.tmp-(exper@trans@n.epi.class*exper@trans@n.age.class)+1
                N0[s,t] <- sum(state[first.index.tmp:last.index.tmp])
              }

              #run experiment to get the next time step
              state <- next.ID.state(state, tran=tmp.trans)

              #stow the previous population size
              for (s in 1:n.subpops) {
                last.index.tmp <- (s*exper@trans@n.epi.class*exper@trans@n.age.class)
                first.index.tmp <- last.index.tmp-(exper@trans@n.epi.class*exper@trans@n.age.class)+1
                NT <- sum(state[first.index.tmp:last.index.tmp])
                growth.rate.each.timestep[s,t] <- log(NT/N0[s,t]) #instantaneous biweekly growth rate
              }

              #print(t)
              #print(dim(state))
              rc [,t] <- state
            }


            result <- new("sim.results.MSIRV.update.demog.vaccine.change.space",
                          data=rc,
                          m.inds = exper@trans@m.inds,
                          s.inds = exper@trans@s.inds,
                          i.inds = exper@trans@i.inds,
                          r.inds = exper@trans@r.inds,
                          v.inds = exper@trans@v.inds,
                          t = exper@t.min+(1:T-1)* exper@step.size,
                          age.class = exper@trans@age.class,
                          births.each.timestep = births.each.timestep,
                          growth.rate.each.timestep = growth.rate.each.timestep,
                          MR1.fail.each.timestep = MR1.fail.each.timestep,
                          MR2.fail.each.timestep = MR2.fail.each.timestep,
                          SIA.fail.each.timestep = SIA.fail.each.timestep,
                          routine.intro = routine.intro,
                          sia.times = sia.times,
                          n.subpops = exper@trans@n.subpops,
                          subpop.class.label = exper@trans@subpop.class.label)

            rc <- new("experiment.result",
                      experiment.def = exper,
                      result = result)


            return(rc)
          }
)


# Run method for experiment.updatedemog.spatial objects

#' @rdname run-methods
#' @aliases run,experiment.updatedemog.spatial,ANY-method
setMethod("run",
          "experiment.updatedemog.spatial",
          function(exper, rescale.WAIFW=T, ...){

            #print(rescale.WAIFW)
            state <- exper@state.t0

            #number of sub-populations
            n.subpops <- exper@trans@n.subpops

            #rescale the WAIFW if a specific R0 is specified
            if (rescale.WAIFW & length(exper@R0)>0) {
              #print("RESCALING!!!")
              exper@trans@waifw <- scaleWAIFW.space(exper@R0,
                                                     state,exper@trans@waifw,
                                                     frequency.dep=exper@trans@frequency.dep,
                                                     suscept.state=exper@trans@s.inds[1])
            }

            #get the number of time steps in the experiment
            T <- round((exper@t.max-exper@t.min)/exper@step.size)+1

            #get seasonal mult
            if (!is.null(exper@season.obj)) {
              mults <- get.seasonal.mult(exper@t0.doy/365+(1:T-1)*exper@step.size,
                                         exper@season.obj)
            } else {
              mults <- rep(1,T)
            }

            #make a temporary transmission object
            tmp.trans <- exper@trans

            #hold the states as we walk through
            rc <- matrix(ncol = T, nrow = nrow(state))
            rc[,1] <- state

            #need output vector for births over time
            births.each.timestep <- growth.rate.each.timestep <- N0 <- matrix(NA, n.subpops, T)
            births.each.timestep[,1] <- tmp.trans@birth.rate

            #output vector for when SIAs were administered - default in this experiment is 0
            sia.times <- routine.intro <- rep(0,T)


            for (t in 2:T) {

              print(t)

              #scale the waifw by seasonality
              if (!is.array(mults)){
                tmp.trans@waifw <- exper@trans@waifw*mults[t]
              } else {
                tmp.trans@waifw <- exper@trans@waifw*mults[,,t]
              }

              #if exper@pop.rescale.each.timestep[t]==NaN then tmp.trans@pop.rescale==NaN and will not rescale
              #otherwise if any number other than NaN it will rescale
              #or if exper@pop.rescale.each.timestep not completed , then numeric(0) and any index [t] is NA
              for (s in 1:n.subpops) {
                if (!is.na(exper@pop.rescale.each.timestep[s,t])) {
                  last.index.tmp <- (s*exper@trans@n.epi.class*exper@trans@n.age.class)
                  first.index.tmp <- last.index.tmp-(exper@trans@n.epi.class*exper@trans@n.age.class)+1
                  state[first.index.tmp:last.index.tmp] <- exper@pop.rescale.each.timestep[s,t]*
                    (state[first.index.tmp:last.index.tmp]/sum(state[first.index.tmp:last.index.tmp]))
                }
              }

              #put in the correct birth rate for that time-step, if it varies
              for (s in 1:n.subpops) {
                if (ncol(exper@births.per.1000.each.timestep)>1) {
                  last.index.tmp <- (s*exper@trans@n.epi.class*exper@trans@n.age.class)
                  first.index.tmp <- last.index.tmp-(exper@trans@n.epi.class*exper@trans@n.age.class)+1
                  tmp.trans@birth.rate[s] = (exper@births.per.1000.each.timestep[s,t]*
                                               sum(state[first.index.tmp:last.index.tmp])/1000)
                } else {
                  tmp.trans@birth.rate[s] <- exper@trans@birth.rate[s]
                }
              }
              births.each.timestep[,t] <-  tmp.trans@birth.rate

              #put in time appropriate survival rate otherwise it uses original surv.matrix set up for trans object and keeps constant over time
              if (!is.na(exper@surv.each.timestep[1,1])) {
                for (s in 1:n.subpops) {
                  last.index.tmp1 <- (s*exper@trans@n.age.class)
                  first.index.tmp1 <- last.index.tmp1-(exper@trans@n.age.class)+1
                  tmp.age.surv.matrix <- ExtractAgeSpecificSurvivalMatrix(tmp.trans=tmp.trans, maternal.obj=exper@maternal.obj,
                                                                          surv.at.timestep.t=exper@surv.each.timestep[first.index.tmp1:last.index.tmp1,t])

                  last.index.tmp2 <- (s*exper@trans@n.epi.class*exper@trans@n.age.class)
                  first.index.tmp2 <- last.index.tmp2-(exper@trans@n.epi.class*exper@trans@n.age.class)+1
                  tmp.trans@age.surv.matrix[first.index.tmp2:last.index.tmp2,] <- tmp.age.surv.matrix
                }
              }

              #stow the previous population size
              for (s in 1:n.subpops) {
                last.index.tmp <- (s*exper@trans@n.epi.class*exper@trans@n.age.class)
                first.index.tmp <- last.index.tmp-(exper@trans@n.epi.class*exper@trans@n.age.class)+1
                N0[s,t] <- sum(state[first.index.tmp:last.index.tmp])
              }

              #run experiment to get the next time step
              state <- next.ID.state(state, tran=tmp.trans)

              #stow the previous population size
              for (s in 1:n.subpops) {
                last.index.tmp <- (s*exper@trans@n.epi.class*exper@trans@n.age.class)
                first.index.tmp <- last.index.tmp-(exper@trans@n.epi.class*exper@trans@n.age.class)+1
                NT <- sum(state[first.index.tmp:last.index.tmp])
                growth.rate.each.timestep[s,t] <- log(NT/N0[s,t]) #instantaneous biweekly growth rate
              }

              #print(t)
              #print(dim(state))
              rc [,t] <- state
            }

            result <- new("sim.results.MSIRV.update.demog.space",
                          data=rc,
                          m.inds = exper@trans@m.inds,
                          s.inds = exper@trans@s.inds,
                          i.inds = exper@trans@i.inds,
                          r.inds = exper@trans@r.inds,
                          v.inds = exper@trans@v.inds,
                          t = exper@t.min+(1:T-1)* exper@step.size,
                          age.class = exper@trans@age.class,
                          births.each.timestep = births.each.timestep,
                          growth.rate.each.timestep = growth.rate.each.timestep,
                          routine.intro = routine.intro,
                          sia.times = sia.times,
                          n.subpops = exper@trans@n.subpops,
                          subpop.class.label = exper@trans@subpop.class.label)

            rc <- new("experiment.result",
                      experiment.def = exper,
                      result = result)


            return(rc)
          }
)


# Run method for experiment.updatedemog.vaccinationchange.vaccinationlimitations.spatial objects

#' @rdname run-methods
#' @aliases run,experiment.updatedemog.vaccinationchange.vaccinationlimitations.spatial,ANY-method
setMethod("run",
          "experiment.updatedemog.vaccinationchange.vaccinationlimitations.spatial",
          function(exper, rescale.WAIFW=T, ...){
            print("THIS EXPERIMENT RUN METHOD HAS NOT YET BEEN CODED")
          }
)


# Run method for experiment.updatedemog.vaccinationchange.spatial.schoolvaccination objects

#' @rdname run-methods
#' @aliases run,experiment.updatedemog.vaccinationchange.spatial,schoolvaccination,ANY-method
setMethod("run",
          "experiment.updatedemog.vaccinationchange.school.spatial",
          function(exper, rescale.WAIFW=T, ...) {

            #print(rescale.WAIFW)
            state <- exper@state.t0

            #number of sub-populations
            n.subpops <- exper@trans@n.subpops

            #rescale the WAIFW if a specific R0 is specified
            if (rescale.WAIFW & length(exper@R0)>0) {
              #print("RESCALING!!!")
              exper@trans@waifw <- scaleWAIFW.space(exper@R0,
                                                    state,exper@trans@waifw,
                                                    frequency.dep=exper@trans@frequency.dep,
                                                    suscept.state=exper@trans@s.inds[1])
            }

            #get the number of time steps in the experiment
            T <- round((exper@t.max-exper@t.min)/exper@step.size)+1

            #get seasonal mult
            if (!is.null(exper@season.obj)) {
              mults <- get.seasonal.mult(exper@t0.doy/365+(1:T-1)*exper@step.size,
                                         exper@season.obj)
            } else {
              mults <- rep(1,T)
            }

            #make a temporary transmission object
            tmp.trans <- exper@trans

            #hold the states as we walk through
            rc <- matrix(ncol = T, nrow = nrow(state))
            rc[,1] <- state

            #need output vector for births over time
            births.each.timestep <- growth.rate.each.timestep <- N0 <- matrix(NA, n.subpops, T)
            births.each.timestep[,1] <- tmp.trans@birth.rate

            #generate the age and space (rows) and year (columns) specific vaccination matrix
            routine <- space.wrapper.get.routine.time.age.specific(time.step= exper@step.size*12,
                                                                   age.classes=exper@trans@age.class,
                                                                   space.time.specific.MR1cov=exper@time.specific.MR1cov,
                                                                   age.min.MR1=exper@time.specific.min.age.MR1,
                                                                   age.max.MR1=exper@time.specific.max.age.MR1,
                                                                   space.time.specific.MR2cov=exper@time.specific.MR2cov,
                                                                   age.min.MR2=exper@time.specific.min.age.MR2,
                                                                   age.max.MR2=exper@time.specific.max.age.MR2,
                                                                   obj.vcdf.MR1=exper@obj.vcdf.MR1,
                                                                   obj.vcdf.MR2=exper@obj.vcdf.MR2,
                                                                   obj.prob.vsucc=exper@obj.prob.vsucc,
                                                                   MR1MR2correlation=F)

            #getting year of routine introductions
            routine.intro <- rep(0, T)
            if (any(colSums(exper@time.specific.MR1cov)!=0)) routine.intro[min(which(colSums(exper@time.specific.MR1cov)>0))*(1/exper@step.size)+1] <- 1
            if (any(colSums(exper@time.specific.MR2cov)!=0)) routine.intro[min(which(colSums(exper@time.specific.MR2cov)>0))*(1/exper@step.size)+1] <- 1

            #specifying that each year (column) of routine should be repeated 24 times
            index.routine.vacc <- c(1,rep(1:ncol(routine$age.time.specific.routine), each=(T-1)/exper@t.max))
            if (length(index.routine.vacc)<T) index.routine.vacc[(length(index.routine.vacc)+1):T] <- index.routine.vacc[length(index.routine.vacc)]

            #generate the age and space (rows) and year (columns) specific vaccination matrix - same size as `routine`
            SIA <- space.wrapper.get.sia.time.age.specific(age.classes=exper@trans@age.class,
                                                           space.time.specific.SIAcov=exper@time.specific.SIAcov,
                                                           age.min.sia=exper@time.specific.min.age.SIA,
                                                           age.max.sia=exper@time.specific.max.age.SIA,
                                                           obj.prob.vsucc=exper@obj.prob.vsucc)

            #getting sia.times as a vector of 0 / 1 to represent when the SIA is to take place
            index.sia.vacc <- rep(NA,T)
            year.sia <- which(colSums(exper@time.specific.SIAcov)!=0)
            index.sia.vacc[(year.sia-1)*(T-1)/exper@t.max + round((exper@sia.timing.in.year*(T-1)/exper@t.max))] <-  year.sia #minus 1 because adding the sia.timing
            sia.times <- ifelse(!is.na(index.sia.vacc), 1, 0)

            #generate the age and space (rows) and year (columns) specific vaccination matrix
            schoolvacc <- space.wrapper.get.schoolvacc.time.age.specific(time.step= exper@step.size*12,
                                                                                         age.classes=exper@trans@age.class,
                                                                                         space.time.specific.cov=exper@time.specific.schoolvacc.cov,
                                                                                         age.min=exper@time.specific.min.age.schoolvacc,
                                                                                         age.max=exper@time.specific.max.age.schoolvacc,
                                                                                         list.obj.vcdf=exper@list.obj.vcdf.schoolvacc,
                                                                                         obj.prob.vsucc=exper@obj.prob.vsucc.schoolvacc)

            #getting school.vaccination.times as a vector of 0 / 1 to represent when the SIA is to take place
            index.school.vacc <- rep(NA,T)
            year.schoolvacc <- which(colSums(exper@time.specific.schoolvacc.cov)!=0)
            year.schoolvacc[(year.schoolvacc-1)*(T-1)/exper@t.max + round((exper@schoolvacc.timing.in.year*(T-1)/exper@t.max))] <-  year.schoolvacc #minus 1 because adding the sia.timing
            schoolvacc.times <- ifelse(!is.na(year.schoolvacc), 1, 0)


            #need output vectors for primary vaccination failure over time
            MR1.fail.each.timestep <- MR2.fail.each.timestep <- SIA.fail.each.timestep <- schoolvacc.fail.each.timestep <- matrix(0, n.subpops, T)
            MR1.fail.each.timestep[,1] <- routine$prop.fail.MR1[,1]
            MR2.fail.each.timestep[,1] <- routine$prop.fail.MR2[,1]
            SIA.fail.each.timestep[,1] <- 0
            schoolvacc.fail.each.timestep[,1] <- 0

            for (t in 2:T) {

              print(t)

              #scale the waifw by seasonality
              if (!is.array(mults)){
                tmp.trans@waifw <- exper@trans@waifw*mults[t]
              } else {
                tmp.trans@waifw <- exper@trans@waifw*mults[,,t]
              }

              #if exper@pop.rescale.each.timestep[t]==NaN then tmp.trans@pop.rescale==NaN and will not rescale
              #otherwise if any number other than NaN it will rescale
              #or if exper@pop.rescale.each.timestep not completed , then numeric(0) and any index [t] is NA
              for (s in 1:n.subpops) {
                if (!is.na(exper@pop.rescale.each.timestep[s,t])) {
                  last.index.tmp <- (s*exper@trans@n.epi.class*exper@trans@n.age.class)
                  first.index.tmp <- last.index.tmp-(exper@trans@n.epi.class*exper@trans@n.age.class)+1
                  state[first.index.tmp:last.index.tmp] <- exper@pop.rescale.each.timestep[s,t]*
                    (state[first.index.tmp:last.index.tmp]/sum(state[first.index.tmp:last.index.tmp]))
                }
              }

              #put in the correct birth rate for that time-step, if it varies
              for (s in 1:n.subpops) {
                if (ncol(exper@births.per.1000.each.timestep)>1) {
                  last.index.tmp <- (s*exper@trans@n.epi.class*exper@trans@n.age.class)
                  first.index.tmp <- last.index.tmp-(exper@trans@n.epi.class*exper@trans@n.age.class)+1
                  tmp.trans@birth.rate[s] = (exper@births.per.1000.each.timestep[s,t]*
                                               sum(state[first.index.tmp:last.index.tmp])/1000)
                } else {
                  tmp.trans@birth.rate[s] <- exper@trans@birth.rate[s]
                }
              }
              births.each.timestep[,t] <-  tmp.trans@birth.rate

              #put in time appropriate survival rate otherwise it uses original surv.matrix set up for trans object and keeps constant over time
              if (!is.na(exper@surv.each.timestep[1,1])) {
                for (s in 1:n.subpops) {
                  last.index.tmp1 <- (s*exper@trans@n.age.class)
                  first.index.tmp1 <- last.index.tmp1-(exper@trans@n.age.class)+1
                  tmp.age.surv.matrix <- ExtractAgeSpecificSurvivalMatrix(tmp.trans=tmp.trans, maternal.obj=exper@maternal.obj,
                                                                          surv.at.timestep.t=exper@surv.each.timestep[first.index.tmp1:last.index.tmp1,t])

                  last.index.tmp2 <- (s*exper@trans@n.epi.class*exper@trans@n.age.class)
                  first.index.tmp2 <- last.index.tmp2-(exper@trans@n.epi.class*exper@trans@n.age.class)+1
                  tmp.trans@age.surv.matrix[first.index.tmp2:last.index.tmp2,] <- tmp.age.surv.matrix
                }
              }

              ##put in correct vaccination coverage
              for (s in 1:n.subpops) {
                last.index.tmp <- (s*exper@trans@n.age.class)
                first.index.tmp <- last.index.tmp-(exper@trans@n.age.class)+1

                routine.vacc.prob <- routine$age.time.specific.routine[first.index.tmp:last.index.tmp,]
                sia.vacc.prob <- SIA$age.time.specific.SIA[first.index.tmp:last.index.tmp,]
                school.vacc.prob <- schoolvacc$age.time.specific.SIA[first.index.tmp:last.index.tmp,]

                if (!is.na(index.sia.vacc[t]) & is.na(index.school.vacc)){ #if SIA and no school vacc
                  tmp.trans@vac.per@pvacc.in.age.class[first.index.tmp:last.index.tmp] <-
                    routine.vacc.prob[,index.routine.vacc[t]] +
                    sia.vacc.prob[,index.sia.vacc[t]] -
                    (routine.vacc.prob[,index.routine.vacc[t]]*sia.vacc.prob[,index.sia.vacc[t]])
                  #stow output
                  MR1.fail.each.timestep[t] <- routine$prop.fail.MR1[s,index.routine.vacc[t]]
                  MR2.fail.each.timestep[t] <- routine$prop.fail.MR2[s,index.routine.vacc[t]]
                  SIA.fail.each.timestep[t] <- SIA$prop.fail.SIA[s,index.sia.vacc[t]]

                } else if ((is.na(index.sia.vacc[t]) & is.na(index.school.vacc))) { #if no SIA AND no school vacc
                  tmp.trans@vac.per@pvacc.in.age.class[first.index.tmp:last.index.tmp] <-
                    routine.vacc.prob[,index.routine.vacc[t]]
                  #stow output
                  MR1.fail.each.timestep[t] <- routine$prop.fail.MR1[s,index.routine.vacc[t]]
                  MR2.fail.each.timestep[t] <- routine$prop.fail.MR2[s,index.routine.vacc[t]]

                } else if (is.na(index.sia.vacc[t]) & !is.na(index.school.vacc)) {#if no SIA and school vacc
                  tmp.trans@vac.per@pvacc.in.age.class[first.index.tmp:last.index.tmp] <-
                    routine.vacc.prob[,index.routine.vacc[t]] +
                    school.vacc.prob[,index.school.vacc[t]] -
                    (routine.vacc.prob[,index.routine.vacc[t]]*school.vacc.prob[,index.school.vacc[t]])
                  #stow output
                  MR1.fail.each.timestep[t] <- routine$prop.fail.MR1[s,index.routine.vacc[t]]
                  MR2.fail.each.timestep[t] <- routine$prop.fail.MR2[s,index.routine.vacc[t]]
                  schoolvacc.fail.each.timestep[t] <- schoolvacc$prop.fail.schoolvacc[s,index.school.vacc[t]]

                } else if (!is.na(index.sia.vacc[t]) & !is.na(index.school.vacc)){ #if SIA and school vacc
                  routine.sia.tmp <-  routine.vacc.prob[,index.routine.vacc[t]] +
                    sia.vacc.prob[,index.sia.vacc[t]] - (routine.vacc.prob[,index.routine.vacc[t]]*sia.vacc.prob[,index.sia.vacc[t]])
                  tmp.trans@vac.per@pvacc.in.age.class[first.index.tmp:last.index.tmp] <-
                    routine.sia.tmp + school.vacc.prob[,index.school.vacc[t]] - (routine.sia.tmp*school.vacc.prob[,index.school.vacc[t]])
                  #stow output
                  MR1.fail.each.timestep[t] <- routine$prop.fail.MR1[s,index.routine.vacc[t]]
                  MR2.fail.each.timestep[t] <- routine$prop.fail.MR2[s,index.routine.vacc[t]]
                  SIA.fail.each.timestep[t] <- SIA$prop.fail.SIA[s,index.sia.vacc[t]]
                  schoolvacc.fail.each.timestep[t] <- schoolvacc$prop.fail.schoolvacc[s,index.school.vacc[t]]
                }
              }

              #stow the previous population size
              for (s in 1:n.subpops) {
                last.index.tmp <- (s*exper@trans@n.epi.class*exper@trans@n.age.class)
                first.index.tmp <- last.index.tmp-(exper@trans@n.epi.class*exper@trans@n.age.class)+1
                N0[s,t] <- sum(state[first.index.tmp:last.index.tmp])
              }

              #run experiment to get the next time step
              state <- next.ID.state(state, tran=tmp.trans)

              #stow the previous population size
              for (s in 1:n.subpops) {
                last.index.tmp <- (s*exper@trans@n.epi.class*exper@trans@n.age.class)
                first.index.tmp <- last.index.tmp-(exper@trans@n.epi.class*exper@trans@n.age.class)+1
                NT <- sum(state[first.index.tmp:last.index.tmp])
                growth.rate.each.timestep[s,t] <- log(NT/N0[s,t]) #instantaneous biweekly growth rate
              }

              #print(t)
              #print(dim(state))
              rc [,t] <- state
            }


            result <- new("sim.results.MSIRV.update.demog.vaccine.change.space.school",
                          data=rc,
                          m.inds = exper@trans@m.inds,
                          s.inds = exper@trans@s.inds,
                          i.inds = exper@trans@i.inds,
                          r.inds = exper@trans@r.inds,
                          v.inds = exper@trans@v.inds,
                          t = exper@t.min+(1:T-1)* exper@step.size,
                          age.class = exper@trans@age.class,
                          births.each.timestep = births.each.timestep,
                          growth.rate.each.timestep = growth.rate.each.timestep,
                          MR1.fail.each.timestep = MR1.fail.each.timestep,
                          MR2.fail.each.timestep = MR2.fail.each.timestep,
                          SIA.fail.each.timestep = SIA.fail.each.timestep,
                          schoolvacc.fail.each.timestep = schoolvacc.fail.each.timestep,
                          routine.intro = routine.intro,
                          sia.times = sia.times,
                          schoolvacc.times = schoolvacc.times,
                          n.subpops = exper@trans@n.subpops,
                          subpop.class.label = exper@trans@subpop.class.label)

            rc <- new("experiment.result",
                      experiment.def = exper,
                      result = result)


            return(rc)
          }
)









