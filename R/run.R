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
                  tmp.trans@vac.per@pvacc.in.age.class <- 1-(exper@prop.inacc[t]+
                                                               (1-exper@prop.inacc[t])*(1-(routine.vacc.prob[index.routine.vacc[t],] +
                                                                                             sia.vacc.prob[index.sia.vacc[t],] +
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










