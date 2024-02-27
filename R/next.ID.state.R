#' Methods for Getting State at Next Time Step (next.ID.state methods)
#'
#'
#'
#' @param state current state system (ID.state.matric)
#' @param tran ID transition object that has all info for supervising the transition
#' @param ... ignored
#'
#' @importFrom methods is
#'
#' @include setClasses.R
#' @return an updated ID.state.matrix
#'
#' @export
#'
#' @docType methods
#' @rdname next.ID.state-methods
#'

setGeneric("next.ID.state",
           function(state, tran, ...) standardGeneric("next.ID.state"))

#' @rdname next.ID.state-methods
#' @aliases next.ID.state,ID.state.matrix,ID.transition.SIR-method
setMethod("next.ID.state",
          c("ID.state.matrix", "ID.transition.SIR"),
          SIRIDTran)

#' @rdname next.ID.state-methods
#' @aliases next.ID.state,ID.state.matrix,ID.transition.SIR.vac-method
setMethod("next.ID.state",
          c("ID.state.matrix", "ID.transition.SIR.vac"),
          function(state, tran) {

            #A bit of a hack, so we do not double deal with vaccination
            if (!is(tran, "ID.transition.MSIRV")) {
              delta.state <-  state[tran@s.inds]*
                (tran@vac.per@pvacc.in.age.class)

              state[tran@s.inds] <- state[tran@s.inds] - delta.state
              state[tran@r.inds] <- state[tran@r.inds] + delta.state

              mid.age.class <- (tran@age.class + c(0, tran@age.class[2:tran@n.age.class-1]))/2
              sia.sv.prob <-tran@sia.vac * pvacsuccess(mid.age.class, tran@sia.vsucc)

              delta.sia <- state[tran@s.inds]*sia.sv.prob
              state[tran@s.inds] <- state[tran@s.inds] - delta.sia
              state[tran@r.inds] <- state[tran@r.inds] + delta.sia
            }

            #rc <- callNextMethod(state, tran)
            rc <- SIRIDTran(state,tran) #avoid callNextMethod overhead
            return(rc)
          })

#' @rdname next.ID.state-methods
#' @aliases next.ID.state,ID.state.matrix,ID.transition.MSIRV-method
setMethod("next.ID.state",
          c("ID.state.matrix", "ID.transition.MSIRV"),
          function (state, tran) {

            tran@age.surv.matrix <- .Call("update_age_surv_MSIRV",
                                          tran@age.surv.matrix,
                                          sz = nrow(tran@age.surv.matrix),
                                          tran@vac.per@pvacc.in.age.class,
                                          tran@v.inds,
                                          tran@m.inds,
                                          tran@s.inds,
                                          tran@i.inds)

            rc <- SIRIDTran(state,tran) #avoid callNextMethod overhead
            # can we change this too?
            #rc <- callNextMethod(state, tran)
            #print(rc[1:5,1])
            return(rc)
          })


# next ID state for MSIRV.space class
#' @rdname next.ID.state-methods
#'
#'
#' @aliases next.ID.state,ID.state.matrix,ID.transition.MSIRV.space-method
setMethod("next.ID.state",
          c("ID.state.matrix", "ID.transition.MSIRV.space"),
          function (state, tran) {

            #get the transition matrix
            tran.matrix <- tran@age.surv.matrix

            #get a spatial index on same scale as epi indexes
            index.loc <-  tran@subpop.class.label[tran@i.inds]
            one.loc <- index.loc==1

            for (n in 1:tran@n.subpops) {
              this.loc <- index.loc==n

              #The vaccination logic
              age.spec.vacc.prob <- tran@vac.per@pvacc.in.age.class[this.loc]

              ##Updating transition matrix to include vaccination
              #M->V transition
              tran.matrix[tran@v.inds[this.loc],tran@m.inds[one.loc]] <-  tran@age.surv.matrix[tran@v.inds[this.loc],tran@m.inds[one.loc]] *
                age.spec.vacc.prob
              #M->M transition
              tran.matrix[tran@m.inds[this.loc],tran@m.inds[one.loc]] <-  tran@age.surv.matrix[tran@m.inds[this.loc],tran@m.inds[one.loc]] *
                (1-age.spec.vacc.prob)
              #M->S transition
              tran.matrix[tran@s.inds[this.loc],tran@m.inds[one.loc]] <-  tran@age.surv.matrix[tran@s.inds[this.loc],tran@m.inds[one.loc]] *
                (1-age.spec.vacc.prob)
              #S->V transition
              tran.matrix[tran@v.inds[this.loc],tran@s.inds[one.loc]] <-  tran@age.surv.matrix[tran@v.inds[this.loc],tran@s.inds[one.loc]] *
                age.spec.vacc.prob
              #S->S transition
              tran.matrix[tran@s.inds[this.loc],tran@s.inds[one.loc]] <-  tran@age.surv.matrix[tran@s.inds[this.loc],tran@s.inds[one.loc]] *
                (1-age.spec.vacc.prob)
              #S->I transition
              tran.matrix[tran@i.inds[this.loc],tran@s.inds[one.loc]] <-  tran@age.surv.matrix[tran@i.inds[this.loc],tran@s.inds[one.loc]] *
                (1-age.spec.vacc.prob)


              ## Calculate the phi matrix
              #define denominator of phi matrix - preventing NAs by setting min to 1
              if (tran@frequency.dep) {denom<-max(sum(state[tran@subpop.class.label==n]),1)} else {denom<-1}

              phi <- tran@waifw%*%(state[tran@i.inds[this.loc],]^tran@exponent)/denom
              phi <- 1 - exp(-phi)
              #hist(phi, breaks=20, xlim=c(0,1), main=n)
              phi <- matrix(phi, nrow=state@n.age.class, ncol=state@n.age.class)
              phi <- t(phi)

              #print(c("phi",range(phi)))

              ##Now that phi is calculated, update the tran matrix
              #people who stay susceptible
              tran.matrix[tran@s.inds[this.loc], tran@s.inds[one.loc]] <-
                tran.matrix[tran@s.inds[this.loc], tran@s.inds[one.loc]] * (1-phi)

              #people who become infected
              tran.matrix[tran@i.inds[this.loc], tran@s.inds[one.loc]] <-
                tran.matrix[tran@i.inds[this.loc], tran@s.inds[one.loc]] * (phi)

              #no one stays infected
              tran.matrix[tran@i.inds, tran@i.inds[one.loc]] <- 0

              #parametric fit
              #prop local susceptibility (age specific)
              age.struct.here <- (state[tran@m.inds[this.loc]]+state[tran@s.inds[this.loc]]+
                                    state[tran@i.inds[this.loc]]+state[tran@r.inds[this.loc]]+
                                    state[tran@v.inds[this.loc]])
              xtj <- (state[tran@s.inds[this.loc]]/pmax(age.struct.here,1))
              xtj <- xtj/length(tran@introduction.rate[this.loc]) #rescale since have multiple age classes
              #print(c("xtj",range(xtj)))
              #print(range(age.struct.here))

              #prop non-local infectious - sum over age classes - broken down by location
              #since need to multiply by each coupling par
              yt <- sapply(split(state[tran@i.inds[!this.loc]],index.loc[!this.loc]),sum)/
                pmax(sapply(split(state[tran@m.inds[!this.loc]]+state[tran@s.inds[!this.loc]]+
                                    state[tran@i.inds[!this.loc]]+state[tran@r.inds[!this.loc]]+
                                    state[tran@v.inds[!this.loc]],index.loc[!this.loc]),sum),1)

              #ADD NEW GUYS TO THE INTRODUCTION.RATE FOR THIS SITE - at the moment distributed evenly over age
              if (!tran@is.stochastic) {
                tran@introduction.rate[this.loc] <-  tran@introduction.rate[this.loc] +
                  (1-exp(-sum(tran@coupling[n,-n]*yt)*xtj))
              } else {
                immigs <- rbinom(length(xtj),1,pmax(1-exp(-sum(tran@coupling[n,-n]*yt)*xtj),0))
                #print(yt); print(xtj); print(tran@coupling[n,-n])
                #print(immigs)
                if (sum(immigs)>0) {
                  tran@introduction.rate[this.loc] <-
                    rbinom(length(tran@i.inds[this.loc]),1,tran@introduction.rate[this.loc])+immigs

                } else {
                  tran@introduction.rate[this.loc] <-
                    rbinom(length(tran@i.inds[this.loc]),1,tran@introduction.rate[this.loc])
                }
              }
            }



            #print(c(unique(tran.matrix)))

            if (!tran@is.stochastic) {

              #loop over sites and multiply matrices
              for (n in 1:tran@n.subpops) {
                here <- which(tran@subpop.class.label==n,arr.ind=TRUE)
                state[here,] <- tran.matrix[here,]%*%state[here,]
              }

              #add in the births to 1,1 for now, assuming that is correct.
              #might go back on this later.
              state[!duplicated(tran@subpop.class.label),1] <-  state[!duplicated(tran@subpop.class.label),1] +
                tran@birth.rate #A BIT OF A HACK

              #crude way of handling introductions and emigrations
              state[tran@i.inds,1] <- state[tran@i.inds,1] + tran@introduction.rate

            } else {

              #mortality probability in each category of the n.age.class * no classes
              mort <- 1-rep(tran@survival.rate, each=state@n.epi.class)


              #loop over and distribute via a multinomial
              newstate <- rep(0,1+length(tran.matrix[,1]))

              for (n in 1:tran@n.subpops) {

                here <- which(tran@subpop.class.label==n,arr.ind=TRUE)
                here.next <- c(here,length(newstate))


                #  print(here)

                for (k in 1:length(tran.matrix[1,])) {

                  #print("living, then dying")
                  #print(sum(tran.matrix[here,k]))
                  #print(mort[((n-1)*tran@n.age.class)+k])

                  newstate[here.next] <- newstate[here.next] +
                    rmultinom(1,state[here[k],1],c(tran.matrix[here,k],mort[((n-1)*tran@n.age.class)+k]))
                }
              }

              #print("prop dead")
              #print(newstate[length(newstate)]/sum(newstate))

              state[,1] <- newstate[1:length(tran.matrix[,1])]

              #print(unique(state[,1]))

              #print("before births")
              #print(range(c(state@.Data)))

              #add in the births to 1,1 for now, assuming that is correct.
              state[!duplicated(tran@subpop.class.label),1] <-  state[!duplicated(tran@subpop.class.label),1] +
                rpois(tran@n.subpops,tran@birth.rate)


              #print("after births")
              #print(range(c(state@.Data)))

              #print(unique(tran@introduction.rate))


              #the stoch dist introduced above for efficiency
              state[tran@i.inds,1] <- state[tran@i.inds,1] + tran@introduction.rate



            }


            return(state)
          })




