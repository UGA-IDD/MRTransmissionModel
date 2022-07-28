#### Getting State at Next Time Step (next.ID.state methods) ####

#Method that takes an ID.state.matrix and an ID.transition object
#and create the next state of the ID.state.matrix
#
#Parameters -
#   state - current system state (ID.state.matrix)
#   tran - ID transition object that has all of the infor supervising the transition
#
#Returns -
#  an updated ID.state.matrix
setGeneric("next.ID.state",
           function(state, tran, ...) standardGeneric("next.ID.state"))

# default version of next.ID.state
setMethod("next.ID.state",
          c("ID.state.matrix", "ID.transition.SIR"),
          SIRIDTran)

# next ID state for ID.transition.SIR.vac class
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

#next ID state for MSIRV class
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

            rc <- callNextMethod(state, tran)
            #print(rc[1:5,1])
            return(rc)
          })
