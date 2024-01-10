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

