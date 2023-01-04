#' Last and final transition function to be called
#'
#' @param state xxx
#' @param tran xxx
#'
#' @importFrom stats rmultinom rpois
#'
#' @return xxx
#' @export
#'

SIRIDTran <- function (state, tran) {

  #define denominator of phi matrix
  if (tran@frequency.dep) {denom<-sum(state)} else {denom<-1}

  ## ## #define the phi matrix
  ## phi <- (tran@waifw%*%(state[tran@i.inds,]^tran@exponent))/denom

  ## phi <- 1 - exp(-phi)

  ## phi <- matrix(phi, nrow=state@n.age.class ,
  ##               ncol=state@n.age.class)

  ## phi <- t(phi)


  ## #get the transition matrix
  ## tran.matrix <- tran@age.surv.matrix

  ## #make susceptible part of matrix
  ## tran.matrix[tran@s.inds, tran@s.inds] <-
  ##     tran.matrix[tran@s.inds, tran@s.inds] * (1-phi)

  ## #make susceptible part of matrix
  ## tran.matrix[tran@i.inds, tran@s.inds] <-
  ##     tran.matrix[tran@i.inds, tran@s.inds] * (phi)

  ## #no one stays infecrted...might get more sophisticated later
  ## tran.matrix[tran@i.inds, tran@i.inds] <- 0


  tran.matrix <- .Call("calc_phi_and_update_tran",
                       tran@waifw,
                       state[,1],
                       tran@s.inds,
                       tran@i.inds,
                       tran@exponent,
                       denom,
                       tran@age.surv.matrix)


  if (!tran@is.stochastic) {

    birthst <- tran@get.births(state[,1],tran)
    #print(tran@get.births)
    #print(tran@birth.rate)
    #print(birthst[1])#xxjnow

    state[,] <- tran.matrix%*%state

    #add in the births to 1,1 for now, assuming that is correct.
    #might go back on this later.
    #state[1,1] <-  state[1,1] + tran@birth.rate #A BIT OF A HACK
    state[,1] <-  state[,1] + birthst #xxj - current births based on function
    #crude way of handling introductions
    state[tran@i.inds,] <- state[tran@i.inds,] +
      tran@introduction.rate

    #print(range(state[tran@i.inds,]))

  } else {
    birthst <- tran@get.births(state[,1],tran)
    #mortality probablity in each category of the n.age.class * no classes
    mort <- 1-rep(tran@survival.rate, each=state@n.epi.class)

    state[,1] <- .Call("do_ID_transition_SIR_stochastic_moves_cl",
                       as.integer(state[,1]),
                       tran.matrix)

    # # do you always want to use this function?
    # # remove the if statement
    # if (is.loaded("do_ID_transition_SIR_stochastic_moves")) {
    #   #C implementation of this.
    #   #state[,1] <- .C("do_ID_transition_SIR_stochastic_moves",
    #   #                as.integer(state[,1]),
    #   #                as.integer(length(state[,1])),
    #   #                as.double(tran.matrix),
    #   #                newstate=as.integer(state[,1]))$newstate
    #
    #
    #   state[,1] <- .Call("do_ID_transition_SIR_stochastic_moves_cl",
    #                    as.integer(state[,1]),
    #                    tran.matrix)
    #
    #   # ALWAYS use the C code
    #
    # } else {
    #   #loop over and distribute via a multinomial
    #   newstate <- rep(0,1+length(tran.matrix[,1]))
    #   for (k in 1:length(tran.matrix[1,])) {
    #     newstate <- newstate+
    #       rmultinom(1,state[k,1],c(tran.matrix[,k],mort[k]))
    #   }
    #
    #   state[,1] <- newstate[1:length(tran.matrix[,1])]
    #
    # }

    #add in the births to 1,1 for now, assuming that is correct.
    #might go back on this later.
    #state[1,1] <-  state[1,1] + rpois(1,tran@birth.rate) #A BIT OF A HACK
    state[,1] <- state[,1] +rpois(length(birthst),birthst)#xxj - current births based on function
    #crude way of handling introductions
    state[tran@i.inds,1] <- state[tran@i.inds,1] +
      rpois(length(tran@i.inds),tran@introduction.rate)
  }
  return(state)
}
