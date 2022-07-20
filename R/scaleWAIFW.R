#' Function to scale the WAIFW matrix to a particular R0 given a particular population
#'
#' @param R0 the reproductive rate to scale to
#' @param state a population at disease free equilibrium
#' @param waifw the matrix to scale
#' @param frequency.dep xxx
#' @param suscept.state xxx
#'
#' @return a scaled vertions of waifw
#' @export
#'

scaleWAIFW <- function(R0, state, waifw, frequency.dep=F,suscept.state = 1) {

  if (frequency.dep) denom <- sum(state[,1]) else denom <- 1

  #print(length(state[,1]))

  #move everyone into susceptible category for DFE
  if (length(state@epi.class)==length(state[,1])) {
    DFE.state <- matrix(0,length(state[,1]),1)
    for (k in  1:state@n.epi.class)
      DFE.state[state@epi.class==suscept.state,1] <-
        DFE.state[state@epi.class==suscept.state,1]+state[state@epi.class==k,1]
  } else {
    epi.class.here <- rep(state@epi.class,2)
    DFE.state <- matrix(0,length(epi.class.here),1)
    for (k in  1:state@n.epi.class) {
      DFE.state[epi.class.here == suscept.state,] <- DFE.state[epi.class.here == suscept.state,]+
        state[epi.class.here == k,]
      # print(k)
    }
    #add the genders together
    DFE.state <- DFE.state[1:length(state@epi.class),]+DFE.state[(length(state@epi.class)+1):(2*length(state@epi.class)),]
    DFE.state <- matrix(DFE.state,length(state@epi.class),1)
  }


  #worked kinda
  #next.gen <- (state[state@epi.class==suscept.state,1]*waifw)/denom

  #plot(DFE.state[state@epi.class==suscept.state,1],type="l")


  #more correct
  next.gen <- DFE.state[state@epi.class==suscept.state,1]*(1-exp(-waifw/denom))

  #get the first eigen value
  cur.R0 <- Re(eigen(next.gen)$value[1])

  #More correct transform
  R.ratio <- R0/cur.R0; #print(R0); #print(cur.R0); #print(R.ratio)
  waifw <- -log(1-R.ratio*(1-exp(-waifw/denom)))*denom

  return(waifw)
}
