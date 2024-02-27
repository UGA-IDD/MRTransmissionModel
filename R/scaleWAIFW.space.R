#' Function to scale the WAIFW matrix to a particular R0 given a particular subpopulation, but for a spatially structured population
#'
#' @param R0 the reproductive rate to scale to
#' @param state a population at disease free equilibrium
#' @param waifw the matrix to scale
#' @param frequency.dep xxx
#' @param suscept.state xxx
#'
#' @return a scaled versions of the age-specific waifw
#' @export
#'

scaleWAIFW.space <- function(R0, state, waifw, frequency.dep=F, suscept.state = 1) {

  #move everyone into susceptible category for DFE
  DFE.state <- matrix(0,length(state[,1]),1)
  for (k in  1:state@n.epi.class)
    DFE.state[state@epi.class==suscept.state,1] <-
      DFE.state[state@epi.class==suscept.state,1]+state[state@epi.class==k,1]

  #state.here <- rowSums(matrix(DFE.state[state@epi.class==suscept.state,1],
  #                     state@n.age.class,length(state[state@epi.class==suscept.state,1])/state@n.age.class))

  #just use first site (not whole pop)
  state.here <- DFE.state[state@epi.class==suscept.state,1][1:state@n.age.class]
  if (frequency.dep) denom <- sum(state.here) else denom <- 1
  #plot(state.here,type="l")

  #more correct
  next.gen <- state.here*(1-exp(-waifw/denom))

  #get the first eigen value
  cur.R0 <- Re(eigen(next.gen)$value[1])

  #More correct transform
  R.ratio <- R0/cur.R0; #print(R0); #print(cur.R0); #print(R.ratio)
  waifw <- -log(1-R.ratio*(1-exp(-waifw/denom)))*denom

  return(waifw)
}
