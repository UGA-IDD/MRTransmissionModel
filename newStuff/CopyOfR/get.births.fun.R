#' This function is an argument of the experiment builder to ensure fertility from susceptible mothers is susceptible
#'
#' @param state state vector at this stage in simulation, length n.epi.class*n.age.class
#' @param tran transition matrix at this stage in simulation
#' @param fert.curve age specific fertility rate used to establish fraction susceptible mothers
#' @param lower.age.boundary corresponding lower age boundaries for the fert.curve
#'
#' @return vector of births of length state
#' @export
#'

get.births.fun <- function(state,tran,
                           fert.curve=c(0.0,36.3,210.6,164.8,68.3,27.5,10.1,2.8,0.0),
                           lower.age.boundary = c(0,15,20,25,30,35,40,45,50,101)){ #xxamy - age.class change from lower.age.boundary = c(0,15,20,25,30,35,40,45,50,100)){

  # Births based on transition object which is changing at each time point
  births <- tran@birth.rate
  #print(births)

  # Establish fertility curve
  n.age.cats <- length(tran@age.class)
  ages.to.use <- (tran@age.class + c(0, tran@age.class[2:n.age.cats-1]))/2
  fertility.category <- findInterval(ages.to.use/12,lower.age.boundary) #tells which fert.cat to pull from for each age cat for all 280 age cats.  so 1 through 180 (12 months * 15 years) is 1 or fertilty rate of 0

  # Get the proportion of susceptible mothers
  #the age specific fertility rate is helping with is the determining the proportion of births that are susceptible
  #otherwise the crude birth rate is used to determine number of births
  totpop <- state[tran@m.inds]+state[tran@s.inds]+state[tran@i.inds]+state[tran@r.inds]+state[tran@v.inds] #total pop per age category (280 age cats)
  prop.susc.mothers <- sum(state[tran@s.inds]*fert.curve[fertility.category])/sum(totpop*fert.curve[fertility.category]) #proportion from 0 to 1
  #print(prop.susc.mothers)

  # Make the new baby vector - rescaling birth rate for change in births over time
  new.babies <- rep(0,length(state))
  new.babies[1] <- (1-prop.susc.mothers)*births #just dividing the births into S class and M class
  new.babies[2] <- prop.susc.mothers*births

  #print(new.babies[1:10])
  return(new.babies)

}
