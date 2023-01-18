## Amy Winter
## 26 June 2022
## Tasks still to do
#1. below are DFE - my DFE function isn't real as it is just DF not DFE.  Consider revising
#2. In setMethod("run", "experiment.updatedemog.vaccinationchange.vaccinationlimitations".. something is wrong b/c prop.inacc option without inefficiency option has lower prop vaccinated by age than with inefficiency option
#3. The correlation code isn't correct. basically in the c code prop susceptible getting vaccinated is the same as the prop population getting vaccinated.  This isn't true if correlation b/w doses.
## the code runs, but it isn't right. I would need new c code to allow for this correctly ...
# 4. spatial code: get.births.here() is likely an issue that needs t be fixed once try to estimate CRS
# 5. good bit of code in next.ID.state() that *works* but I need to review better



#Creates a state matrix at the disease free equilibrium
#from a transition object and a population size
#
#Parameters -
#     pop - the population size
#     tran - the transition object
#     ... - other arguments to matrix creation
#
#Return -
#     an ID.state.matrix
create.DFE.ID.state.matrix <- function(pop, tran, ...) {
  I <- diag(tran@n.age.class)
  
  #birth vector
  B <- c(1,rep(0,tran@n.age.class-1))
  #B <- c(tran@birth.rate,rep(0,tran@n.age.class-1)) #debug
  
  #make the population distribution matrix
  dist <- solve(I-tran@age.surv.matrix[tran@s.inds,tran@s.inds])%*%B
  dist <- dist/sum(dist)
  
  
  
  rc <- create.ID.state.matrix(tran@n.age.class,
                               tran@n.epi.class, ...)
  rc[tran@s.inds,1] <- dist*pop
  
  return(rc)
  
}

#Creates a state matrix at the disease free equilibrium
#from a transition object and a population size
#for a growing population by iterating forwards to stability
#
#Parameters -
#     pop - the population size
#     tran - the transition object
#     ... - other arguments to matrix creation
#
#Return -
#     an ID.state.matrix
create.DFE.ID.state.matrix.growing <- function(pop, tran, ...) {
  I <- diag(tran@n.age.class)
  
  #birth vector
  B <- c(tran@birth.rate,rep(0,tran@n.age.class-1))
  #print(B[1])
  
  #make the population distribution matrix for stable to kickoff
  dist <- solve(I-tran@age.surv.matrix[tran@v.inds,tran@v.inds])%*%B
  dist <- dist/sum(dist)
  
  #create the state matrix
  rc <- create.ID.state.matrix(tran@n.age.class,
                               tran@n.epi.class, ...)
  #iterate to get stability
  popstart <- dist*pop
  for (t in 1:30000) {
    popstart1 <- tran@age.surv.matrix[tran@v.inds,tran@v.inds]%*%popstart
    popstart <- popstart1 + (B*sum(popstart)/pop)
  }
  
  #introduce
  dist <- popstart/sum(popstart)
  rc[tran@s.inds,1] <- dist*pop
  
  return(rc)
  
}


#Create disease equilibrium starting state from a
#DFE state
#
#Parameters -
#   state - the DFE state
#   tran - the scaled transmission matrix
#
#Return -
#   a population with the number of susceptibles at approximate
#   disease equilibrium (ignoring seasonality and all that jazz)
#
create.DE.ID.state.matrix <- function (state, tran) {
  
  #for a first pass, just worrk about the age to age beta
  foi <- diag(tran@waifw)*state[tran@s.inds,1]
  cum.foi <- cumsum(foi)
  
  psuscept <- exp(-cumsum(cum.foi))
  
  state[tran@r.inds,1] <- state[tran@s.inds,1]*(1-psuscept)
  state[tran@s.inds,1] <- state[tran@s.inds,1]*(psuscept)
  
  return(state)
}