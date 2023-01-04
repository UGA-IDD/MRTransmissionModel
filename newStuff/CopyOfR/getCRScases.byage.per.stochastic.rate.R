#' Function to get the total number of CRS pregnancies by age
#'
#' @param sim.res a sim.result.MSIRV object
#' @param crs_rate_gestational_age_nbiweeks vector; gestational ages (in biweeks) that coincides with crs_rate
#' @param crs_rate vector; probability crs case given rubella infection in pregnancy and live birth that coincides with crs_rate_gestational_age_nbiweeks
#' @param trans xxx
#' @param no.gens.in.year xxx
#' @param births vector; number of births per time period taken from the simulation
#' @param asfr.over.time matrix; age-specific fertility rate in every age class and year (nrow=7, ncol=length of sim in years)
#' @param repro.age.sex.dist.over.time xxx
#'
#' @return an object with age classes set by fertility age classes
#' @export
#'

getCRScases.byage.per.stochastic.rate <- function(sim.res, crs_rate_gestational_age_nbiweeks=c(6,8),
                                                  crs_rate=c(0.4, 0.15),
                                                  trans, no.gens.in.year, births,
                                                  asfr.over.time, repro.age.sex.dist.over.time) {

  crs.age.time <- matrix(0,nrow=35,ncol=ncol(sim.res@.Data))
  for (a in 1:length(crs_rate)){

    if (a==1){
      # Get a vector providing risk for every age class of susceptibles
      CRSrisk <- getCumRiskTrimester(sim.res=sim.res, nbiweeks=crs_rate_gestational_age_nbiweeks[a])
      #plot(colSums(CRSrisk),type="l")

      # Get the number of individuals who become sick in nbiweeks
      SickInFollowingTrim <- sim.res[sim.res@s.inds,]*CRSrisk
      #plot(colSums(SickInFollowingTrim), type="l")

    } else { #if a>1
      # Get a vector providing risk for every age class of susceptibles
      CRSrisk.y <- getCumRiskTrimester(sim.res=sim.res, nbiweeks=crs_rate_gestational_age_nbiweeks[(a-1)])
      CRSrisk.o <- getCumRiskTrimester(sim.res=sim.res, nbiweeks=crs_rate_gestational_age_nbiweeks[a])
      CRSrisk <- CRSrisk.o-CRSrisk.y

      # Get the number of individuals who become sick in nbiweeks
      SickInFollowingTrim <- sim.res[sim.res@s.inds,]*CRSrisk
    }

    # Sum the 15-19 year age groups (in months) into one year age groups
    index.15.to.20 <- findInterval((sim.res@age.class-1)/12,c(15,16,17,18,19,20,101))+1
    tmp.15.to.20 <- matrix(NA, nrow=length(unique(index.15.to.20)), ncol=ncol(SickInFollowingTrim))
    for (i in 1:length(unique(index.15.to.20))){
      tmp.15.to.20[i,] <- apply(SickInFollowingTrim, 2, function(x) sum(x[index.15.to.20==i]))
    }
    SickInFollowingTrim.repro.ages <- rbind(tmp.15.to.20[-c(1,7),], SickInFollowingTrim[241:270,])

    #getting fertility indexing
    lower.age.boundary <- c(15,20,25,30,35,40,45,50) # xxamy - age.class change from lower.age.boundary <- c(0,15,20,25,30,35,40,45,49,100) # for asfr
    repro.ages <- 15:49
    fertility.category <- findInterval(repro.ages,lower.age.boundary)

    # Multiply to get number of pregnancies with CRS
    crs.age.time.tmp <- crs_rate[a]*(1/no.gens.in.year)*asfr.over.time[fertility.category,]*
      SickInFollowingTrim.repro.ages*repro.age.sex.dist.over.time[fertility.category,]

    #Adding to the number of CRS cases based on each gestational age group associated with a rate
    crs.age.time <- crs.age.time+crs.age.time.tmp

  }

  # Adjust the CRS totals given that births based on CBR not ASFR in the simulation
  #getting population size of women of childbearing age
  pop.age <- matrix(NA, (max(trans@age.class))/12, ncol(sim.res)) #xxamy - age.class change from pop.age <- matrix(NA, (max(trans@age.class)-1)/12, ncol(sim.res))
  for (t in 1:ncol(sim.res)){
    state <- sim.res[,t]
    vec0 <- GetNumber.per.AgeGroup(trans, state, tval=1)
    pop.age[,t] <- GetNumber.per.AgeYear(vec=vec0, age.classes=trans@age.class)
  }
  row.names(pop.age) <- 0:100
  repro.ages.index <- findInterval(0:100, seq(15,50,5))+1
  repro.age.pop <- matrix(NA, nrow=length(unique(repro.ages.index)), ncol=ncol(pop.age))
  for (i in 1:length(unique(repro.ages.index))){
    repro.age.pop[i,] <- apply(pop.age, 2, function(x) sum(x[repro.ages.index==i]))
  }
  repro.age.pop <- repro.age.pop[-c(1,9),] #narrowing down to only ages 15-49
  #getting births based on age-specific fertility rates
  births.asfr <- colSums(repro.age.pop*asfr.over.time/no.gens.in.year*repro.age.sex.dist.over.time)
  #getting scalar (births based on sim divided by births based on asfr)
  scale <- births/births.asfr
  #adjusting crs given scalar
  crs.adj.age.time <- crs.age.time*scale

  return(crs.adj.age.time)
}
