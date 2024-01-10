#' Function to get the VIMC output needed from each result object
#'
#' @param sim.res a sim.result.MSIRV object
#' @param t.max numeric; the t.max used to create sim.res
#' @param setup country specific object used to create sim.res
#' @param year numeric; used to create sim.res
#' @param crs_rate_gestational_age_nbiweeks vector; gestational ages (in biweeks) that coincides with crs_rate_agespecific
#' @param crs_rate_agespecific vector; probability crs case given rubella infection in pregnancy and live birth that coincides with crs_rate_gestational_age_nbiweeks
#' @param crs_rate_overall numeric; average CRS rate over age groups per crs_rate_gestational_age_nbiweeks
#' @param fetal_death_rate numeric; prob fetal death given rubella infection in pregnancy
#' @param infant_death_rate numeric; prob infant death given rubella infection in pregnancy
#' @param scenario xxx
#' @param run.number xxx
#' @param R0 xxx
#'
#' @importFrom stats smooth.spline
#'
#' @return a large list of matrices by age and time in biweeks (population, rubella cases, rubella seroprevalence, CRS cases, CRS deaths, CRS dalys, CRS rate)
#' @export
#'

getOutput_201910gavi_v4 <- function(sim.res, t.max, setup, year,
                                    crs_rate_gestational_age_nbiweeks,
                                    crs_rate_agespecific,
                                    crs_rate_overall,
                                    fetal_death_rate,
                                    infant_death_rate,
                                    scenario, run.number, R0){

  tmpr <- sim.res
  epi.state <- tmpr@experiment.def@trans@epi.class
  no.gens.in.year <- 1/tmpr@experiment.def@step.size

  #Population size by age cohort at mid-year (week 13 in a year)
  pop.1yage.1ytime <- GetPop.per.Year(res=tmpr@result@.Data, trans=tmpr@experiment.def@trans, no.gens.in.year)

  #Rubella cases by age cohort and year (sum over each)
  rubella.cases.1yage.1ytime <- GetRubellaCases.per.Year(res=tmpr@result@.Data, trans=tmpr@experiment.def@trans,
                                                         epi.state=epi.state, no.gens.in.year)

  #Seroprevalence by age cohort and year (mid-year seroprev)
  rubella.seropos.1yage.1ytime <- GetRubellaSeroprevalence.per.Year(res=tmpr@result@.Data, trans=tmpr@experiment.def@trans,
                                                                    epi.state=epi.state, no.gens.in.year)

  #CRS cases by womens repro age (15-49) and year
  ##setting up inputs to estimates CRS
  asfr.index <- which(colnames(setup$asfr.1950.2100)==as.character(year)):which(colnames(setup$asfr.1950.2100)==as.character(year+t.max-1))
  asfr.over.time <- matrix(apply(setup$asfr.1950.2100[,asfr.index], 1, function(x) rep(x, each=no.gens.in.year)), ncol=(no.gens.in.year*length(asfr.index)), byrow=T)
  asfr.over.time <- cbind(asfr.over.time[,1], asfr.over.time) #repeat first time point
  repro.age.sex.dist.over.time <- matrix(apply(setup$repro.age.sex.dist.1950.2100[,asfr.index], 1, function(x) rep(x, each=no.gens.in.year)), ncol=(no.gens.in.year*length(asfr.index)), byrow=T)
  repro.age.sex.dist.over.time <- cbind(repro.age.sex.dist.over.time[,1], repro.age.sex.dist.over.time) #repeat first row
  ##getting crs cases by biweek
  crs.1yage.biweek <- as.matrix(getCRScases.byage.per.stochastic.rate(sim.res=tmpr@result,
                                                                      crs_rate_gestational_age_nbiweeks=crs_rate_gestational_age_nbiweeks,
                                                                      crs_rate=crs_rate_agespecific,
                                                                      trans=tmpr@experiment.def@trans, no.gens.in.year,
                                                                      births=tmpr@result@births.each.timestep, asfr.over.time=asfr.over.time,
                                                                      repro.age.sex.dist.over.time=repro.age.sex.dist.over.time))

  ##getting crs cases by year
  crs.1yage.1ytime <- GetSumInYear(mat=crs.1yage.biweek, no.gens.in.year=no.gens.in.year)

  #checking for NAs because it means something went wrong - very important because I change all NAs to zero below
  if(sum(is.na(rubella.cases.1yage.1ytime)) > 0 | sum(is.na(rubella.seropos.1yage.1ytime)) > 0 |
     sum(is.na(crs.1yage.1ytime)) > 0) {
    stop(paste("NAs produced by the model \n",
               setup$iso3code, "no. rubella cases NAs:", sum(is.na(rubella.cases.1yage.1ytime)),
               "\n", setup$iso3code, "no. rubella seroprev NAs", sum(is.na(rubella.seropos.1yage.1ytime)),
               "\n",setup$iso3code, "no. CRS cases NAs", sum(is.na(crs.1yage.1ytime)),
               "\n", scenario, ", run number:", run.number, ", R0:", R0))
  }

  #Averaging over two year time points because biannual dynamics
  #Rubella Cases
  last.col <- ncol(rubella.cases.1yage.1ytime)
  rubella.cases.1yage.1ytime.orig <- rubella.cases.1yage.1ytime
  rubella.year.totals <- colSums(rubella.cases.1yage.1ytime)
  rubella.smooth.year.totals <- sapply(seq(2,(length(rubella.year.totals))), function(x) mean(rubella.year.totals[x:(x-1)]))
  rubella.cases.1yage.1ytime <- t(apply(rubella.cases.1yage.1ytime[,-1], 1, function(x) x/(rubella.year.totals[-1])*rubella.smooth.year.totals))
  rubella.cases.1yage.1ytime[is.na(rubella.cases.1yage.1ytime)] <- 0
  rubella.cases.1yage.1ytime <- cbind(rubella.cases.1yage.1ytime.orig[,1], rubella.cases.1yage.1ytime) #add back col 1
  #CRS Cases
  crs.1yage.1ytime.orig <- crs.1yage.1ytime
  crs.year.totals <- colSums(crs.1yage.1ytime)
  crs.smooth.year.totals <- sapply(seq(2,(length(crs.year.totals))), function(x) mean(crs.year.totals[x:(x-1)]))
  crs.1yage.1ytime <- t(apply(crs.1yage.1ytime[,-1], 1, function(x) x/(crs.year.totals[-1])*crs.smooth.year.totals))
  crs.1yage.1ytime[is.na(crs.1yage.1ytime)] <- 0
  crs.1yage.1ytime <- cbind(crs.1yage.1ytime.orig[,1], crs.1yage.1ytime) #add back col 1
  #Rubella Seroprevalence
  rubella.seropos.1yage.1ytime.orig <- rubella.seropos.1yage.1ytime
  rubella.year.totals <- colSums(rubella.seropos.1yage.1ytime)
  rubella.smooth.year.totals <- sapply(seq(2,(length(rubella.year.totals))), function(x) mean(rubella.year.totals[x:(x-1)]))
  rubella.seropos.1yage.1ytime <- t(apply(rubella.seropos.1yage.1ytime[,-1], 1, function(x) x/(rubella.year.totals[-1])*rubella.smooth.year.totals))
  rubella.seropos.1yage.1ytime[is.na(rubella.seropos.1yage.1ytime)] <- 0
  rubella.seropos.1yage.1ytime <- cbind(rubella.seropos.1yage.1ytime.orig[,1], rubella.seropos.1yage.1ytime) #add back col 1
  #IF vaccine introduced - only average pre-vaccine, keep the post-vaccine raw output
  if(length(which(colSums(tmpr@result@.Data[tmpr@result@v.inds,])!=0))!=0) {
    # don't smooth after vaccine introduction
    index.vaccine.intro <- min(which(colSums(tmpr@result@.Data[tmpr@result@v.inds,])!=0))
    year.vaccine.intro <- ceiling(index.vaccine.intro/no.gens.in.year)
    rubella.cases.1yage.1ytime <- cbind(rubella.cases.1yage.1ytime[,1:(year.vaccine.intro-1)], rubella.cases.1yage.1ytime.orig[,year.vaccine.intro:last.col])
    crs.1yage.1ytime <- cbind(crs.1yage.1ytime[,1:(year.vaccine.intro-1)], crs.1yage.1ytime.orig[,year.vaccine.intro:last.col])
    rubella.seropos.1yage.1ytime <- cbind(rubella.seropos.1yage.1ytime[,1:(year.vaccine.intro-1)], rubella.seropos.1yage.1ytime.orig[,year.vaccine.intro:last.col])
  }

  #crs incidence rates
  births <- tmpr@result@births.each.timestep
  sidx = seq.int(from=1, to=length(births), by=no.gens.in.year)
  eidx = c((sidx-1)[2:length(sidx)], length(births))
  births.year = sapply(1:length(sidx), function(i) sum(births[sidx[i]:eidx[i]]))[1:(t.max)]
  crs.rate <- 100000*colSums(crs.1yage.1ytime)/births.year

  #Deaths from CRS
  fetal.deaths <- crs.1yage.1ytime/crs_rate_overall*fetal_death_rate
  child.deaths <- crs.1yage.1ytime*infant_death_rate
  deaths.1yage.1ytime <- fetal.deaths+child.deaths

  #CRS DALYs
  le.data <- c(62, 66, 74, 79)
  gbd.2010 <- c(29.2, 27.8, 22.9, 19.2)
  le.1yrtime <- setup$e0.1950.2100[which(1950:2100==year):length(setup$e0.1950.2100)]
  daly.percase.2010 <- predict(smooth.spline(le.data, gbd.2010), le.1yrtime)
  daly.percase.2010$y[daly.percase.2010$x<le.data[1]] <- gbd.2010[1]
  daly.percase.2010$y[daly.percase.2010$x>le.data[4]] <- gbd.2010[4]
  dalys.1yage.1ytime <- t(apply(crs.1yage.1ytime,1,function(x) x*daly.percase.2010$y[1:(t.max)]))

  #Reducing to only data between 2000 and 2100
  pop.1yage.1ytime <- round(pop.1yage.1ytime[,21:t.max],3) #want 2000-2100
  rubella.cases.1yage.1ytime <- round(rubella.cases.1yage.1ytime[,21:t.max],3) #want 2000-2100
  crs.1yage.1ytime <- round(crs.1yage.1ytime[,21:t.max],3) #want 2000-2100
  deaths.1yage.1ytime <- round(deaths.1yage.1ytime[,21:t.max],3) #want 2000-2100
  dalys.1yage.1ytime <- round(dalys.1yage.1ytime[,21:t.max],3) #want 2000-2100
  rubella.seropos.1yage.1ytime <- round(rubella.seropos.1yage.1ytime[,1:41],3) #want 1980-2020
  crs.rate <- round(crs.rate[1:t.max],3)

  return(list(pop.1yage.1ytime=pop.1yage.1ytime,
              rubella.cases.1yage.1ytime=rubella.cases.1yage.1ytime,
              crs.1yage.1ytime=crs.1yage.1ytime,
              deaths.1yage.1ytime=deaths.1yage.1ytime,
              dalys.1yage.1ytime=dalys.1yage.1ytime,
              rubella.seropos.1yage.1ytime=rubella.seropos.1yage.1ytime,
              crs.rate=crs.rate))
}

