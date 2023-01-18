
## Revised functions 
## 29 November 2022


## Revised Functions - These functions already exist but need to be updated ####
## This revised functions should not break what already exists

#Class definition for transition matrix. This matrix knows how to shift the
#ID.state.matrix from one time step to another, and contains all the fun stuff
#we need to do the transition (i.e., the matrix)
setClass("ID.transition.SIR",
         representation(n.epi.class="numeric", #number of different epi classes
                        epi.class="numeric", #the epi class of each row, between 1 and n.epi.classes
                        epi.class.label = "character", #epi class labels
                        s.inds = "numeric", #indexes of the susceptible states
                        i.inds = "numeric", #indexes of the infectious states
                        r.inds = "numeric", #indexes of the recovered states
                        n.age.class = "numeric", #the number of age classes
                        age.class = "numeric",#the actual age classes defined by upper age in class
                        aging.rate = "numeric", #percent that age out of each age class at each time step
                        survival.rate = "ANY", #the percent who survive in this age class at each time step -- xxamy changed to "ANY"
                        birth.rate = "ANY", #the birth rate for this transition -- xxamy changed to "ANY"
                        waifw = "matrix", #who acquires infection from whom matrix
                        age.surv.matrix = "matrix", #matrix governing age/survival
                        introduction.rate = "numeric", #an vector of possibly 0 case importation rates in each age class
                        exponent = "numeric",        #exponent on infecteds to correct for discretization
                        frequency.dep = "logical",	#true/false indicating if freq. dependence or density dep. is implemented
                        is.stochastic = "logical",   #true/false indicating if stochasticity desired
                        get.births = "function" #a function that takes state and birth rate; allowing births to be shaped by pop struct
         ))

#Class definition for transition matrix. This matrix knows how to shift the
#ID.state.matrix from one time step to another, and contains all the fun stuff
#we need to do the transition (i.e., the matrix)
setClass("nMx",
         representation(rate.years="numeric", #year associated with each time-specific rates --xxamy added this new object rate.years
                        rates="data.frame", #age specific death rates per 1 (rows) by time (columns)
                        mid.age="numeric" #mid-age associated with age-specific rates
         ))


#Function to pull in Age Specific Death Rates and return survivorship for running out transients
#
#Paramaters -
#   age.classes - a set of age classes for which tran is being built
#   generation time - generation time
#   nMx - object of age specific death rates over time and mid-age per rate
#   year - year to pull age specific death rates that will be used to run out the transients (should coincide with DFE year)
#
#Returns-
#   returns the age profile of survivorship in units of the generation time
create.surv.prob.over.age <- function(age.classes, generation.time, nMx=NULL, year=1990, check=F){
  
  rate.years <- nMx@rate.years #seq(1950,2100,5) #xxamy - revised what feeds into rate.years here
  index <- min(which(findInterval(rate.years, year)==1))
  rates.overtime <- nMx@rates
  mid.age <- nMx@mid.age
  generation.time.year <- generation.time/12   # put generation.time in terms of years rather than months (0.5 in months) now (0.041667 in years)
  
  yearly.mortality <- rates.overtime[,index]
  
  # Fit a smooth spline to this
  fit <- smooth.spline(mid.age, log(yearly.mortality), df=length(mid.age))
  
  # Convert the fit to our ages
  survs <- 1-exp(predict(fit, age.classes/12)$y) 
  # gives me only the $y or estiamted m.x.n for age.classes/12 (puts age.classes in years rather than months)
  # for each age predict the survival rate based on the smooth.spline
  # exp to get back to m.x.n from it's previous log state and take 1-m.x.n to get survival rate
  #plot(survs)
  
  # Convert yearly rates to generation time rates
  survs <- 1+(log(survs)/(1/generation.time.year)) 
  
  # Need to adjust last age class to deplete
  survs[length(survs)] <-  0.5^(12*1/generation.time)
  
  if (check) {
    plot(mid.age, yearly.mortality, type="b", xlim=c(0,100) ,ylim=range(yearly.mortality[mid.age<100]),pch=19)
    points(age.classes/12,(1-survs^(1/generation.time.year)), col=2,type="b")
    points(fit$x,exp(fit$y),col=4, type="l")
  }
  
  return(survs)
}


#Function to pull Age Specific Death Rates and return survivorship over time
#
#Parameters -
#   age.classes - a set of age classes for which tran is being built
#   generation.time - generation time
#   nMx - object of age specific death rates over time and mid-age per rate
#   nMx.years - years to interpolate age specific death rates 
#
#Returns-
#   returns the age profile (rows) of survivorship (in units of the generation time) by time (cols) as a matrix
create.surv.prob.over.age.time <- function(age.classes, generation.time, nMx=NULL, nMx.years, check=F){
  
  mid.age <- nMx@mid.age
  generation.time.year <- generation.time/12   # put generation.time in terms of years rather than months (0.5 in months) now (0.041667 in years)
  no.gens.in.year <- 12*1/generation.time
  
  rate.years <- nMx@rate.years #seq(1950,2100,5) #xxamy - revised what feeds into rate.years here
  indexes <- findInterval(nMx.years, rate.years)
  rates.years <- data.frame(nMx@rates[,indexes])
  rates.generations <- rates.years[,rep(seq_len(ncol(rates.years)), each=no.gens.in.year)]
  rates.generations <- cbind(rates.generations[,1], rates.generations) #adding additional column to match the length of the experiment time steps
  
  survs <- matrix(NA, nrow=length(age.classes), ncol=ncol(rates.generations))
  
  for (t in 1:ncol(rates.generations)){
    yearly.mortality <- rates.generations[,t]
    
    # Fit a smooth spline to this
    fit <- smooth.spline(mid.age, log(yearly.mortality), df=length(mid.age))
    
    # Convert the fit to our ages
    survs[,t] <- 1-exp(predict(fit, age.classes/12)$y) 
    # gives me only the $y or estiamted m.x.n for age.classes/12 (puts age.classes in years rather than months)
    # for each age predict the survival rate based on the smooth.spline
    # exp to get back to m.x.n from it's previous log state and take 1-m.x.n to get survival rate
    #plot(survs)
    
    # Convert yearly rates to generation time rates
    survs[,t] <- 1+(log(survs[,t])/(1/generation.time.year)) 
    
    # Need to adjust last age class to deplete
    survs[,t][length(survs[,t])] <-  0.5^no.gens.in.year
  }
  
  if (check) {
    plot(mid.age, yearly.mortality, type="b", xlim=c(0,50) ,ylim=range(yearly.mortality[mid.age<60]),pch=19)
    points(age.classes/12,(1-survs[,t]^(1/generation.time.year)), col=2,type="b")
    points(fit$x,exp(fit$y),col=4, type="l")
  }
  return(survs)
  
}


### Function to get population size, age structure, and fertility from 1950 to 2105 
### the last five years I just keep all 2100 rates constant
### based on unpd wpp2017 package estimates
###
##' @param - uncode - UN country code
##' @param - age.classes  - in years
##' @param - if.missing.use.region - boolean (MHL, XK, TUV)
###
### Returns demography for each country
### Use to be getDemography.original
getDemography.wpp2017 <- function(uncode=NA, age.classes=c(1:101), if.missing.use.region=F){
  
  time.by5 <- seq(1950,2100,5)
  time <- seq(1950,2100,1)
  
  #UNPD data
  data(pop)
  
  if (any(pop$country_code==uncode)){
    cc <- uncode 
  } else if (if.missing.use.region) {
    if (uncode == 999) {
      region <- "Europe"
      cc <- pop$country_code[tolower(pop$name)==tolower(region)]
      print(paste("No UNPD data for uncode ", uncode, "; using region ", region, " instead", sep=""))
    }
    if (uncode != 999) {
      region <- countrycode::countrycode(uncode, origin="un", destination="continent")
      cc <- pop$country_code[tolower(pop$name)==tolower(region)]
      print(paste("No UNPD data for uncode ", uncode, "; using region ", region, " instead", sep=""))
    }
  } else {
    print(paste("No UNPD data for uncode ", uncode, "; NAs produced", sep=""))
    pop.total.1950.2100=NA
    pop.age.byageclasses.1950.2100=NA
    tfr.1950.2100=NA
    births.1950.2100=NA
    cbr.1950.2100=NA
    e0.1950.2100=NA
    asfr.1950.2100=NA
    repro.age.sex.dist.1950.2100=NA
    asdr.1950.2100.by5=NA
    asdr.object=NA
  }
  
  #Population total in years 1950 to 2100
  data(pop)
  data(popproj)
  pop.total.1950.2100.by5 <-  as.numeric(cbind(pop[pop$country_code==cc,3:ncol(pop)],
                                               popproj[popproj$country_code==cc,3:ncol(popproj)]))
  f <- smooth.spline(time.by5, pop.total.1950.2100.by5)
  pop.total.1950.2100 <- predict(f,time)$y 
  
  #Population Age Structure by Age over time
  data(popFprojMed)
  data(popMprojMed)
  popF.age.1950.2100.by5 <- cbind(popF[popF$country_code==cc,4:ncol(popF)], popFprojMed[popFprojMed$country_code==cc,4:ncol(popFprojMed)])
  popM.age.1950.2100.by5 <- cbind(popM[popM$country_code==cc,4:ncol(popM)], popMprojMed[popMprojMed$country_code==cc,4:ncol(popMprojMed)])
  pop.age.1950.2100.by5 <- popF.age.1950.2100.by5+popM.age.1950.2100.by5
  ages.by5 <- popF[popF$country_code==cc,3]
  mid.age.by5 <- seq(2.5, 102.5, 5)
  popF.age.1950.2100 <- popM.age.1950.2100 <- matrix(NA, length(ages.by5), length(time))
  for (age in 1:length(ages.by5)){
    #new with wpp2017 because there are NA in older ages and oldest cohorts
    time.by5.index <- which(!is.na(popF.age.1950.2100.by5[age,]))
    time.by5.tmp <- time.by5[time.by5.index] 
    time.index <- ((time.by5.index[1]-1)*5+1):length(time)
    time.tmp <- time[time.index]
    f <- smooth.spline(time.by5.tmp, popF.age.1950.2100.by5[age,time.by5.index])
    popF.age.1950.2100[age,time.index] <- predict(f,time.tmp)$y 
    m <- smooth.spline(time.by5.tmp, popM.age.1950.2100.by5[age,time.by5.index])
    popM.age.1950.2100[age,time.index] <- predict(m,time.tmp)$y 
  }
  colnames(popF.age.1950.2100) <- colnames(popM.age.1950.2100) <- time
  rownames(popF.age.1950.2100) <- rownames(popM.age.1950.2100) <- ages.by5
  #forcing 100+ to be 0 if NA - new with wpp2017 - then using na.spline in the zoo package to fill them in
  popF.age.1950.2100[21,is.na(popF.age.1950.2100[21,])] <- 0
  popM.age.1950.2100[21,is.na(popM.age.1950.2100[21,])] <- 0
  popF.age.byageclasses.1950.2100 <- apply(popF.age.1950.2100, 2, function(x) predict(smooth.spline(mid.age.by5, na.spline(x)), age.classes)$y)
  popM.age.byageclasses.1950.2100 <- apply(popM.age.1950.2100, 2, function(x) predict(smooth.spline(mid.age.by5, na.spline(x)), age.classes)$y)
  #adjust for varying bin width!
  popF.age.byageclasses.1950.2100 <- popF.age.byageclasses.1950.2100*diff(c(0,age.classes)) 
  popM.age.byageclasses.1950.2100 <- popM.age.byageclasses.1950.2100*diff(c(0,age.classes))
  pop.tmp <- popF.age.byageclasses.1950.2100+popM.age.byageclasses.1950.2100
  #rescale the age distibution to the estimated population each year (pop.total.1950.2100)
  pop.age.byageclasses.1950.2100 <- matrix(NA, nrow(pop.tmp), ncol(pop.tmp))
  for (t in 1:ncol(pop.tmp)){
    pop.age.byageclasses.1950.2100[,t] <- pop.total.1950.2100[t]*pop.tmp[,t]/sum(pop.tmp[,t])
  }
  rownames(pop.age.byageclasses.1950.2100) <- age.classes
  colnames(pop.age.byageclasses.1950.2100) <- time
  
  #Sex Distribution by age and over time
  pop.tmp <- popF.age.byageclasses.1950.2100+popM.age.byageclasses.1950.2100
  sex.dist <- matrix(NA, nrow(pop.tmp), ncol(pop.tmp))
  for (c in 1:ncol(pop.tmp)){
    sex.dist[,c] <- popF.age.byageclasses.1950.2100[,c]/pop.tmp[,c]
  }
  #Sex Distribution grouped by 5 year reproductive ages over time
  repro.ages.index <- findInterval(age.classes, seq(15,50,5))+1
  repro.age.sex.dist.1950.2100 <- matrix(NA, nrow=length(unique(repro.ages.index)), ncol=ncol(sex.dist))
  for (i in 1:length(unique(repro.ages.index))){
    repro.age.sex.dist.1950.2100[i,] <- apply(sex.dist, 2, function(x) mean(x[repro.ages.index==i]))
  }
  repro.age.sex.dist.1950.2100 <- repro.age.sex.dist.1950.2100[-c(1,9),]
  rownames(repro.age.sex.dist.1950.2100) <- c("15-19","20-24","25-29","30-34","35-39","40-44","45-49")
  colnames(repro.age.sex.dist.1950.2100) <- time
  
  #TFR over time
  data(tfr)
  data(tfrprojMed)
  mid.time.by5 <- seq(1952, 2097, 5) #TFR is given over a range, therefore assume it is the mid-period
  tfr.1950.2100.by5 <-  as.numeric(cbind(tfr[tfr$country_code==cc,3:15],
                                         tfrprojMed[tfrprojMed$country_code==cc,3:ncol(tfrprojMed)]))
  f <- smooth.spline(mid.time.by5, tfr.1950.2100.by5)
  tfr.1950.2100 <- predict(f,time)$y 
  names(tfr.1950.2100) <- time
  
  #Number of women of reproductive age (in five year age groups 15 to 50) over time
  repro.ages <- seq(15,45,5)
  popF.15to50.1950.2100 <- matrix(NA,7,length(1950:2100))
  for (c in 1:ncol(pop.age.byageclasses.1950.2100)){
    for (a in 1:length(repro.ages)){
      index <- which(age.classes>=repro.ages[a] & age.classes<(repro.ages[a]+5))
      popF.15to50.1950.2100[a,c] <- sum(pop.age.byageclasses.1950.2100[index,c]*sex.dist[index,c])
    }
  }
  colnames(popF.15to50.1950.2100) <- time
  rownames(popF.15to50.1950.2100) <- repro.ages
  
  #ASFR overtime
  data(percentASFR)
  p.asfr <- (percentASFR[percentASFR$country_code==cc,])
  p.asfr.1950.2100 <- matrix(NA, nrow=nrow(p.asfr), ncol=(length(1950:2100)))
  for (t in 1:(ncol(p.asfr)-3)){
    p.asfr.1950.2100[,(t*5-4):(t*5)] <- p.asfr[,(t+3)]/100
  }
  p.asfr.1950.2100[,ncol(p.asfr.1950.2100)] <-  p.asfr[,ncol(p.asfr)]/100
  asfr.1950.2100 <- matrix(NA, nrow=nrow(p.asfr), ncol=(length(1950:2100)))
  for (t in 1:ncol(asfr.1950.2100)){
    asfr.1950.2100[,t] <- tfr.1950.2100[t]/5*p.asfr.1950.2100[,t]
  }
  colnames(asfr.1950.2100) <- time
  
  #Births over time
  births.1950.2100 <- rep(NA, length(1950:2100))
  for (t in 1:length(births.1950.2100)){
    births.1950.2100[t] <- sum(asfr.1950.2100[,t]*popF.15to50.1950.2100[,t])
  }
  names(births.1950.2100) <- time
  
  #Crude Birth Rates over time
  cbr.1950.2100 <- rep(NA, length(1950:2100))
  cbr.1950.2100 <- births.1950.2100/pop.total.1950.2100*1000
  names(cbr.1950.2100) <- time
  
  #Age Specific Death Rate over time
  data(mxF)
  asdr.1950.2100.by5 <- (mxF[mxF$country_code==cc,4:ncol(mxF)])
  asdr.1950.2100.by5 <- cbind(asdr.1950.2100.by5, "2100-2105" = asdr.1950.2100.by5[,ncol(asdr.1950.2100.by5)]) #revised to 2105 xxamy
  asdr.maxage <- 5*(length((mxF[mxF$country_code==cc,3]))-2)
  asdr.object <- new("nMx",
                     rate.years = seq(1950, 2095, 5), #xxamy - added rate.years 
                     rates = asdr.1950.2100.by5,
                     mid.age = c(0.5, seq(2.5, (asdr.maxage+2.5), 5)))
  
  #Life Expectancy over time
  data(e0F)
  data(e0Fproj)
  mid.time.by5 <- seq(1952, 2097, 5) #e0 is given over a range, therefore assume it is the mid-period
  e0.1950.2100.by5 <-  as.numeric(cbind(e0F[e0F$country_code==cc,3:15],
                                        e0Fproj[e0Fproj$country_code==cc,3:ncol(e0Fproj)]))
  f <- smooth.spline(mid.time.by5, e0.1950.2100.by5)
  e0.1950.2100 <- predict(f,time)$y 
  names(e0.1950.2100) <- time
  
  
  return(list(pop.total.1950.2100=pop.total.1950.2100, 
              pop.age.byageclasses.1950.2100=pop.age.byageclasses.1950.2100, 
              tfr.1950.2100=tfr.1950.2100,
              e0.1950.2100=e0.1950.2100,
              asfr.1950.2100=asfr.1950.2100,
              repro.age.sex.dist.1950.2100=repro.age.sex.dist.1950.2100,
              births.1950.2100=births.1950.2100, 
              cbr.1950.2100=cbr.1950.2100,
              asdr.1950.2100.by5=asdr.1950.2100.by5,
              asdr.object=asdr.object)) 
}


### Function to get population size, age structure, and fertility from 1950 to 2105 
### the last five years I just keep all 2100 rates constant
### based on unpd wpp2019 package estimates
###
##' @param - uncode - UN country code
##' @param - age.classes  - in years
##' @param - if.missing.use.region - boolean (MHL, XK, TUV)
###
### Returns demography for each country
getDemography.wpp2019 <- function(uncode=NA, age.classes=c(1:101), if.missing.use.region=F){
  
  cc <- uncode
  
  time.by5 <- seq(1950,2100,5)
  time <- seq(1950,2100,1)
  
  #UNPD data
  data(pop)
  
  if (any(pop$country_code==uncode)){
    cc <- uncode 
  } else if (if.missing.use.region) {
    if (uncode == 999) {
      region <- "Europe"
      cc <- pop$country_code[tolower(pop$name)==tolower(region)]
      print(paste("No UNPD data for uncode ", uncode, "; using region ", region, " instead", sep=""))
    }
    if (uncode != 999) {
      region <- countrycode::countrycode(uncode, origin="un", destination="continent")
      cc <- pop$country_code[tolower(pop$name)==tolower(region)]
      print(paste("No UNPD data for uncode ", uncode, "; using region ", region, " instead", sep=""))
    }
  } else {
    print(paste("No UNPD data for uncode ", uncode, "; NAs produced", sep=""))
    pop.total.1950.2100=NA
    pop.age.byageclasses.1950.2100=NA
    tfr.1950.2100=NA
    births.1950.2100=NA
    cbr.1950.2100=NA
    e0.1950.2100=NA
    asfr.1950.2100=NA
    repro.age.sex.dist.1950.2100=NA
    asdr.1950.2100.by5=NA
    asdr.object=NA
  }
  
  #Population total in years 1950 to 2100
  data(pop)
  data(popproj)
  
  pop.total.1950.2100.by5 <-  as.numeric(cbind(pop[pop$country_code==cc,3:ncol(pop)],
                                               popproj[popproj$country_code==cc,3:ncol(popproj)]))
  f <- smooth.spline(time.by5, pop.total.1950.2100.by5)
  pop.total.1950.2100 <- predict(f,time)$y 
  
  #Population Age Structure by Age over time
  data(popFprojMed)
  data(popMprojMed)
  popF.age.1950.2100.by5 <- cbind(popF[popF$country_code==cc,4:ncol(popF)], popFprojMed[popFprojMed$country_code==cc,4:ncol(popFprojMed)])
  popM.age.1950.2100.by5 <- cbind(popM[popM$country_code==cc,4:ncol(popM)], popMprojMed[popMprojMed$country_code==cc,4:ncol(popMprojMed)])
  pop.age.1950.2100.by5 <- popF.age.1950.2100.by5+popM.age.1950.2100.by5
  ages.by5 <- popF[popF$country_code==cc,3]
  mid.age.by5 <- seq(2.5, 102.5, 5)
  popF.age.1950.2100 <- popM.age.1950.2100 <- matrix(NA, length(ages.by5), length(time))
  for (age in 1:length(ages.by5)){
    #new with wpp2017 because there are NA in older ages and oldest cohorts
    time.by5.index <- which(!is.na(popF.age.1950.2100.by5[age,]))
    time.by5.tmp <- time.by5[time.by5.index] 
    time.index <- ((time.by5.index[1]-1)*5+1):length(time)
    time.tmp <- time[time.index]
    f <- smooth.spline(time.by5.tmp, popF.age.1950.2100.by5[age,time.by5.index])
    popF.age.1950.2100[age,time.index] <- predict(f,time.tmp)$y 
    m <- smooth.spline(time.by5.tmp, popM.age.1950.2100.by5[age,time.by5.index])
    popM.age.1950.2100[age,time.index] <- predict(m,time.tmp)$y 
  }
  colnames(popF.age.1950.2100) <- colnames(popM.age.1950.2100) <- time
  rownames(popF.age.1950.2100) <- rownames(popM.age.1950.2100) <- ages.by5
  #forcing 100+ to be 0 if NA - new with wpp2017 - then using na.spline in the zoo package to fill them in
  popF.age.1950.2100[21,is.na(popF.age.1950.2100[21,])] <- 0
  popM.age.1950.2100[21,is.na(popM.age.1950.2100[21,])] <- 0
  popF.age.byageclasses.1950.2100 <- apply(popF.age.1950.2100, 2, function(x) predict(smooth.spline(mid.age.by5, na.spline(x)), age.classes)$y)
  popM.age.byageclasses.1950.2100 <- apply(popM.age.1950.2100, 2, function(x) predict(smooth.spline(mid.age.by5, na.spline(x)), age.classes)$y)
  #new wpp2019 - if prediction creates negative numbers than force to 0
  popF.age.byageclasses.1950.2100[popF.age.byageclasses.1950.2100<0] <- 0
  popM.age.byageclasses.1950.2100[popM.age.byageclasses.1950.2100<0] <- 0
  #adjust for varying bin width!
  popF.age.byageclasses.1950.2100 <- popF.age.byageclasses.1950.2100*diff(c(0,age.classes)) 
  popM.age.byageclasses.1950.2100 <- popM.age.byageclasses.1950.2100*diff(c(0,age.classes))
  pop.tmp <- popF.age.byageclasses.1950.2100+popM.age.byageclasses.1950.2100
  #rescale the age distibution to the estimated population each year (pop.total.1950.2100)
  pop.age.byageclasses.1950.2100 <- matrix(NA, nrow(pop.tmp), ncol(pop.tmp))
  for (t in 1:ncol(pop.tmp)){
    pop.age.byageclasses.1950.2100[,t] <- pop.total.1950.2100[t]*pop.tmp[,t]/sum(pop.tmp[,t])
  }
  rownames(pop.age.byageclasses.1950.2100) <- age.classes
  colnames(pop.age.byageclasses.1950.2100) <- time
  
  #Sex Distribution by age and over time
  pop.tmp <- popF.age.byageclasses.1950.2100+popM.age.byageclasses.1950.2100
  sex.dist <- matrix(NA, nrow(pop.tmp), ncol(pop.tmp))
  for (c in 1:ncol(pop.tmp)){
    sex.dist[,c] <- popF.age.byageclasses.1950.2100[,c]/pop.tmp[,c]
  }
  #to correct for NAs when dividing by pop.tmp=0 - doesnt really matter because no populaton in those cells though
  sex.dist[is.na(sex.dist)] <- 0.5
  
  #Sex Distribution grouped by 5 year reproductive ages over time
  repro.ages.index <- findInterval(age.classes, seq(15,50,5))+1
  repro.age.sex.dist.1950.2100 <- matrix(NA, nrow=length(unique(repro.ages.index)), ncol=ncol(sex.dist))
  for (i in 1:length(unique(repro.ages.index))){
    repro.age.sex.dist.1950.2100[i,] <- apply(sex.dist, 2, function(x) mean(x[repro.ages.index==i]))
  }
  repro.age.sex.dist.1950.2100 <- repro.age.sex.dist.1950.2100[-c(1,9),]
  rownames(repro.age.sex.dist.1950.2100) <- c("15-19","20-24","25-29","30-34","35-39","40-44","45-49")
  colnames(repro.age.sex.dist.1950.2100) <- time
  
  #TFR over time
  data(tfr)
  data(tfrprojMed)
  mid.time.by5 <- seq(1952, 2097, 5) #TFR is given over a range, therefore assume it is the mid-period
  tfr.1950.2100.by5 <-  as.numeric(cbind(tfr[tfr$country_code==cc,3:16],
                                         tfrprojMed[tfrprojMed$country_code==cc,3:ncol(tfrprojMed)]))
  f <- smooth.spline(mid.time.by5, tfr.1950.2100.by5)
  tfr.1950.2100 <- predict(f,time)$y 
  names(tfr.1950.2100) <- time
  
  #Number of women of reproductive age (in five year age groups 15 to 50) over time
  repro.ages <- seq(15,45,5)
  popF.15to50.1950.2100 <- matrix(NA,7,length(1950:2100))
  for (c in 1:ncol(pop.age.byageclasses.1950.2100)){
    for (a in 1:length(repro.ages)){
      index <- which(age.classes>=repro.ages[a] & age.classes<(repro.ages[a]+5))
      popF.15to50.1950.2100[a,c] <- sum(pop.age.byageclasses.1950.2100[index,c]*sex.dist[index,c])
    }
  }
  colnames(popF.15to50.1950.2100) <- time
  rownames(popF.15to50.1950.2100) <- repro.ages
  
  #ASFR overtime
  data(percentASFR)
  p.asfr <- (percentASFR[percentASFR$country_code==cc,])
  p.asfr.1950.2100 <- matrix(NA, nrow=nrow(p.asfr), ncol=(length(1950:2100)))
  for (t in 1:(ncol(p.asfr)-3)){
    p.asfr.1950.2100[,(t*5-4):(t*5)] <- p.asfr[,(t+3)]/100
  }
  p.asfr.1950.2100[,ncol(p.asfr.1950.2100)] <-  p.asfr[,ncol(p.asfr)]/100
  asfr.1950.2100 <- matrix(NA, nrow=nrow(p.asfr), ncol=(length(1950:2100)))
  for (t in 1:ncol(asfr.1950.2100)){
    asfr.1950.2100[,t] <- tfr.1950.2100[t]/5*p.asfr.1950.2100[,t]
  }
  colnames(asfr.1950.2100) <- time
  
  #Births over time
  births.1950.2100 <- rep(NA, length(1950:2100))
  for (t in 1:length(births.1950.2100)){
    births.1950.2100[t] <- sum(asfr.1950.2100[,t]*popF.15to50.1950.2100[,t])
  }
  names(births.1950.2100) <- time
  
  #Crude Birth Rates over time
  cbr.1950.2100 <- rep(NA, length(1950:2100))
  cbr.1950.2100 <- births.1950.2100/pop.total.1950.2100*1000
  names(cbr.1950.2100) <- time
  
  #Age Specific Death Rate over time
  data(mxF)
  asdr.1950.2100.by5 <- (mxF[mxF$country_code==uncode,4:ncol(mxF)])
  asdr.1950.2100.by5 <- cbind(asdr.1950.2100.by5, "2100-2105" = asdr.1950.2100.by5[,ncol(asdr.1950.2100.by5)]) #revised to 2105 xxamy
  asdr.maxage <- 5*(length((mxF[mxF$country_code==uncode,3]))-2)
  asdr.object <- new("nMx",
                     rate.years = seq(1950, 2095, 5), #xxamy - added rate.years 
                     rates = asdr.1950.2100.by5,
                     mid.age = c(0.5, seq(2.5, (asdr.maxage+2.5), 5)))
  
  #Life Expectancy over time
  data(e0F)
  data(e0Fproj)
  mid.time.by5 <- seq(1952, 2097, 5) #e0 is given over a range, therefore assume it is the mid-period
  e0.1950.2100.by5 <-  as.numeric(cbind(e0F[e0F$country_code==cc,3:16],
                                        e0Fproj[e0Fproj$country_code==cc,3:ncol(e0Fproj)]))
  f <- smooth.spline(mid.time.by5, e0.1950.2100.by5)
  e0.1950.2100 <- predict(f,time)$y 
  names(e0.1950.2100) <- time
  
  return(list(pop.total.1950.2100=pop.total.1950.2100, 
              pop.age.byageclasses.1950.2100=pop.age.byageclasses.1950.2100, 
              tfr.1950.2100=tfr.1950.2100,
              e0.1950.2100=e0.1950.2100,
              asfr.1950.2100=asfr.1950.2100,
              repro.age.sex.dist.1950.2100=repro.age.sex.dist.1950.2100,
              births.1950.2100=births.1950.2100, 
              cbr.1950.2100=cbr.1950.2100,
              asdr.1950.2100.by5=asdr.1950.2100.by5,
              asdr.object=asdr.object)) 
}


#Class to define an experiment.
setClass("experiment",
         representation(name = "character", #a human readable name for the experiment
                        R0 = "numeric", #if included the R0 we scale to
                        state.t0 = "ID.state.matrix",#the initial state
                        t.min = "numeric", #the starting time for the experiment in years
                        t.max = "numeric", #the ending time for the experiment in years
                        t0.doy = "numeric", #the starting day of year, should default to 0
                        step.size = "numeric", #the step size in years, 1/no.gens.per.year
                        trans = "ID.transition.SIR", #defines the transitions for this experiment
                        season.obj = "seasonal.force", #seasonal forcing function
                        births.per.1000.each.timestep = "ANY", #the estimated crude birth rate per time step, length is T -- xxamy changed to "ANY"
                        description = "character" #describe the experiment
         ))

#Class for an experiment with birth and death rates updated each time step and population rescaling
setClass("experiment.updatedemog",
         representation(trans = "ID.transition.SIR.vac", #the transitions for this experiment
                        surv.each.timestep = "matrix", #the estimated age-specific survival rates over the length of the experiment, ncol is T
                        pop.rescale.each.timestep = "ANY", #a vector of 0's if not call for population rescale, if !=0 then rescale by this pop number --xxamy changed to "ANY"
                        maternal.obj = "maternal.exp.decay"),  #maternal antibody object
         contains=c("experiment")) 

#Class for an experiment with changing vaccination coverage over time
setClass("experiment.updatedemog.vaccinationchange",
         representation(time.specific.MR1cov = "ANY", #vector of values of the coverage values --xxamy changed to "ANY"
                        time.specific.MR2cov = "ANY", #vector of values of the coverage values --xxamy changed to "ANY"
                        time.specific.SIAcov = "ANY", #vector of values of the coverage values --xxamy changed to "ANY"
                        time.specific.min.age.MR1 = "ANY", #--xxamy changed to "ANY"
                        time.specific.max.age.MR1 = "ANY", #--xxamy changed to "ANY"
                        time.specific.min.age.MR2 = "ANY", #--xxamy changed to "ANY"
                        time.specific.max.age.MR2 = "ANY", #--xxamy changed to "ANY"
                        time.specific.min.age.SIA = "ANY", #--xxamy changed to "ANY"
                        time.specific.max.age.SIA = "ANY", #--xxamy changed to "ANY"
                        obj.vcdf.MR1 = "vaccine.cdf.byage",
                        obj.vcdf.MR2 = "vaccine.cdf.byage",
                        obj.prob.vsucc = "prob.vsucc.byage",
                        sia.timing.in.year = "ANY"), #time of the SIA in years (e.g., 0.5 is July of the year) --xxamy changed to "ANY"
         contains=c("experiment.updatedemog"))

#Class for an experiment 
setClass("experiment.updatedemog.vaccinationchange.vaccinationlimitations",
         representation(MR1MR2correlation = "logical", 
                        MR1SIAcorrelation = "logical", 
                        MR2SIAcorrelation = "logical", 
                        SIAinefficient = "logical",
                        SIAinacc = "logical",
                        prop.inacc = "ANY"), #--xxamy changed to "ANY"
         contains=c("experiment.updatedemog.vaccinationchange"))

# Holds the results of a simulation with an experiment.updatedemog object
setClass("sim.results.MSIRV.update.demog",
         representation(births.each.timestep = "ANY", #the number of births per time step as output from simulation --xxamy changed to "ANY"
                        growth.rate.each.timestep = "ANY" #--xxamy changed to "ANY"
         ),
         contains="sim.results.MSIRV")

# Holds the results output for a simulation with an experiment.updatedemog.vaccinationchange object
setClass("sim.results.MSIRV.update.demog.vaccine.change",
         representation(MR1.fail.each.timestep = "ANY", #the number of births per time step as output from simulation --xxamy changed to "ANY"
                        MR2.fail.each.timestep = "ANY", #--xxamy changed to "ANY"
                        SIA.fail.each.timestep = "ANY" #--xxamy changed to "ANY"
         ),
         contains="sim.results.MSIRV.update.demog")

setMethod("run",
          "experiment.updatedemog.vaccinationchange",
          function(exper, rescale.WAIFW=T, ...) {
            
            #print(rescale.WAIFW)
            state <- exper@state.t0
            
            #rescale the WAIFW if a specific R0 is specified
            if (rescale.WAIFW & length(exper@R0)>0) {
              #print("RESCALING!!!")
              exper@trans@waifw <- scale.WAIFW(exper@R0,
                                               state,exper@trans@waifw,
                                               frequency.dep=exper@trans@frequency.dep,
                                               suscept.state=exper@trans@s.inds[1])
            }
            
            #get the number of time steps in the experiment
            T <- round((exper@t.max-exper@t.min)/exper@step.size)+1
            
            if (!is.null(exper@season.obj)) {
              mults <- get.seasonal.mult(exper@t0.doy/365+(1:T-1)*exper@step.size,
                                         exper@season.obj)
            } else {
              mults <- rep(1,T)
            }
            
            #make a temporary transmission object
            tmp.trans <- exper@trans
            
            #hold the states as we walk through
            rc <- matrix(ncol = T, nrow = nrow(state))
            rc[,1] <- state
            
            #need output vector for births over time
            births.each.timestep <- growth.rate.each.timestep <- rep(NA, T)
            births.each.timestep[1] <- tmp.trans@birth.rate
            
            #generate the age and time specific vaccination matrix
            routine <- get.routine.time.age.specific(time.step= exper@step.size*12,
                                                     age.classes=exper@trans@age.class,
                                                     time.specific.MR1cov=exper@time.specific.MR1cov,
                                                     age.min.MR1=exper@time.specific.min.age.MR1, 
                                                     age.max.MR1=exper@time.specific.max.age.MR1, 
                                                     time.specific.MR2cov=exper@time.specific.MR2cov,
                                                     age.min.MR2=exper@time.specific.min.age.MR2, 
                                                     age.max.MR2=exper@time.specific.max.age.MR2, 
                                                     obj.vcdf.MR1=exper@obj.vcdf.MR1,
                                                     obj.vcdf.MR2=exper@obj.vcdf.MR2,
                                                     obj.prob.vsucc=exper@obj.prob.vsucc,
                                                     MR1MR2correlation=F)
            routine.intro <- rep(0, T)
            if (any(exper@time.specific.MR1cov!=0)) routine.intro[min(which(exper@time.specific.MR1cov>0))*(1/exper@step.size)+1] <- 1
            if (any(exper@time.specific.MR2cov!=0)) routine.intro[min(which(exper@time.specific.MR2cov>0))*(1/exper@step.size)+1] <- 1
            index.routine.vacc <- c(1,rep(1:nrow(routine$age.time.specific.routine), each=(T-1)/exper@t.max))
            SIA <- get.sia.time.age.specific(age.classes=exper@trans@age.class,
                                             time.specific.SIAcov=exper@time.specific.SIAcov,
                                             age.min.sia=exper@time.specific.min.age.SIA,
                                             age.max.sia=exper@time.specific.max.age.SIA,
                                             obj.prob.vsucc=exper@obj.prob.vsucc)
            index.sia.vacc <- rep(NA,T)
            year.sia <- which(exper@time.specific.SIAcov!=0)
            index.sia.vacc[(year.sia-1)*(T-1)/exper@t.max + round((exper@sia.timing.in.year*(T-1)/exper@t.max))] <-  year.sia #minus 1 because adding the sia.timing
            sia.times <- ifelse(!is.na(index.sia.vacc), 1, 0)
            
            #need output vectors for primary vaccination failure over time
            MR1.fail.each.timestep <- MR2.fail.each.timestep <- SIA.fail.each.timestep <- rep(0, T)
            MR1.fail.each.timestep[1] <- routine$prop.fail.MR1[1]
            MR2.fail.each.timestep[1] <- routine$prop.fail.MR2[1]
            SIA.fail.each.timestep[1] <- 0 
            
            for (t in 2:T) {
              
              #scale the waifw by seasonality
              if (!is.array(mults)){
                tmp.trans@waifw <- exper@trans@waifw*mults[t] 
              } else {
                tmp.trans@waifw <- exper@trans@waifw*mults[,,t]
              }
              
              #if exper@pop.rescale.each.timestep[t]==NaN then tmp.trans@pop.rescale==NaN and will not rescale
              #otherwise if any number other than NaN it will rescale
              #or if exper@pop.rescale.each.timestep not completed , then numeric(0) and any index [t] is NA
              if (!is.na(exper@pop.rescale.each.timestep[t])) {state <- exper@pop.rescale.each.timestep[t] *(state/sum(state))}
              
              #put in the correct birth rate for that time-step, if it varies
              if (length(exper@births.per.1000.each.timestep)>1) {
                # Read 'birth rate'
                tmp.trans@birth.rate = (exper@births.per.1000.each.timestep[t]*sum(state)/1000)
              } else {
                tmp.trans@birth.rate <- exper@trans@birth.rate*sum(rc[,t-1])/sum(exper@state.t0)
              }
              births.each.timestep[t] <-  tmp.trans@birth.rate
              
              #put in time appropriate survival rate otherwise it uses original surv.matrix set up for trans object and keeps constant over time
              if (!is.na(exper@surv.each.timestep[1,1])) {
                tmp.trans@age.surv.matrix <- ExtractAgeSpecificSurvivalMatrix(tmp.trans=tmp.trans, maternal.obj=exper@maternal.obj, 
                                                                              surv.at.timestep.t=exper@surv.each.timestep[,t])
              }
              
              ##put in correct vaccination coverage
              routine.vacc.prob <- routine$age.time.specific.routine
              sia.vacc.prob <- SIA$age.time.specific.SIA
              if (!is.na(index.sia.vacc[t])){ #if SIA
                tmp.trans@vac.per@pvacc.in.age.class <- 
                  routine.vacc.prob[index.routine.vacc[t],] + 
                  sia.vacc.prob[index.sia.vacc[t],] - #xxamy changed from addition to subtraction
                  (routine.vacc.prob[index.routine.vacc[t],]*sia.vacc.prob[index.sia.vacc[t],]) 
                #stow output
                MR1.fail.each.timestep[t] <- routine$prop.fail.MR1[index.routine.vacc[t]]
                MR2.fail.each.timestep[t] <- routine$prop.fail.MR2[index.routine.vacc[t]]
                SIA.fail.each.timestep[t] <- SIA$prop.fail.SIA[index.sia.vacc[t]]
                
              } else { #if no SIA
                tmp.trans@vac.per@pvacc.in.age.class <- routine.vacc.prob[index.routine.vacc[t],]
                #stow output
                MR1.fail.each.timestep[t] <- routine$prop.fail.MR1[index.routine.vacc[t]]
                MR2.fail.each.timestep[t] <- routine$prop.fail.MR2[index.routine.vacc[t]]
              }
              
              
              N0 <- sum(state)
              state <- next.ID.state(state, tmp.trans)
              NT <- sum(state)
              growth.rate.each.timestep[t] <- log(NT/N0) #instantaneous biweekly growth rate
              
              #print(t)
              #print(dim(state))
              rc [,t] <- state
              
            }
            
            
            rc <- new("sim.results.MSIRV.update.demog.vaccine.change",
                      data=rc,
                      m.inds = exper@trans@m.inds,
                      s.inds = exper@trans@s.inds,
                      i.inds = exper@trans@i.inds,
                      r.inds = exper@trans@r.inds,
                      v.inds = exper@trans@v.inds,
                      t = exper@t.min+(1:T-1)* exper@step.size,
                      age.class = exper@trans@age.class,
                      births.each.timestep = births.each.timestep,
                      growth.rate.each.timestep = growth.rate.each.timestep,
                      MR1.fail.each.timestep = MR1.fail.each.timestep,
                      MR2.fail.each.timestep = MR2.fail.each.timestep, 
                      SIA.fail.each.timestep = SIA.fail.each.timestep,
                      routine.intro = routine.intro,
                      sia.times = sia.times)
            
            rc <- new("experiment.result",
                      experiment.def = exper,
                      result = rc)
            
            
            return(rc)
          }
)

setMethod("run",
          "experiment.updatedemog.vaccinationchange.vaccinationlimitations",
          function(exper, rescale.WAIFW=T, ...) {
            
            #print(rescale.WAIFW)
            state <- exper@state.t0
            
            #rescale the WAIFW if a specific R0 is specified
            if (rescale.WAIFW & length(exper@R0)>0) {
              #print("RESCALING!!!")
              exper@trans@waifw <- scale.WAIFW(exper@R0,
                                               state,exper@trans@waifw,
                                               frequency.dep=exper@trans@frequency.dep,
                                               suscept.state=exper@trans@s.inds[1])
            }
            
            #get the number of time steps in the experiment
            T <- round((exper@t.max-exper@t.min)/exper@step.size)+1
            
            if (!is.null(exper@season.obj)) {
              mults <- get.seasonal.mult(exper@t0.doy/365+(1:T-1)*exper@step.size,
                                         exper@season.obj)
            } else {
              mults <- rep(1,T)
            }
            
            #make a temporary transmission object
            tmp.trans <- exper@trans
            
            #hold the states as we walk through
            rc <- matrix(ncol = T, nrow = nrow(state))
            rc[,1] <- state
            
            #need output vector for births over time
            births.each.timestep <- growth.rate.each.timestep <- rep(NA, T)
            births.each.timestep[1] <- tmp.trans@birth.rate
            
            #there are a couple options that are not yet coded - stop experiment if these are selected
            if ((exper@MR1SIAcorrelation & exper@SIAinacc) | (exper@MR2SIAcorrelation & exper@SIAinacc))  {
              stop("error: SIA limitation can either be correlation with a routine dose OR inaccessible population (with or without inefficiency), not both")
            }
            if ((exper@MR1SIAcorrelation & exper@SIAinefficient) | (exper@MR2SIAcorrelation & exper@SIAinefficient))  {
              stop("error: SIA limitation can either be correlation with a routine dose OR SIA inefficiency population (with or without inaccessible), not both")
            }
            if (!exper@MR1SIAcorrelation & exper@MR2SIAcorrelation){
              stop("error: MR1SIAcorrelation = FALSE and MR2SIAcorrelation= T option isn't yet coded or available")
            }
            
            #generate the age and time specific vaccination matrices for routine and SIAs
            routine <- get.routine.time.age.specific(time.step= exper@step.size*12,
                                                     age.classes=exper@trans@age.class,
                                                     time.specific.MR1cov=exper@time.specific.MR1cov,
                                                     age.min.MR1=exper@time.specific.min.age.MR1, 
                                                     age.max.MR1=exper@time.specific.max.age.MR1, 
                                                     time.specific.MR2cov=exper@time.specific.MR2cov,
                                                     age.min.MR2=exper@time.specific.min.age.MR2, 
                                                     age.max.MR2=exper@time.specific.max.age.MR2, 
                                                     obj.vcdf.MR1=exper@obj.vcdf.MR1,
                                                     obj.vcdf.MR2=exper@obj.vcdf.MR2,
                                                     obj.prob.vsucc=exper@obj.prob.vsucc,
                                                     MR1MR2correlation=exper@MR1MR2correlation)
            routine.intro <- rep(0, T)
            if (any(exper@time.specific.MR1cov!=0)) routine.intro[min(which(exper@time.specific.MR1cov>0))*(1/exper@step.size)+1] <- 1
            if (any(exper@time.specific.MR2cov!=0)) routine.intro[min(which(exper@time.specific.MR2cov>0))*(1/exper@step.size)+1] <- 1
            index.routine.vacc <- c(1,rep(1:nrow(routine$age.time.specific.routine), each=(T-1)/exper@t.max))
            SIA <- get.sia.time.age.specific(age.classes=exper@trans@age.class,
                                             time.specific.SIAcov=exper@time.specific.SIAcov,
                                             age.min.sia=exper@time.specific.min.age.SIA,
                                             age.max.sia=exper@time.specific.max.age.SIA,
                                             obj.prob.vsucc=exper@obj.prob.vsucc)
            index.sia.vacc <- rep(NA,T)
            year.sia <- which(exper@time.specific.SIAcov!=0)
            index.sia.vacc[(year.sia-1)*(T-1)/exper@t.max + round((exper@sia.timing.in.year*(T-1)/exper@t.max))] <-  year.sia #minus 1 because adding the sia.timing
            sia.times <- ifelse(!is.na(index.sia.vacc), 1, 0)
            
            #need output vectors for primary vaccination failure over time
            MR1.fail.each.timestep <- MR2.fail.each.timestep <- SIA.fail.each.timestep <- rep(0, T)
            MR1.fail.each.timestep[1] <- routine$prop.fail.MR1[1]
            MR2.fail.each.timestep[1] <- routine$prop.fail.MR2[1]
            SIA.fail.each.timestep[1] <- 0 
            
            for (t in 2:T) {
              
              #scale the waifw by seasonality
              if (!is.array(mults)){
                tmp.trans@waifw <- exper@trans@waifw*mults[t] 
              } else {
                tmp.trans@waifw <- exper@trans@waifw*mults[,,t]
              }
              
              #if exper@pop.rescale.each.timestep[t]==NaN then tmp.trans@pop.rescale==NaN and will not rescale
              #otherwise if any number other than NaN it will rescale
              #or if exper@pop.rescale.each.timestep not completed , then numeric(0) and any index [t] is NA
              if (!is.na(exper@pop.rescale.each.timestep[t])) {state <- exper@pop.rescale.each.timestep[t] *(state/sum(state))}
              
              #put in the correct birth rate for that time-step, if it varies
              if (length(exper@births.per.1000.each.timestep)>1) {
                # Read 'birth rate'
                tmp.trans@birth.rate = (exper@births.per.1000.each.timestep[t]*sum(state)/1000)
              } else {
                tmp.trans@birth.rate <- exper@trans@birth.rate*sum(rc[,t-1])/sum(exper@state.t0)
              }
              births.each.timestep[t] <-  tmp.trans@birth.rate
              
              #put in time appropriate survival rate otherwise it uses original surv.matrix set up for trans object and keeps constant over time
              if (!is.na(exper@surv.each.timestep[1,1])) {
                tmp.trans@age.surv.matrix <- ExtractAgeSpecificSurvivalMatrix(tmp.trans=tmp.trans, maternal.obj=exper@maternal.obj, 
                                                                              surv.at.timestep.t=exper@surv.each.timestep[,t])
              }
              
              ##put in correct vaccination coverage
              tmp.trans@vac.per@pvacc.in.age.class <- rep(0, exper@trans@n.age.class) #start with clean slate each time step
              routine.vacc.prob <- routine$age.time.specific.routine
              sia.vacc.prob <- SIA$age.time.specific.SIA
              
              if (!is.na(index.sia.vacc[t])){ #if SIA
                
                #IF SIA independent MR1 or MR2
                if (!exper@MR1SIAcorrelation & !exper@MR2SIAcorrelation & !exper@SIAinacc & !exper@SIAinefficient){ 
                  tmp.trans@vac.per@pvacc.in.age.class <- 
                    routine.vacc.prob[index.routine.vacc[t],] + 
                    sia.vacc.prob[index.sia.vacc[t],] - #xxamy changed from addition to subtraction
                    (routine.vacc.prob[index.routine.vacc[t],]*sia.vacc.prob[index.sia.vacc[t],]) 
                }
                
                #IF SIA and MR1 correlated, SIA & MR2 NOT correlated
                if (exper@MR1SIAcorrelation & !exper@MR2SIAcorrelation & !exper@SIAinacc & !exper@SIAinefficient) {
                  MR1.age.range <- exper@time.specific.min.age.MR1[index.routine.vacc[t]]:exper@time.specific.max.age.MR1[index.routine.vacc[t]]
                  MR2.age.range <- exper@time.specific.min.age.MR2[index.routine.vacc[t]]:exper@time.specific.max.age.MR2[index.routine.vacc[t]]
                  if (exper@time.specific.SIAcov[index.sia.vacc[t]]<exper@time.specific.MR1cov[index.routine.vacc[t]]){
                    #non MR1 or MR2 ages
                    tmp.trans@vac.per@pvacc.in.age.class <- sia.vacc.prob[index.sia.vacc[t],]
                    #MR1 ages
                    tmp.trans@vac.per@pvacc.in.age.class[MR1.age.range] <- 
                      routine.vacc.prob[index.routine.vacc[t],MR1.age.range] + routine$one.minus.ve1[index.routine.vacc[t]]*sia.vacc.prob[index.sia.vacc[t],MR1.age.range]
                    #MR2 ages
                    tmp.trans@vac.per@pvacc.in.age.class[MR2.age.range] <- 
                      routine.vacc.prob[index.routine.vacc[t],MR2.age.range] + 
                      sia.vacc.prob[index.sia.vacc[t],MR2.age.range] + 
                      (routine.vacc.prob[index.routine.vacc[t],MR2.age.range]*sia.vacc.prob[index.sia.vacc[t],MR2.age.range]) 
                  }
                  if (exper@time.specific.SIAcov[index.sia.vacc[t]]>exper@time.specific.MR1cov[index.routine.vacc[t]]){
                    #non MR1 or MR2 ages
                    tmp.trans@vac.per@pvacc.in.age.class <- sia.vacc.prob[index.sia.vacc[t],]
                    #MR1 ages
                    tmp.trans@vac.per@pvacc.in.age.class[MR1.age.range] <- 
                      routine.vacc.prob[index.routine.vacc[t],MR1.age.range] + routine$prop.fail.MR1.byage[index.routine.vacc[t],MR1.age.range] + 
                      (sia.vacc.prob[index.sia.vacc[t],MR1.age.range] - routine.vacc.prob[index.routine.vacc[t],MR1.age.range])
                    #MR2 ages
                    tmp.trans@vac.per@pvacc.in.age.class[MR2.age.range] <- 
                      routine.vacc.prob[index.routine.vacc[t],MR2.age.range] + 
                      sia.vacc.prob[index.sia.vacc[t],MR2.age.range] + 
                      (routine.vacc.prob[index.routine.vacc[t],MR2.age.range]*sia.vacc.prob[index.sia.vacc[t],MR2.age.range]) 
                  }
                }
                
                #IF SIA and MR1 & SIA and MR2 correlated
                if (exper@MR1SIAcorrelation & exper@MR2SIAcorrelation & !exper@SIAinacc & !exper@SIAinefficient) {
                  MR1.age.range <- exper@time.specific.min.age.MR1[index.routine.vacc[t]]:exper@time.specific.max.age.MR1[index.routine.vacc[t]]
                  MR2.age.range <- exper@time.specific.min.age.MR2[index.routine.vacc[t]]:exper@time.specific.max.age.MR2[index.routine.vacc[t]]
                  if (exper@time.specific.SIAcov[index.sia.vacc[t]]<exper@time.specific.MR1cov[index.routine.vacc[t]]){
                    #non MR1 or MR2 ages
                    tmp.trans@vac.per@pvacc.in.age.class <- sia.vacc.prob[index.sia.vacc[t],]
                    #MR1 ages
                    tmp.trans@vac.per@pvacc.in.age.class[MR1.age.range] <- 
                      routine.vacc.prob[index.routine.vacc[t],MR1.age.range] + routine$one.minus.ve1[index.routine.vacc[t]]*sia.vacc.prob[index.sia.vacc[t],MR1.age.range]
                    #MR2 ages
                    tmp.trans@vac.per@pvacc.in.age.class[MR2.age.range] <- 
                      routine.vacc.prob[index.routine.vacc[t],MR2.age.range] + routine$one.minus.ve1[index.routine.vacc[t]]*sia.vacc.prob[index.sia.vacc[t],MR2.age.range]
                  }
                  if (exper@time.specific.SIAcov[index.sia.vacc[t]]>exper@time.specific.MR1cov[index.routine.vacc[t]]){
                    #non MR1 or MR2 ages
                    tmp.trans@vac.per@pvacc.in.age.class <- sia.vacc.prob[index.sia.vacc[t],]
                    #MR1 ages
                    tmp.trans@vac.per@pvacc.in.age.class[MR1.age.range] <- 
                      routine.vacc.prob[index.routine.vacc[t],MR1.age.range] + routine$prop.fail.MR1.byage[index.routine.vacc[t],MR1.age.range] + 
                      (sia.vacc.prob[index.sia.vacc[t],MR1.age.range] - routine.vacc.prob[index.routine.vacc[t],MR1.age.range])
                    #MR2 ages
                    tmp.trans@vac.per@pvacc.in.age.class[MR2.age.range] <- 
                      routine.vacc.prob[index.routine.vacc[t],MR2.age.range] + routine$prop.fail.MR1.byage[index.routine.vacc[t],MR2.age.range] + 
                      (sia.vacc.prob[index.sia.vacc[t],MR2.age.range] - routine.vacc.prob[index.routine.vacc[t],MR2.age.range])
                  }
                }
                
                #If inaccessible population, but not SIA inefficiency
                if (!exper@MR1SIAcorrelation & !exper@MR2SIAcorrelation & exper@SIAinacc & !exper@SIAinefficient) {
                  tmp.trans@vac.per@pvacc.in.age.class <- 1-(exper@prop.inacc[t]+
                                                               (1-exper@prop.inacc[t])*(1-(routine.vacc.prob[index.routine.vacc[t],] + 
                                                                                             sia.vacc.prob[index.sia.vacc[t],] + 
                                                                                             (routine.vacc.prob[index.routine.vacc[t],]*sia.vacc.prob[index.sia.vacc[t],]))))
                  #tmp.trans@vac.per@pvacc.in.age.class <- 1-(exper@prop.inacc[t]+(1-exper@prop.inacc[t])*(1-routine.vacc.prob[index.routine.vacc[t],])*(1-sia.vacc.prob[index.sia.vacc[t],]))
                }
                
                #If inaccessible population AND SIA inefficiency
                if (!exper@MR1SIAcorrelation & !exper@MR2SIAcorrelation & exper@SIAinacc & exper@SIAinefficient) {
                  # with SIA inaccessible & with SIA inefficiency -- VIMC 2017-2021 versions 
                  z <- routine.vacc.prob[index.routine.vacc[t],] #prob. successful routine given access
                  ro <- (1-exper@prop.inacc[t]) #prob. accessible
                  m <- (1-(sia.vacc.prob[index.sia.vacc[t],]*(1/ro)*(1-0.1)))^(1/(1-0.1)) #prob. not successful campaign vaccination given access, and inefficiency of 0.1
                  m[is.na(m)] <- 0 #if NaN then negative number was raised, which means everyone in accessible population was vaccinated by campaign, therefore prop. not vaccinated =0
                  age.spec.vacc.prob <- 1-((1-ro)+(ro*m)) #ages over 36 months, where only eligible for campaign
                  max.age.routine <- exper@time.specific.max.age.MR2[index.routine.vacc[t]]
                  age.spec.vacc.prob[1:max.age.routine] <- 1-((1-ro)+(ro*m*(1-(z/ro))))[1:max.age.routine] #ages 1:36 months where routine takes place
                  #if(sum(age.spec.vacc.prob)>0) print(age.spec.vacc.prob[1:max.age.routine])
                  tmp.trans@vac.per@pvacc.in.age.class <- age.spec.vacc.prob
                }
                
                #stow output
                MR1.fail.each.timestep[t] <- routine$prop.fail.MR1[index.routine.vacc[t]]
                MR2.fail.each.timestep[t] <- routine$prop.fail.MR2[index.routine.vacc[t]]
                SIA.fail.each.timestep[t] <- SIA$prop.fail.SIA[index.sia.vacc[t]]
                
              } else { #if no SIA
                tmp.trans@vac.per@pvacc.in.age.class <- routine.vacc.prob[index.routine.vacc[t],]
                #stow output
                MR1.fail.each.timestep[t] <- routine$prop.fail.MR1[index.routine.vacc[t]]
                MR2.fail.each.timestep[t] <- routine$prop.fail.MR2[index.routine.vacc[t]]
              }
              
              N0 <- sum(state)
              state <- next.ID.state(state, tmp.trans)
              NT <- sum(state)
              growth.rate.each.timestep[t] <- log(NT/N0) #instantaneous biweekly growth rate
              
              #print(t)
              #print(dim(state))
              rc [,t] <- state
              
            }
            
            
            rc <- new("sim.results.MSIRV.update.demog.vaccine.change",
                      data=rc,
                      m.inds = exper@trans@m.inds,
                      s.inds = exper@trans@s.inds,
                      i.inds = exper@trans@i.inds,
                      r.inds = exper@trans@r.inds,
                      v.inds = exper@trans@v.inds,
                      t = exper@t.min+(1:T-1)* exper@step.size,
                      age.class = exper@trans@age.class,
                      births.each.timestep = births.each.timestep,
                      growth.rate.each.timestep = growth.rate.each.timestep,
                      MR1.fail.each.timestep = MR1.fail.each.timestep,
                      MR2.fail.each.timestep = MR2.fail.each.timestep, 
                      SIA.fail.each.timestep = SIA.fail.each.timestep,
                      routine.intro = routine.intro,
                      sia.times = sia.times)
            
            rc <- new("experiment.result",
                      experiment.def = exper,
                      result = rc)
            
            
            return(rc)
          }
)

#Function to set up experiment and run the transients out
#
#Parameters -
#     uncode - UN country code
#     generation.time - the desired generation time in months
#     age.classes - vector - the upper limit of the age classes we are interested in months
#     maternal.decay.rt -rate of maternal decay of immunity to rubella
#     exponent - numeric - exponent for the infected
#     tot.pop - numeric - population you want to scale you whole experiment by at the beginning, NULL or numeric
#     yr.births.per.1000.acrossyears - vector of lenth(t.max) - crude birth rate per 1000 per year
#     coverage - proportion - vaccination coverage
#     targeted.intro - boolean - FALSE means space introduction out over all age classes, TRUE means to concentrate introductions
#     intro.rate - numeric - rate at which new infections are introduced
#     R0 - numeric - basic reproductive number of rubella assumed, used to rescale the waifw
#     t.max - numeric - time in year that plan to run out transients
#     get.births - vector of births with length of 5*n.age.classes
#     seasonal.amp - numeric - amplitude of the seasonal cosine function
#     flat.WAIFW - boolean - if false then default waifw is polymod GB
#     country.specific.WAIFW - boolean - if true then overrides polymod or flat waifw
#     vynnycky.waifw - boolean - true if want emilia's waifw <13 and 15-50 (beta_young, beta_old, and 0.7*beta_old)
#     vynnycky.waifw.betas - vector of two betas (beta_young and beta_old)
#     asdr.object - nMx object with rate and mid-age (age specific death rates)
#     year - numeric - year to pull DFE demography (population and age structure)
#     use_montagu_demog - boolean - T = use montagu demography or F = use UNPD demography
#
#Returns -
#     tran object and state.t0 after transients run out
EX.Country.part1 <- function(uncode,
                             generation.time = 0.5, #generation time in months
                             age.classes = c(1:240, seq(241,720,12)),
                             maternal.decay.rt=0.95, #based on Metcalf, 2012 paper
                             exponent = 0.97,
                             frequency.dep=T,
                             tot.pop=NULL,
                             yr.births.per.1000.acrossyears,
                             intro.rate=1e-10,
                             targeted.intro=F,
                             R0=5,
                             t.max = 20,
                             get.births=NULL,
                             seasonal.amp=0.35,#metcalf et al 2012 used 0.35
                             flat.WAIFW=F,
                             country.specific.WAIFW=F,
                             vynnycky.waifw=F,
                             vynnycky.waifw.betas = c(2,1),
                             asdr.object,
                             year=1990,
                             use_montagu_demog=F,
                             routine.vac=0,
                             routine.vac.age.index=12) { 
  
  tmp <- Get.CountryX.Starting.Pop.MSIRV(uncode=uncode, 
                                         generation.time=generation.time,
                                         age.classes=age.classes,
                                         maternal.decay.rt=maternal.decay.rt,
                                         exponent=exponent,
                                         frequency.dep = frequency.dep,
                                         is.stochastic = FALSE, # Make NOT stochastic for the transient ridding stage
                                         tot.pop=tot.pop, 
                                         yr.births.per.1000 = yr.births.per.1000.acrossyears[1], 
                                         intro.rate=intro.rate,
                                         flat.WAIFW=flat.WAIFW,
                                         asdr.object=asdr.object,
                                         targeted.intro=targeted.intro,
                                         year=year,
                                         get.births=get.births,
                                         use_montagu_demog=use_montagu_demog,
                                         routine.vac=routine.vac,
                                         routine.vac.age.index=routine.vac.age.index)
  
  EX <- new("experiment.updatedemog")
  EX@name <- paste(uncode, "Frequency Dependent", sep=" ") 
  EX@state.t0 = tmp$state
  EX@trans = tmp$tran
  EX@t.min = 0
  EX@t0.doy = 0
  EX@step.size = 1/(12*1/generation.time) #1/no.gens.per.year
  EX@season.obj = new("seasonal.cosine", amplitude = seasonal.amp)
  EX@R0 <- R0 
  EX@t.max <-t.max 
  EX@maternal.obj <- new("maternal.exp.decay", decay.rt=maternal.decay.rt)
  
  #change the WAIFW to country-specific if requested (based on prem et al.)
  if (country.specific.WAIFW & use_montagu_demog) EX@trans@waifw <- get.prem.WAIFW(age.class.boundries=age.classes/12, 
                                                                                   uncode, other.contact.matrix=F,
                                                                                   bandwidth=c(3,3), 
                                                                                   adjustment_start_year=T, year=1980)
  
  if (country.specific.WAIFW & !use_montagu_demog) EX@trans@waifw <- get.prem.WAIFW(age.class.boundries=age.classes/12, 
                                                                                    uncode, other.contact.matrix=F,
                                                                                    bandwidth=c(3,3), 
                                                                                    adjustment_start_year=F, year=1980)
  
  if (vynnycky.waifw) EX@trans@waifw <- get.vynnycky.WAIFW(age.classes/12, 
                                                           beta_young=vynnycky.waifw.betas[1],
                                                           beta_old=vynnycky.waifw.betas[2])
  
  # To run out the transients, keep the birth rate constant across time
  EX@births.per.1000.each.timestep <- c(yr.births.per.1000.acrossyears[1], rep(yr.births.per.1000.acrossyears, each=12*1/generation.time))*generation.time/12
  
  # Setting up survival over time as NULL because not changing over time
  EX@surv.each.timestep <-  matrix(NA,1,1)
  
  # Rescale the pop size (take the popuation distibution from state.t0 and multiply by POP) if you want to
  if (!is.null(tot.pop)) EX@state.t0[,1] <- tot.pop*EX@state.t0[,1]/sum(EX@state.t0[,1])
  
  # seed an infected and run out to equilibrium
  EX@state.t0[3+5*11,1] <-sum(EX@state.t0[,1])*0.00001 # seeding an infected person into i.inds because ind is 3rd epi class out of 5 therefore every 8 is an infected class - and want it intereted at age 11 or 88th row
  EX@state.t0[2+5*11,1] <-EX@state.t0[2+5*11,1]-sum(EX@state.t0[,1])*0.00001 # to keep population the same size we are pulling this 100 from susceptibles at age 11
  
  # This condition was initiating a bunch of vaccinated - change to recovered
  EX@state.t0[EX@trans@r.inds,1][EX@trans@age.class>60] <-EX@state.t0[EX@trans@s.inds,1][EX@trans@age.class>60]
  EX@state.t0[EX@trans@s.inds,1][EX@trans@age.class>60] <-0
  
  # Get rid of negatives
  EX@state.t0[which(EX@state.t0[,1]<0),1] <- 0
  
  # No rescaling of population - 
  EX@pop.rescale.each.timestep <- 0 #--xxamy - added this lines because when changed pop.rescale.each.timestep to "ANY" is became NULL by default but we need it to be a 0
  
  # Run out transients - plot if you want to be sure...
  res <- run(EX, rescale.WAIFW=T) 
  #plot(res@result)
  
  # Reset
  age.struct.1991 <- GetNumber.per.AgeGroup(state=res@experiment.def@state.t0, trans=EX@trans)
  age.struc.sim <- GetNumber.per.AgeGroup(state=res@result@.Data[,ncol(res@result@.Data)], trans=EX@trans)
  prop.struc.sim <- res@result[,ncol(res@result)]/rep(age.struc.sim, each=5)
  new.state <- rep(age.struct.1991, each=5)*prop.struc.sim
  EX@state.t0[,1] <- new.state
  #EX@state.t0[,1] <- sum(res@result[,1])*res@result[,ncol(res@result)]/sum(res@result[,ncol(res@result)]) 
  EX@trans@waifw <- res@experiment.def@trans@waifw
  
  res <- run(EX, rescale.WAIFW=T)
  #plot(res@result)
  
  # Reset
  age.struct.1991 <- GetNumber.per.AgeGroup(state=res@experiment.def@state.t0, trans=EX@trans)
  age.struc.sim <- GetNumber.per.AgeGroup(state=res@result@.Data[,ncol(res@result@.Data)], trans=EX@trans)
  prop.struc.sim <- res@result[,ncol(res@result)]/rep(age.struc.sim, each=5)
  new.state <- rep(age.struct.1991, each=5)*prop.struc.sim
  EX@state.t0[,1] <- new.state
  #EX@state.t0[,1] <- sum(res@result[,1])*res@result[,ncol(res@result)]/sum(res@result[,ncol(res@result)]) 
  EX@trans@waifw <- res@experiment.def@trans@waifw
  
  res <- run(EX, rescale.WAIFW=T)
  #plot(res@result)
  
  # Reset
  age.struct.1991 <- GetNumber.per.AgeGroup(state=res@experiment.def@state.t0, trans=EX@trans)
  age.struc.sim <- GetNumber.per.AgeGroup(state=res@result@.Data[,ncol(res@result@.Data)], trans=EX@trans)
  prop.struc.sim <- res@result[,ncol(res@result)]/rep(age.struc.sim, each=5)
  new.state <- rep(age.struct.1991, each=5)*prop.struc.sim
  EX@state.t0[,1] <- new.state
  #EX@state.t0[,1] <- sum(res@result[,1])*res@result[,ncol(res@result)]/sum(res@result[,ncol(res@result)]) 
  EX@trans@waifw <- res@experiment.def@trans@waifw
  
  return(EX)
}

