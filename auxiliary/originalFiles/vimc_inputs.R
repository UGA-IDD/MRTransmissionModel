
# grant-specific inputs
# these functions priority


## Functions to create inputs based on VIMC (i.e., Montagu) vaccination and demography data
## Authors: Amy Winter

### Function to get population size, age structure, fertility and death rates from 1950 to 2105
### the last five years I just keep all 2100 rates constant
### based on data pulled from VIMC's Montagu
###
##' @param - uncode - UN country code
###
### Returns demography for each country
getDemography <- function(uncode=NA){

  cc <- uncode
  time <- seq(1950,2105,1) #revised to 2105 xxamy
  age.classes <- 0:100

  #Population
  df <- read.csv("./data/montagu_dem_data/202110gavi_v3/pop.csv")
  pop.total.1950.2100 <- df[which(df$country_code_numeric==cc),"value"]
  pop.total.1950.2100 <- c(pop.total.1950.2100, rep(pop.total.1950.2100[length(pop.total.1950.2100)], 5))/1000

  #Age Structure
  df <- read.csv("./data/montagu_dem_data/202110gavi_v3/pop.age.csv")
  pop.age.byageclasses.1950.2100 <- matrix(df[which(df$country_code_numeric==cc & df$year>1949),"value"], nrow=101, ncol=151)
  add <- matrix(rep(pop.age.byageclasses.1950.2100[,151], each=5), ncol=5, byrow=T)
  pop.age.byageclasses.1950.2100 <- cbind(pop.age.byageclasses.1950.2100, add)/1000
  rownames(pop.age.byageclasses.1950.2100) <- age.classes
  colnames(pop.age.byageclasses.1950.2100) <- time

  #Sex Distribution by age and over time
  dff <- read.csv("./data/montagu_dem_data/202110gavi_v3/pop.age.fem.csv")
  pop.age.byageclasses.1950.2100.f <- matrix(dff[which(dff$country_code_numeric==cc & dff$year>1949),"value"], nrow=101, ncol=151)
  add <- matrix(rep(pop.age.byageclasses.1950.2100.f[,151], each=5), ncol=5, byrow=T)
  pop.age.byageclasses.1950.2100.f <- cbind(pop.age.byageclasses.1950.2100.f, add)/1000
  sex.dist <- pop.age.byageclasses.1950.2100.f/pop.age.byageclasses.1950.2100
  #Sex Distribution grouped by 5 year reproductive ages over time
  repro.ages.index <- findInterval(age.classes, seq(15,50,5))
  repro.age.sex.dist.1950.2100 <- matrix(NA, nrow=length(unique(repro.ages.index)), ncol=ncol(sex.dist))
  for (i in 1:length(unique(repro.ages.index))){
    index <- unique(repro.ages.index)[i]
    repro.age.sex.dist.1950.2100[i,] <- apply(sex.dist, 2, function(x) mean(x[repro.ages.index==index]))
  }
  repro.age.sex.dist.1950.2100 <- repro.age.sex.dist.1950.2100[-c(1,length(unique(repro.ages.index))),]
  rownames(repro.age.sex.dist.1950.2100) <- c("15-19","20-24","25-29","30-34","35-39","40-44","45-49")
  colnames(repro.age.sex.dist.1950.2100) <- time

  #Life expectancy
  df <- read.csv("./data/montagu_dem_data/202110gavi_v3/e0.csv")
  e0.1950.2100 <- df[which(df$country_code_numeric==cc),"value"]
  e0.1950.2100 <- c(e0.1950.2100, rep(e0.1950.2100[length(e0.1950.2100)], 5))

  #age-specific fertility rate
  df <- read.csv("./data/montagu_dem_data/202110gavi_v3/asfr.csv")
  asfr.1950.2100.by5 <- matrix(df[which(df$country_code_numeric==cc),"value"], nrow=length(unique(df$age_from)), ncol=length(unique(df$year)))
  #asfr.1950.2100 <- matrix(apply(asfr.1950.2100.by5, 1, function(x) rep(x, each=5)), ncol=(5*ncol(asfr.1950.2100.by5)), byrow=T)
  #smooth over years instead of just repreating as above line does - xxamy
  asfr.1950.2100 <- matrix(NA, 7, length(1950:2099))
  for (p in 1:nrow(asfr.1950.2100)){
    asfr.1950.2100[p,] <- predict(smooth.spline(seq(1952.5,2097.5,5), asfr.1950.2100.by5[p,]), 1950:2099)$y
  }
  add <- matrix(rep(asfr.1950.2100[,150], each=6), ncol=6, byrow=T)
  asfr.1950.2100 <- cbind(asfr.1950.2100, add)
  colnames(asfr.1950.2100) <- time
  asfr.1950.2100[which(asfr.1950.2100<0)] <- 0 #if less than 0, force to 0 (issue with Armenia)

  #crude birth rate
  df <- read.csv("./data/montagu_dem_data/202110gavi_v3/cbr.csv")
  cbr.1950.2100 <- df[which(df$country_code_numeric==cc),"value"]
  cbr.1950.2100 <- c(cbr.1950.2100, rep(cbr.1950.2100[length(cbr.1950.2100)], 6))*1000
  names(cbr.1950.2100) <- time

  #Age Specific Death Rates
  df <- read.csv("./data/montagu_dem_data/202110gavi_v3/asdr.csv")
  asdr.1950.2100.by5 <- matrix(df[which(df$country_code_numeric==cc),"value"],  nrow=length(unique(df$age_from)), ncol=length(unique(df$year)))
  asdr.1950.2100.by5 <- cbind(asdr.1950.2100.by5, asdr.1950.2100.by5[,length(unique(df$year))])
  mid.age <- unique(df$age_from) + (unique(df$age_to)-unique(df$age_from))/2
  mid.age[1] <- 0.5
  asdr.object <- new("nMx",
                     rates = as.data.frame(asdr.1950.2100.by5),
                     mid.age = mid.age)

  #Age Specific Life Expectancy
  #df <- read.csv("./data/montagu_dem_data/201710gavi_v5/le.age.csv")
  #asle.1950.2100.by5 <- matrix(df[which(df$country_code_numeric==cc),"value"], nrow=22, ncol=30)
  #asle.1950.2100 <- matrix(apply(asle.1950.2100.by5, 1, function(x) rep(x, each=5)), ncol=(5*ncol(asle.1950.2100.by5)), byrow=T)
  #add <- matrix(rep(asle.1950.2100[,150], each=6), ncol=6, byrow=T)
  #asle.1950.2100 <- cbind(asle.1950.2100, add)
  #asle.1950.2100a <- matrix(rep(asle.1950.2100[2,], each=4), nrow=4, byrow=F) #constant over four year ages
  #asle.1950.2100b <- matrix(apply(asle.1950.2100[3:21,], 2, function(x) rep(x, each=5)), nrow=(5*(nrow(asle.1950.2100.by5)-3)), byrow=F) #constant over five year ages
  #asle.1950.2100 <- rbind(asle.1950.2100[1,], asle.1950.2100a, asle.1950.2100b, asle.1950.2100[22,])
  #colnames(asle.1950.2100) <- time
  #rownames(asle.1950.2100) <- age.classes

  #Remainig from original function
  tfr.1950.2100=NA
  births.1950.2100=NA


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


### Function to pull in vaccination data from Montagu 202110sgavi_v3 touchstone vaccine scenarios
###
##' @param - iso3code - iso3code
###
### Returns age, year, and coverage for SIA, RCV1, and RCV2 for each vaccination scenario
getMontaguCoverage_202110gavi_v3 <- function(iso3code){

  #Also Estimating inaccessible at the same time, because it depends on routine estimates
  #Inaccessible from DPT1 numbers over time
  dtp <- read.csv("./data/dtp1_estimates.csv")
  dtp.c <- dtp[dtp$ISO_code==iso3code,]
  dtp.c <- dtp.c[order(dtp.c$Year),]
  dtp.1980.to.2100 <- dtp.c$vacc_coverage[1:length(1980:2017)]
  #filling in zeros using zoo package
  dtp.1980.to.2100 <- zoo::na.locf(dtp.1980.to.2100, fromLast=T, na.rm=F) #fills in using value after NA
  dtp.1980.to.2100 <- zoo::na.locf(dtp.1980.to.2100) #fills in using previous value before NA
  dtp.1980.to.2100 <- c(dtp.1980.to.2100, rep(dtp.1980.to.2100[length(dtp.1980.to.2100)], length(2018:2101))) #add one more year 1980to2101
  if (length(dtp.1980.to.2100)==0){
    dtp.1980.to.2100 <- rep(1, 122)
  }

  #CAMPAIGN DEFAULT SCENARIO
  data <- read.csv("./data/montagu_coverage_data/202110gavi_v3/rubella-campaign-default.csv")
  #unique(data$activity_type)
  #unique(data$vaccine) #check vaccine types
  if (length(unique(data$activity_type))>1) print("campaign default has other activity_types") #check activity types

  dfx <- data[data$activity_type=="campaign" & data$country_code==iso3code & data$coverage!=0,]
  sia.coverage.1980.2100 <- sia.min.age.1980.2100 <- sia.max.age.1980.2100 <- rep(0, (2101-1980))
  year.indeces <- as.numeric(dfx$year)-1980+1
  sia.coverage.1980.2100[year.indeces] <- as.numeric(dfx$coverage)
  sia.coverage.1980.2100 <- ifelse(sia.coverage.1980.2100>1, 1, sia.coverage.1980.2100) #cant have coverage over 1.
  sia.min.age.1980.2100[year.indeces] <- ifelse(dfx$age_first==0, 1, as.numeric(dfx$age_first*12))
  sia.max.age.1980.2100[year.indeces] <- as.numeric(dfx$age_last*12)
  if (nrow(dfx)==0){  #if no SIA, it needs to be in this form
    sia.coverage.1980.2100 <- sia.min.age.1980.2100 <- sia.max.age.1980.2100 <- rep(0, (2101-1980))
  }
  rcv1.coverage.1980.2100 <- rep(0, (2101-1980))
  rcv2.coverage.1980.2100 <- rep(0, (2101-1980))
  campaign_default <- new("cov.estimates",
                          sia.cov = sia.coverage.1980.2100,
                          sia.min.age = sia.min.age.1980.2100,
                          sia.max.age = sia.max.age.1980.2100,
                          MR1.cov = rcv1.coverage.1980.2100,
                          MR2.cov = rcv2.coverage.1980.2100,
                          inaccessible.prop = (1-dtp.1980.to.2100))


  #CAMPAIGN BESTCASE SCENARIO
  data <- read.csv("./data/montagu_coverage_data/202110gavi_v3/rubella-campaign-ia2030_target.csv")
  #unique(data$activity_type) #check activity types
  #unique(data$vaccine) #check vaccine types
  if (length(unique(data$activity_type))>1) print("campaign bestcase has other activity_types") #check activity types

  dfx <- data[data$activity_type=="campaign" & data$country_code==iso3code & data$coverage!=0,]
  sia.coverage.1980.2100 <- sia.min.age.1980.2100 <- sia.max.age.1980.2100 <- rep(0, (2101-1980))
  year.indeces <- as.numeric(dfx$year)-1980+1
  sia.coverage.1980.2100[year.indeces] <- as.numeric(dfx$coverage)
  sia.coverage.1980.2100 <- ifelse(sia.coverage.1980.2100>1, 1, sia.coverage.1980.2100) #cant have coverage over 1.
  sia.min.age.1980.2100[year.indeces] <- ifelse(dfx$age_first==0, 1, as.numeric(dfx$age_first*12))
  sia.max.age.1980.2100[year.indeces] <- as.numeric(dfx$age_last*12)
  if (nrow(dfx)==0){  #if no SIA, it needs to be in this form
    sia.coverage.1980.2100 <- sia.min.age.1980.2100 <- sia.max.age.1980.2100 <- rep(0, (2101-1980))
  }
  rcv1.coverage.1980.2100 <- rep(0, (2101-1980))
  rcv2.coverage.1980.2100 <- rep(0, (2101-1980))
  campaign_bestcase <- new("cov.estimates",
                           sia.cov = sia.coverage.1980.2100,
                           sia.min.age = sia.min.age.1980.2100,
                           sia.max.age = sia.max.age.1980.2100,
                           MR1.cov = rcv1.coverage.1980.2100,
                           MR2.cov = rcv2.coverage.1980.2100,
                           inaccessible.prop = (1-dtp.1980.to.2100))



  #RCV1 DEFAULT SCENARIO
  data <- read.csv("./data/montagu_coverage_data/202110gavi_v3/rubella-rcv1-default.csv")
  #unique(data$activity_type) #check activity types
  #unique(data$vaccine) #check vaccine types
  if (any(unique(data$vaccine)=="RCV2")) print("RCV2 listed in rcv1_default")
  dfx.tmp <- data[data$vaccine=="Rubella" & data$activity_type=="routine" & data$coverage!=0,]
  #table(dfx.tmp$age_last) #check ages for RCV1
  if (any(unique(dfx.tmp$age_last)!=0)) print("RCV1 is not at 0 years old in rcv1_default")

  dfx <- data[data$activity_type=="campaign" & data$country_code==iso3code & data$coverage!=0,]
  sia.coverage.1980.2100 <- sia.min.age.1980.2100 <- sia.max.age.1980.2100 <- rep(0, (2101-1980))
  year.indeces <- as.numeric(dfx$year)-1980+1
  sia.coverage.1980.2100[year.indeces] <- as.numeric(dfx$coverage)
  sia.coverage.1980.2100 <- ifelse(sia.coverage.1980.2100>1, 1, sia.coverage.1980.2100) #cant have coverage over 1.
  sia.min.age.1980.2100[year.indeces] <- ifelse(dfx$age_first==0, 1, as.numeric(dfx$age_first*12))
  sia.max.age.1980.2100[year.indeces] <- as.numeric(dfx$age_last*12)
  if (nrow(dfx)==0){  #if no SIA, it needs to be in this form
    sia.coverage.1980.2100 <- sia.min.age.1980.2100 <- sia.max.age.1980.2100 <- rep(0, (2101-1980))
  }
  dfx <- data[data$activity_type=="routine" & data$country_code==iso3code & data$coverage!=0,]
  if (nrow(dfx)!=0){
    year.range <- range(dfx$year)
    if (any(diff(dfx$year)!=1)){ #issue in IRQ
      rcv1.coverage <- c(dfx$coverage[1:which(diff(dfx$year)!=1)],rep(0,(diff(dfx$year)[which(diff(dfx$year)!=1)]-1)), dfx$coverage[(which(diff(dfx$year)!=1)+1):length(dfx$coverage)])
      rcv1.coverage.1980.2100 <- c(rep(0, (year.range[1]-1980)), rcv1.coverage, rep(0, (2100-year.range[2])))
    } else {
      rcv1.coverage <- as.numeric(dfx$coverage)
      rcv1.coverage.1980.2100 <- c(rep(0, (year.range[1]-1980)), rcv1.coverage, rep(0, (2100-year.range[2])))
    }
  } else {
    rcv1.coverage.1980.2100 <- rep(0, (2101-1980))
  }
  rcv2.coverage.1980.2100 <- rep(0, (2101-1980))
  #if dtp1 < RCV1 - force DTP1 to be RCV1
  max.rcv <- rcv1.coverage.1980.2100
  rcv.extra1 <- c(max.rcv, max.rcv[length(max.rcv)])
  dtp.1980.to.2100.tmp <- dtp.1980.to.2100
  dtp.1980.to.2100.tmp[dtp.1980.to.2100<rcv.extra1] <- rcv.extra1[dtp.1980.to.2100<rcv.extra1]
  rcv1_default <- new("cov.estimates",
                      sia.cov = sia.coverage.1980.2100,
                      sia.min.age = sia.min.age.1980.2100,
                      sia.max.age = sia.max.age.1980.2100,
                      MR1.cov = rcv1.coverage.1980.2100,
                      MR2.cov = rcv2.coverage.1980.2100,
                      inaccessible.prop = (1-dtp.1980.to.2100))


  #RCV1 BESTCASE SCENARIO
  data <- read.csv("./data/montagu_coverage_data/202110gavi_v3/rubella-rcv1-ia2030_target.csv")
  #unique(data$activity_type) #check activity types
  #unique(data$vaccine) #check vaccine types
  if (any(unique(data$vaccine)=="RCV2")) print("RCV2 listed in rcv1_bestcase")
  dfx.tmp <- data[data$vaccine=="Rubella" & data$activity_type=="routine" & data$coverage!=0,]
  #table(dfx.tmp$age_last) #check ages for RCV1
  if (any(unique(dfx.tmp$age_last)!=0)) print("RCV1 is not at 0 years old in rcv1_bestcase")

  dfx <- data[data$activity_type=="campaign" & data$country_code==iso3code & data$coverage!=0,]
  sia.coverage.1980.2100 <- sia.min.age.1980.2100 <- sia.max.age.1980.2100 <- rep(0, (2101-1980))
  year.indeces <- as.numeric(dfx$year)-1980+1
  sia.coverage.1980.2100[year.indeces] <- as.numeric(dfx$coverage)
  sia.coverage.1980.2100 <- ifelse(sia.coverage.1980.2100>1, 1, sia.coverage.1980.2100) #cant have coverage over 1.
  sia.min.age.1980.2100[year.indeces] <- ifelse(dfx$age_first==0, 1, as.numeric(dfx$age_first*12))
  sia.max.age.1980.2100[year.indeces] <- as.numeric(dfx$age_last*12)
  if (nrow(dfx)==0){  #if no SIA, it needs to be in this form
    sia.coverage.1980.2100 <- sia.min.age.1980.2100 <- sia.max.age.1980.2100 <- rep(0, (2101-1980))
  }
  dfx <- data[data$activity_type=="routine" & data$country_code==iso3code & data$coverage!=0,]
  if (nrow(dfx)!=0){
    year.range <- range(dfx$year)
    if (any(diff(dfx$year)!=1)){
      rcv1.coverage <- c(dfx$coverage[1:which(diff(dfx$year)!=1)],rep(0,(diff(dfx$year)[which(diff(dfx$year)!=1)]-1)), dfx$coverage[(which(diff(dfx$year)!=1)+1):length(dfx$coverage)])
      rcv1.coverage.1980.2100 <- c(rep(0, (year.range[1]-1980)), rcv1.coverage, rep(0, (2100-year.range[2])))
    } else {
      rcv1.coverage <- as.numeric(dfx$coverage)
      rcv1.coverage.1980.2100 <- c(rep(0, (year.range[1]-1980)), rcv1.coverage, rep(0, (2100-year.range[2])))
    }
  } else {
    rcv1.coverage.1980.2100 <- rep(0, (2101-1980))
  }
  rcv2.coverage.1980.2100 <- rep(0, (2101-1980))
  #if dtp1 < RCV1 - force DTP1 to be RCV1
  max.rcv <- rcv1.coverage.1980.2100
  rcv.extra1 <- c(max.rcv, max.rcv[length(max.rcv)])
  dtp.1980.to.2100.tmp <- dtp.1980.to.2100
  dtp.1980.to.2100.tmp[dtp.1980.to.2100<rcv.extra1] <- rcv.extra1[dtp.1980.to.2100<rcv.extra1]
  rcv1_bestcase <- new("cov.estimates",
                       sia.cov = sia.coverage.1980.2100,
                       sia.min.age = sia.min.age.1980.2100,
                       sia.max.age = sia.max.age.1980.2100,
                       MR1.cov = rcv1.coverage.1980.2100,
                       MR2.cov = rcv2.coverage.1980.2100,
                       inaccessible.prop = (1-dtp.1980.to.2100))



  #RCV2 DEFAULT SCENARIO  (this one has campaign)
  data <- read.csv("./data/montagu_coverage_data/202110gavi_v3/rubella-rcv2-default.csv")
  #unique(data$activity_type) #check activity types
  #unique(data$vaccine) #check vaccine types
  dfx.tmp <- data[data$vaccine=="Rubella" & data$activity_type=="routine" & data$coverage!=0,]
  #table(dfx.tmp$age_last) #check ages for RCV1
  if (any(unique(dfx.tmp$age_last)!=0)) print("RCV1 is not at 0 years old in rcv2_default")
  dfx.tmp2 <- data[data$vaccine=="RCV2" & data$activity_type=="routine" & data$coverage!=0,]
  #table(dfx.tmp2$age_last) #check ages for RCV2
  if (any(unique(dfx.tmp2$age_last)!=2)) print("RCV2 is not at 2 years old in rcv2_deafult")

  dfx <- data[data$activity_type=="campaign" & data$country_code==iso3code & data$coverage!=0,]
  sia.coverage.1980.2100 <- sia.min.age.1980.2100 <- sia.max.age.1980.2100 <- rep(0, (2101-1980))
  year.indeces <- as.numeric(dfx$year)-1980+1
  sia.coverage.1980.2100[year.indeces] <- as.numeric(dfx$coverage)
  sia.coverage.1980.2100 <- ifelse(sia.coverage.1980.2100>1, 1, sia.coverage.1980.2100) #cant have coverage over 1.
  sia.min.age.1980.2100[year.indeces] <- ifelse(dfx$age_first==0, 1, as.numeric(dfx$age_first*12))
  sia.max.age.1980.2100[year.indeces] <- as.numeric(dfx$age_last*12)
  if (nrow(dfx)==0){  #if no SIA, it needs to be in this form
    sia.coverage.1980.2100 <- sia.min.age.1980.2100 <- sia.max.age.1980.2100 <- rep(0, (2101-1980))
  }
  dfx <- data[data$vaccine=="Rubella" & data$activity_type=="routine" & data$country_code==iso3code & data$coverage!=0,]
  if (nrow(dfx)!=0){
    year.range <- range(dfx$year)
    if (any(diff(dfx$year)!=1)){
      rcv1.coverage <- c(dfx$coverage[1:which(diff(dfx$year)!=1)],rep(0,(diff(dfx$year)[which(diff(dfx$year)!=1)]-1)), dfx$coverage[(which(diff(dfx$year)!=1)+1):length(dfx$coverage)])
      rcv1.coverage.1980.2100 <- c(rep(0, (year.range[1]-1980)), rcv1.coverage, rep(0, (2100-year.range[2])))
    } else {
      rcv1.coverage <- as.numeric(dfx$coverage)
      rcv1.coverage.1980.2100 <- c(rep(0, (year.range[1]-1980)), rcv1.coverage, rep(0, (2100-year.range[2])))
    }
  } else {
    rcv1.coverage.1980.2100 <- rep(0, (2101-1980))
  }
  dfx <- data[data$vaccine=="RCV2" & data$activity_type=="routine" & data$country_code==iso3code & data$coverage!=0,]
  if (nrow(dfx)!=0){
    year.range <- range(dfx$year)
    if (any(diff(dfx$year)!=1)){
      rcv2.coverage <- c(dfx$coverage[1:which(diff(dfx$year)!=1)],rep(0,(diff(dfx$year)[which(diff(dfx$year)!=1)]-1)), dfx$coverage[(which(diff(dfx$year)!=1)+1):length(dfx$coverage)])
      rcv2.coverage.1980.2100 <- c(rep(0, (year.range[1]-1980)), rcv2.coverage, rep(0, (2100-year.range[2])))
    } else {
      rcv2.coverage <- as.numeric(dfx$coverage)
      rcv2.coverage.1980.2100 <- c(rep(0, (year.range[1]-1980)), rcv2.coverage, rep(0, (2100-year.range[2])))
    }
    rcv2.coverage.1980.2100[min(which(rcv1.coverage.1980.2100!=0))] <- 0 #the first year of RCV introduction cannot include RCV2! talk to VIMC
    max.rcv <- apply(cbind(rcv1.coverage.1980.2100, rcv2.coverage.1980.2100),1,max) #should be RCV1! forcing just in case
  } else {
    rcv2.coverage.1980.2100 <- rep(0, (2101-1980))
    max.rcv <- rcv1.coverage.1980.2100
  }
  #if dtp1 < RCV1 - force minimum acceptible to be RCV1
  rcv.extra1 <- c(max.rcv, max.rcv[length(max.rcv)])
  dtp.1980.to.2100.tmp <- dtp.1980.to.2100
  dtp.1980.to.2100.tmp[dtp.1980.to.2100<rcv.extra1] <- rcv.extra1[dtp.1980.to.2100<rcv.extra1]
  rcv2_default <- new("cov.estimates",
                      sia.cov = sia.coverage.1980.2100,
                      sia.min.age = sia.min.age.1980.2100,
                      sia.max.age = sia.max.age.1980.2100,
                      MR1.cov = rcv1.coverage.1980.2100,
                      MR2.cov = rcv2.coverage.1980.2100,
                      inaccessible.prop = (1-dtp.1980.to.2100))


  #RCV2 BESTCASE SCENARIO (this one has campaign)
  data <- read.csv("./data/montagu_coverage_data/202110gavi_v3/rubella-rcv2-ia2030_target.csv")
  #unique(data$activity_type) #check activity types
  #unique(data$vaccine) #check vaccine types
  dfx.tmp <- data[data$vaccine=="Rubella" & data$activity_type=="routine" & data$coverage!=0,]
  #table(dfx.tmp$age_last) #check ages for RCV1
  if (any(unique(dfx.tmp$age_last)!=0)) print("RCV1 is not at 0 years old in rcv2_default")
  dfx.tmp2 <- data[data$vaccine=="RCV2" & data$activity_type=="routine" & data$coverage!=0,]
  #table(dfx.tmp2$age_last) #check ages for RCV2
  if (any(unique(dfx.tmp2$age_last)!=2)) print("RCV2 is not at 2 years old in rcv2_deafult")

  dfx <- data[data$activity_type=="campaign" & data$country_code==iso3code & data$coverage!=0,]
  sia.coverage.1980.2100 <- sia.min.age.1980.2100 <- sia.max.age.1980.2100 <- rep(0, (2101-1980))
  year.indeces <- as.numeric(dfx$year)-1980+1
  sia.coverage.1980.2100[year.indeces] <- as.numeric(dfx$coverage)
  sia.coverage.1980.2100 <- ifelse(sia.coverage.1980.2100>1, 1, sia.coverage.1980.2100) #cant have coverage over 1.
  sia.min.age.1980.2100[year.indeces] <- ifelse(dfx$age_first==0, 1, as.numeric(dfx$age_first*12))
  sia.max.age.1980.2100[year.indeces] <- as.numeric(dfx$age_last*12)
  if (nrow(dfx)==0){  #if no SIA, it needs to be in this form
    sia.coverage.1980.2100 <- sia.min.age.1980.2100 <- sia.max.age.1980.2100 <- rep(0, (2101-1980))
  }
  dfx <- data[data$vaccine=="Rubella" & data$activity_type=="routine" & data$country_code==iso3code & data$coverage!=0,]
  if (nrow(dfx)!=0){
    year.range <- range(dfx$year)
    if (any(diff(dfx$year)!=1)){
      rcv1.coverage <- c(dfx$coverage[1:which(diff(dfx$year)!=1)],rep(0,(diff(dfx$year)[which(diff(dfx$year)!=1)]-1)), dfx$coverage[(which(diff(dfx$year)!=1)+1):length(dfx$coverage)])
      rcv1.coverage.1980.2100 <- c(rep(0, (year.range[1]-1980)), rcv1.coverage, rep(0, (2100-year.range[2])))
    } else {
      rcv1.coverage <- as.numeric(dfx$coverage)
      rcv1.coverage.1980.2100 <- c(rep(0, (year.range[1]-1980)), rcv1.coverage, rep(0, (2100-year.range[2])))
    }
  } else {
    rcv1.coverage.1980.2100 <- rep(0, (2101-1980))
  }
  dfx <- data[data$vaccine=="RCV2" & data$activity_type=="routine" & data$country_code==iso3code & data$coverage!=0,]
  if (nrow(dfx)!=0){
    year.range <- range(dfx$year)
    if (any(diff(dfx$year)!=1)){
      rcv2.coverage <- c(dfx$coverage[1:which(diff(dfx$year)!=1)],rep(0,(diff(dfx$year)[which(diff(dfx$year)!=1)]-1)), dfx$coverage[(which(diff(dfx$year)!=1)+1):length(dfx$coverage)])
      rcv2.coverage.1980.2100 <- c(rep(0, (year.range[1]-1980)), rcv2.coverage, rep(0, (2100-year.range[2])))
    } else {
      rcv2.coverage <- as.numeric(dfx$coverage)
      rcv2.coverage.1980.2100 <- c(rep(0, (year.range[1]-1980)), rcv2.coverage, rep(0, (2100-year.range[2])))
    }
    rcv2.coverage.1980.2100[min(which(rcv1.coverage.1980.2100!=0))] <- 0 #the first year of RCV introduction cannot include RCV2! talk to VIMC
    max.rcv <- apply(cbind(rcv1.coverage.1980.2100, rcv2.coverage.1980.2100),1,max) #should be RCV1! forcing just in case
  } else {
    rcv2.coverage.1980.2100 <- rep(0, (2101-1980))
    max.rcv <- rcv1.coverage.1980.2100
  }
  #if dtp1 < RCV1 - force minimum acceptible to be RCV1
  rcv.extra1 <- c(max.rcv, max.rcv[length(max.rcv)])
  dtp.1980.to.2100.tmp <- dtp.1980.to.2100
  dtp.1980.to.2100.tmp[dtp.1980.to.2100<rcv.extra1] <- rcv.extra1[dtp.1980.to.2100<rcv.extra1]
  rcv2_bestcase <- new("cov.estimates",
                       sia.cov = sia.coverage.1980.2100,
                       sia.min.age = sia.min.age.1980.2100,
                       sia.max.age = sia.max.age.1980.2100,
                       MR1.cov = rcv1.coverage.1980.2100,
                       MR2.cov = rcv2.coverage.1980.2100,
                       inaccessible.prop = (1-dtp.1980.to.2100))



  #RCV1 & RCV2 DEFAULT SCENARIO (this one has no campaign)
  data <- read.csv("./data/montagu_coverage_data/202110gavi_v3/rubella-rcv1-rcv2-default.csv")
  #unique(data$activity_type) #check activity types
  #unique(data$vaccine) #check vaccine types
  dfx.tmp <- data[data$vaccine=="Rubella" & data$activity_type=="routine" & data$coverage!=0,]
  #table(dfx.tmp$age_last) #check ages for RCV1
  if (any(unique(dfx.tmp$age_last)!=0)) print("RCV1 is not at 0 years old in rcv2_default")
  dfx.tmp2 <- data[data$vaccine=="RCV2" & data$activity_type=="routine" & data$coverage!=0,]
  #table(dfx.tmp2$age_last) #check ages for RCV2
  if (any(unique(dfx.tmp2$age_last)!=2)) print("RCV2 is not at 2 years old in rcv2_deafult")

  dfx <- data[data$activity_type=="campaign" & data$country_code==iso3code & data$coverage!=0,]
  sia.coverage.1980.2100 <- sia.min.age.1980.2100 <- sia.max.age.1980.2100 <- rep(0, (2101-1980))
  year.indeces <- as.numeric(dfx$year)-1980+1
  sia.coverage.1980.2100[year.indeces] <- as.numeric(dfx$coverage)
  sia.coverage.1980.2100 <- ifelse(sia.coverage.1980.2100>1, 1, sia.coverage.1980.2100) #cant have coverage over 1.
  sia.min.age.1980.2100[year.indeces] <- ifelse(dfx$age_first==0, 1, as.numeric(dfx$age_first*12))
  sia.max.age.1980.2100[year.indeces] <- as.numeric(dfx$age_last*12)
  if (nrow(dfx)==0){  #if no SIA, it needs to be in this form
    sia.coverage.1980.2100 <- sia.min.age.1980.2100 <- sia.max.age.1980.2100 <- rep(0, (2101-1980))
  }
  dfx <- data[data$vaccine=="Rubella" & data$activity_type=="routine" & data$country_code==iso3code & data$coverage!=0,]
  if (nrow(dfx)!=0){
    year.range <- range(dfx$year)
    if (any(diff(dfx$year)!=1)){
      rcv1.coverage <- c(dfx$coverage[1:which(diff(dfx$year)!=1)],rep(0,(diff(dfx$year)[which(diff(dfx$year)!=1)]-1)), dfx$coverage[(which(diff(dfx$year)!=1)+1):length(dfx$coverage)])
      rcv1.coverage.1980.2100 <- c(rep(0, (year.range[1]-1980)), rcv1.coverage, rep(0, (2100-year.range[2])))
    } else {
      rcv1.coverage <- as.numeric(dfx$coverage)
      rcv1.coverage.1980.2100 <- c(rep(0, (year.range[1]-1980)), rcv1.coverage, rep(0, (2100-year.range[2])))
    }
  } else {
    rcv1.coverage.1980.2100 <- rep(0, (2101-1980))
  }
  dfx <- data[data$vaccine=="RCV2" & data$activity_type=="routine" & data$country_code==iso3code & data$coverage!=0,]
  if (nrow(dfx)!=0){
    year.range <- range(dfx$year)
    if (any(diff(dfx$year)!=1)){
      rcv2.coverage <- c(dfx$coverage[1:which(diff(dfx$year)!=1)],rep(0,(diff(dfx$year)[which(diff(dfx$year)!=1)]-1)), dfx$coverage[(which(diff(dfx$year)!=1)+1):length(dfx$coverage)])
      rcv2.coverage.1980.2100 <- c(rep(0, (year.range[1]-1980)), rcv2.coverage, rep(0, (2100-year.range[2])))
    } else {
      rcv2.coverage <- as.numeric(dfx$coverage)
      rcv2.coverage.1980.2100 <- c(rep(0, (year.range[1]-1980)), rcv2.coverage, rep(0, (2100-year.range[2])))
    }
    rcv2.coverage.1980.2100[min(which(rcv1.coverage.1980.2100!=0))] <- 0 #the first year of RCV introduction cannot include RCV2! talk to VIMC
    max.rcv <- apply(cbind(rcv1.coverage.1980.2100, rcv2.coverage.1980.2100),1,max) #should be RCV1! forcing just in case
  } else {
    rcv2.coverage.1980.2100 <- rep(0, (2101-1980))
    max.rcv <- rcv1.coverage.1980.2100
  }
  #if dtp1 < RCV1 - force minimum acceptible to be RCV1
  rcv.extra1 <- c(max.rcv, max.rcv[length(max.rcv)])
  dtp.1980.to.2100.tmp <- dtp.1980.to.2100
  dtp.1980.to.2100.tmp[dtp.1980.to.2100<rcv.extra1] <- rcv.extra1[dtp.1980.to.2100<rcv.extra1]
  rcv1_rcv2_default <- new("cov.estimates",
                           sia.cov = sia.coverage.1980.2100,
                           sia.min.age = sia.min.age.1980.2100,
                           sia.max.age = sia.max.age.1980.2100,
                           MR1.cov = rcv1.coverage.1980.2100,
                           MR2.cov = rcv2.coverage.1980.2100,
                           inaccessible.prop = (1-dtp.1980.to.2100))



  #RCV1 & RCV2 BESTCASE SCENARIO (this one has no campaign)
  data <- read.csv("./data/montagu_coverage_data/202110gavi_v3/rubella-rcv1-rcv2-ia2030_target.csv")
  #unique(data$activity_type) #check activity types
  #unique(data$vaccine) #check vaccine types
  dfx.tmp <- data[data$vaccine=="Rubella" & data$activity_type=="routine" & data$coverage!=0,]
  #table(dfx.tmp$age_last) #check ages for RCV1
  if (any(unique(dfx.tmp$age_last)!=0)) print("RCV1 is not at 0 years old in rcv2_default")
  dfx.tmp2 <- data[data$vaccine=="RCV2" & data$activity_type=="routine" & data$coverage!=0,]
  #table(dfx.tmp2$age_last) #check ages for RCV2
  if (any(unique(dfx.tmp2$age_last)!=2)) print("RCV2 is not at 2 years old in rcv2_deafult")

  dfx <- data[data$activity_type=="campaign" & data$country_code==iso3code & data$coverage!=0,]
  sia.coverage.1980.2100 <- sia.min.age.1980.2100 <- sia.max.age.1980.2100 <- rep(0, (2101-1980))
  year.indeces <- as.numeric(dfx$year)-1980+1
  sia.coverage.1980.2100[year.indeces] <- as.numeric(dfx$coverage)
  sia.coverage.1980.2100 <- ifelse(sia.coverage.1980.2100>1, 1, sia.coverage.1980.2100) #cant have coverage over 1.
  sia.min.age.1980.2100[year.indeces] <- ifelse(dfx$age_first==0, 1, as.numeric(dfx$age_first*12))
  sia.max.age.1980.2100[year.indeces] <- as.numeric(dfx$age_last*12)
  if (nrow(dfx)==0){  #if no SIA, it needs to be in this form
    sia.coverage.1980.2100 <- sia.min.age.1980.2100 <- sia.max.age.1980.2100 <- rep(0, (2101-1980))
  }
  dfx <- data[data$vaccine=="Rubella" & data$activity_type=="routine" & data$country_code==iso3code & data$coverage!=0,]
  if (nrow(dfx)!=0){
    year.range <- range(dfx$year)
    if (any(diff(dfx$year)!=1)){
      rcv1.coverage <- c(dfx$coverage[1:which(diff(dfx$year)!=1)],rep(0,(diff(dfx$year)[which(diff(dfx$year)!=1)]-1)), dfx$coverage[(which(diff(dfx$year)!=1)+1):length(dfx$coverage)])
      rcv1.coverage.1980.2100 <- c(rep(0, (year.range[1]-1980)), rcv1.coverage, rep(0, (2100-year.range[2])))
    } else {
      rcv1.coverage <- as.numeric(dfx$coverage)
      rcv1.coverage.1980.2100 <- c(rep(0, (year.range[1]-1980)), rcv1.coverage, rep(0, (2100-year.range[2])))
    }
  } else {
    rcv1.coverage.1980.2100 <- rep(0, (2101-1980))
  }
  dfx <- data[data$vaccine=="RCV2" & data$activity_type=="routine" & data$country_code==iso3code & data$coverage!=0,]
  if (nrow(dfx)!=0){
    year.range <- range(dfx$year)
    if (any(diff(dfx$year)!=1)){
      rcv2.coverage <- c(dfx$coverage[1:which(diff(dfx$year)!=1)],rep(0,(diff(dfx$year)[which(diff(dfx$year)!=1)]-1)), dfx$coverage[(which(diff(dfx$year)!=1)+1):length(dfx$coverage)])
      rcv2.coverage.1980.2100 <- c(rep(0, (year.range[1]-1980)), rcv2.coverage, rep(0, (2100-year.range[2])))
    } else {
      rcv2.coverage <- as.numeric(dfx$coverage)
      rcv2.coverage.1980.2100 <- c(rep(0, (year.range[1]-1980)), rcv2.coverage, rep(0, (2100-year.range[2])))
    }
    rcv2.coverage.1980.2100[min(which(rcv1.coverage.1980.2100!=0))] <- 0 #the first year of RCV introduction cannot include RCV2! talk to VIMC
    max.rcv <- apply(cbind(rcv1.coverage.1980.2100, rcv2.coverage.1980.2100),1,max) #should be RCV1! forcing just in case
  } else {
    rcv2.coverage.1980.2100 <- rep(0, (2101-1980))
    max.rcv <- rcv1.coverage.1980.2100
  }
  #if dtp1 < RCV1 - force minimum acceptible to be RCV1
  rcv.extra1 <- c(max.rcv, max.rcv[length(max.rcv)])
  dtp.1980.to.2100.tmp <- dtp.1980.to.2100
  dtp.1980.to.2100.tmp[dtp.1980.to.2100<rcv.extra1] <- rcv.extra1[dtp.1980.to.2100<rcv.extra1]
  rcv1_rcv2_bestcase <- new("cov.estimates",
                            sia.cov = sia.coverage.1980.2100,
                            sia.min.age = sia.min.age.1980.2100,
                            sia.max.age = sia.max.age.1980.2100,
                            MR1.cov = rcv1.coverage.1980.2100,
                            MR2.cov = rcv2.coverage.1980.2100,
                            inaccessible.prop = (1-dtp.1980.to.2100))


  #NO VACCINATION SCENARIO
  data <- read.csv("./data/montagu_coverage_data/202110gavi_v3/rubella-routine-no-vaccination.csv")
  #unique(data$activity_type) #check activity types
  #unique(data$vaccine) #check vaccine types
  #dfx <- data[data$vaccine=="Rubella" & data$activity_type=="routine" & data$coverage!=0,]
  #table(dfx$age_last) #check ages
  #dfx <- data[data$vaccine=="RCV2" & data$activity_type=="routine" & data$coverage!=0,]
  #table(dfx$age_last) #check ages
  dfx <- data[data$activity_type=="campaign" & data$country_code==iso3code & data$coverage!=0,]
  sia.coverage.1980.2100 <- sia.min.age.1980.2100 <- sia.max.age.1980.2100 <- rep(0, (2101-1980))
  year.indeces <- as.numeric(dfx$year)-1980+1
  sia.coverage.1980.2100[year.indeces] <- as.numeric(dfx$coverage)
  sia.coverage.1980.2100 <- ifelse(sia.coverage.1980.2100>1, 1, sia.coverage.1980.2100) #cant have coverage over 1.
  sia.min.age.1980.2100[year.indeces] <- ifelse(dfx$age_first==0, 1, as.numeric(dfx$age_first*12))
  sia.max.age.1980.2100[year.indeces] <- as.numeric(dfx$age_last*12)
  if (nrow(dfx)==0){  #if no SIA, it needs to be in this form
    sia.coverage.1980.2100 <- sia.min.age.1980.2100 <- sia.max.age.1980.2100 <- rep(0, (2101-1980))
  }
  dfx <- data[data$vaccine=="Rubella" & data$activity_type=="routine" & data$country_code==iso3code & data$coverage!=0,]
  if (nrow(dfx)!=0){
    year.range <- range(dfx$year)
    if (any(diff(dfx$year)!=1)){
      rcv1.coverage <- c(dfx$coverage[1:which(diff(dfx$year)!=1)],rep(0,(diff(dfx$year)[which(diff(dfx$year)!=1)]-1)), dfx$coverage[(which(diff(dfx$year)!=1)+1):length(dfx$coverage)])
      rcv1.coverage.1980.2100 <- c(rep(0, (year.range[1]-1980)), rcv1.coverage, rep(0, (2100-year.range[2])))
    } else {
      rcv1.coverage <- as.numeric(dfx$coverage)
      rcv1.coverage.1980.2100 <- c(rep(0, (year.range[1]-1980)), rcv1.coverage, rep(0, (2100-year.range[2])))
    }
  } else {
    rcv1.coverage.1980.2100 <- rep(0, (2101-1980))
  }
  dfx <- data[data$vaccine=="RCV2" & data$activity_type=="routine" & data$country_code==iso3code & data$coverage!=0,]
  if (nrow(dfx)!=0){
    year.range <- range(dfx$year)
    if (any(diff(dfx$year)!=1)){
      rcv2.coverage <- c(dfx$coverage[1:which(diff(dfx$year)!=1)],rep(0,(diff(dfx$year)[which(diff(dfx$year)!=1)]-1)), dfx$coverage[(which(diff(dfx$year)!=1)+1):length(dfx$coverage)])
      rcv2.coverage.1980.2100 <- c(rep(0, (year.range[1]-1980)), rcv2.coverage, rep(0, (2100-year.range[2])))
    } else {
      rcv2.coverage <- as.numeric(dfx$coverage)
      rcv2.coverage.1980.2100 <- c(rep(0, (year.range[1]-1980)), rcv2.coverage, rep(0, (2100-year.range[2])))
    }
    rcv2.coverage.1980.2100[min(which(rcv1.coverage.1980.2100!=0))] <- 0 #the first year of RCV introduction cannot include RCV2! talk to VIMC
    max.rcv <- apply(cbind(rcv1.coverage.1980.2100, rcv2.coverage.1980.2100),1,max) #should be RCV1! forcing just in case
  } else {
    rcv2.coverage.1980.2100 <- rep(0, (2101-1980))
    max.rcv <- rcv1.coverage.1980.2100
  }
  #if dtp1 < RCV1 - force minimum acceptible to be RCV1
  rcv.extra1 <- c(max.rcv, max.rcv[length(max.rcv)])
  dtp.1980.to.2100.tmp <- dtp.1980.to.2100
  dtp.1980.to.2100.tmp[dtp.1980.to.2100<rcv.extra1] <- rcv.extra1[dtp.1980.to.2100<rcv.extra1]
  no_vaccination <- new("cov.estimates",
                        sia.cov = sia.coverage.1980.2100,
                        sia.min.age = sia.min.age.1980.2100,
                        sia.max.age = sia.max.age.1980.2100,
                        MR1.cov = rcv1.coverage.1980.2100,
                        MR2.cov = rcv2.coverage.1980.2100,
                        inaccessible.prop = (1-dtp.1980.to.2100))


  return(list(campaign_default=campaign_default,
              campaign_bestcase=campaign_bestcase,
              rcv1_default=rcv1_default,
              rcv1_bestcase=rcv1_bestcase,
              rcv2_default=rcv2_default,
              rcv2_bestcase=rcv2_bestcase,
              rcv1_rcv2_default=rcv1_rcv2_default,
              rcv1_rcv2_bestcase=rcv1_rcv2_bestcase,
              no_vaccination = no_vaccination))
}


### Function to get all the things you need to run
### setup coverage and demography for country (Montagu for touchstone 202110gavi_v3)
###
##'  @param - country iso3code
###
### Returns demography and coverage for each country
setupCountry_202110gavi_v3 <- function(country){

  iso3code <- country #correcting for XK b/c no match
  if (iso3code=="XK") {
    uncode <- 999
  } else {
    uncode <- as.numeric(countrycode::countrycode(iso3code, origin="iso3c", destination="un"))
  }
  if (is.na(uncode)) stop("uncode not matched to iso3code, see setupCountry_202110gavi_v3()")

  #Demography
  demog <- getDemography(uncode=uncode)
  pop.total.1950.2100 <- demog$pop.total.1950.2100*1000
  pop.age.byageclasses.1950.2100 <- demog$pop.age.byageclasses.1950.2100*1000
  asfr.1950.2100 <- demog$asfr.1950.2100
  e0.1950.2100 <- demog$e0.1950.2100
  fert <- c(0,asfr.1950.2100[,(2015-1950+1)]*1000,0)
  cbr.1950.2100 <- demog$cbr.1950.2100
  asdr.1950.2100.by5 <- demog$asdr.1950.2100.by5
  asdr.object <- demog$asdr.object
  repro.age.sex.dist.1950.2100 <- demog$repro.age.sex.dist.1950.2100

  #Coverage estiamtes from Montagu, and inaccessible population
  coverage <- getMontaguCoverage_202110gavi_v3(iso3code)

  #set up the births function
  get.births.here <- function(state,tran,
                              fert.curve=fert,
                              lower.age.boundary = c(0,15,20,25,30,35,40,45,50,101)) { #xxamy - age.class change from lower.age.boundary = c(0,15,20,25,30,35,40,45,49,100)) {
    return(get.births.fun(state=state,tran=tran,fert.curve=fert.curve,lower.age.boundary=lower.age.boundary))

  }

  return(list(uncode=uncode,
              iso3code=iso3code,
              pop.total.1950.2100=pop.total.1950.2100,
              pop.age.byageclasses.1950.2100=pop.age.byageclasses.1950.2100,
              e0.1950.2100=e0.1950.2100,
              asfr.1950.2100=asfr.1950.2100,
              repro.age.sex.dist.1950.2100=repro.age.sex.dist.1950.2100,
              cbr.1950.2100=cbr.1950.2100,
              yr.agespecificbirths.per.1000=fert,
              asdr.1950.2100.by5=asdr.1950.2100.by5,
              asdr.object=asdr.object,
              get.births.here=get.births.here,
              coverage = coverage))
}
