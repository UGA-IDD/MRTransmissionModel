#' Function to pull in vaccination data from Montagu 202110sgavi_v3 touchstone vaccine scenarios
#'
#' @param iso3code iso3code
#'
#' @importFrom utils read.csv
#'
#' @return age, year, and coverage for SIA, RCV1, and RCV2 for each vaccination scenario
#' @export
#'

getMontaguCoverage_202110gavi_v3 <- function(iso3code){

  #Also Estimating inaccessible at the same time, because it depends on routine estimates
  #Inaccessible from DPT1 numbers over time
  dtp <- readRDS("vimc_private/dtp1_estimates.RDS")
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
  data <- read.csv("vimc_private/montagu_coverage_data/202110gavi_v3/rubella-campaign-default.csv")
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
  data <- read.csv("vimc_private/montagu_coverage_data/202110gavi_v3/rubella-campaign-ia2030_target.csv")
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
  data <- read.csv("vimc_private/montagu_coverage_data/202110gavi_v3/rubella-rcv1-default.csv")
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
  data <- read.csv("vimc_private/montagu_coverage_data/202110gavi_v3/rubella-rcv1-ia2030_target.csv")
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
  data <- read.csv("vimc_private/montagu_coverage_data/202110gavi_v3/rubella-rcv2-default.csv")
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
  data <- read.csv("vimc_private/montagu_coverage_data/202110gavi_v3/rubella-rcv2-ia2030_target.csv")
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
  data <- read.csv("vimc_private/montagu_coverage_data/202110gavi_v3/rubella-rcv1-rcv2-default.csv")
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
  data <- read.csv("vimc_private/montagu_coverage_data/202110gavi_v3/rubella-rcv1-rcv2-ia2030_target.csv")
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
  data <- read.csv("vimc_private/montagu_coverage_data/202110gavi_v3/rubella-routine-no-vaccination.csv")
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
