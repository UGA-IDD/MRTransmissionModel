## Functions to create data inputs based on WHO vaccination data and UNPD demography data
## Authors: Amy Winter

### Function to get all the things you need to run
### setup coverage and demography for country (NOT Montagu)
###
##'  @param - country - character - country name 
###
### Returns demography and coverage for each country
setupCountry.Dec2021 <- function(country){ 
  
  uncode <- countrycode::countrycode(country, origin="country.name", destination="un")
  iso3code <- countrycode::countrycode(country, origin="country.name", destination="iso3c")
  
  demog <- getDemography.wpp2019(uncode=uncode)
  pop.total.1950.2100 <- demog$pop.total.1950.2100*1000
  pop.age.byageclasses.1950.2100 <- demog$pop.age.byageclasses.1950.2100*1000
  tfr.1950.2100 <- demog$tfr.1950.2100
  asfr.1950.2100 <- demog$asfr.1950.2100
  e0.1950.2100 <- demog$e0.1950.2100
  fert <- c(0,asfr.1950.2100[,(2015-1950+1)]*1000,0)
  births.1950.2100 <- demog$births.1950.2100
  cbr.1950.2100 <- demog$cbr.1950.2100
  asdr.1950.2100.by5 <- demog$asdr.1950.2100.by5
  asdr.object <- demog$asdr.object
  repro.age.sex.dist.1950.2100 <- demog$repro.age.sex.dist.1950.2100
  
  # MCV routine coverage
  df <- read_xlsx("./data/MCV_WUENIC_coverage_estimates_DownloadedDec2021.xlsx", sheet="Sheet1")
  MCV1.coverage <- df %>%
    dplyr::filter(COVERAGE_CATEGORY=="WUENIC",
                  CODE==iso3code,
                  ANTIGEN=="MCV1") %>%
    dplyr::select(COVERAGE) %>% 
    dplyr::mutate(COVERAGE=COVERAGE/99) %>%
    unlist() %>% as.numeric %>% rev
  MCV1.coverage <- c(MCV1.coverage, rep(mean(MCV1.coverage[39:41]), length(2021:2100)))
  MCV1.coverage[is.na(MCV1.coverage)] <- 0 #is NA assume 0
  MCV2.coverage <- df %>%
    dplyr::filter(COVERAGE_CATEGORY=="WUENIC",
                  CODE==iso3code,
                  ANTIGEN=="MCV2") %>%
    dplyr::select(COVERAGE) %>% 
    dplyr::mutate(COVERAGE=COVERAGE/99) %>%
    unlist() %>% as.numeric %>% rev
  MCV2.coverage <- c(rep(0, length(1980:1999)), MCV2.coverage, rep(mean(MCV2.coverage[19:21], na.rm=T), length(2021:2100)))
  MCV2.coverage[is.na(MCV2.coverage)] <- 0 #is NA assume 0
  
  # RCV routine Coverage
  df <- read_xlsx("./data/RCV_WUENIC_coverage_estimates_DownloadedDec2021.xlsx", sheet="Sheet1")
  RCV1.coverage <- df %>%
    dplyr::filter(COVERAGE_CATEGORY=="WUENIC",
                  CODE==iso3code,
                  ANTIGEN=="RCV1") %>%
    dplyr::select(COVERAGE) %>% 
    dplyr::mutate(COVERAGE=COVERAGE/99) %>%
    unlist() %>% as.numeric %>% rev
  RCV1.coverage <- c(RCV1.coverage, rep(mean(RCV1.coverage[39:41]), length(2021:2100)))
  RCV1.coverage[is.na(RCV1.coverage)] <- 0 #is NA assume 0
  
  # SIA Coverage - FOR MEASLES
  SIA.measles <- getWHO.measles.SIAs(iso3code, uncode, upper.limit=0.97) ##xxamy - need to make sure this works for all countries
  measlesSIA.coverage <- age.min.sia.measles <- age.max.sia.measles <- rep(0, length(MCV1.coverage))
  measlesSIA.coverage[c(1980:2100) %in% SIA.measles$year] <- SIA.measles$coverage
  age.min.sia.measles[c(1980:2100) %in% SIA.measles$year] <- SIA.measles$age.min
  age.max.sia.measles[c(1980:2100) %in% SIA.measles$year] <- SIA.measles$age.max
  
  #SIA Coverage - For RUBELLA - may need updating
  SIA.rubella <- getWHO.rubella.SIAs(iso3code, uncode, upper.limit=0.97)
  rubellaSIA.coverage <- age.min.sia.rubella <- age.max.sia.rubella <- rep(0, length(MCV1.coverage))
  rubellaSIA.coverage[c(1980:2100) %in% SIA.rubella$year] <- SIA.rubella$coverage
  age.min.sia.rubella[c(1980:2100) %in% SIA.rubella$year] <- SIA.rubella$age.min
  age.max.sia.rubella[c(1980:2100) %in% SIA.rubella$year] <- SIA.rubella$age.max
  
  #set up the births function 
  get.births.here <- function(state,tran,
                              fert.curve=fert,
                              lower.age.boundary = c(0,15,20,25,30,35,40,45,50,101)) { #xxamy - age.class change from lower.age.boundary = c(0,15,20,25,30,35,40,45,49,100)) {
    return(get.births.fun(state=state,tran=tran,fert.curve=fert.curve,lower.age.boundary=lower.age.boundary))
    
  }
  
  #Inaccessible from DPT1 numbers over time
  dtp <- read.csv("./data/dtp1_estimates.csv")
  dtp.c <- dtp[dtp$ISO_code==iso3code,]
  dtp.c <- dtp.c[order(dtp.c$Year),]
  dtp.1980.to.2100 <- dtp.c$vacc_coverage[1:length(1980:2017)]
  #filling in zeros using zoo package
  dtp.1980.to.2100 <- zoo::na.locf(dtp.1980.to.2100, fromLast=T, na.rm=F) #fills in using value after NA
  dtp.1980.to.2100 <- zoo::na.locf(dtp.1980.to.2100) #fills in using previous value before NA
  dtp.1980.to.2100 <- c(dtp.1980.to.2100, rep(dtp.1980.to.2100[length(dtp.1980.to.2100)], length(2018:2101))) #add one more year 1980to2101
  
  if (length(dtp.1980.to.2100)==0){ #Palestine and Kosovo dont have DTP1 estimates
    dtp.1980.to.2100 <- rep(1, (length(RCV1.coverage.1980to2100)+1))
  }
  
  #Force minimum acceptible to be MCV1, if dtp1 < MCV1 - xxxamy
  mcv.extra1 <- c(MCV1.coverage, MCV1.coverage[length(MCV1.coverage)])
  dtp.1980.to.2100[dtp.1980.to.2100<mcv.extra1] <- mcv.extra1[dtp.1980.to.2100<mcv.extra1]
  
  
  return(list(uncode=uncode,
              iso3code=iso3code,
              pop.total.1950.2100=pop.total.1950.2100,
              pop.age.byageclasses.1950.2100=pop.age.byageclasses.1950.2100,
              tfr.1950.2100=tfr.1950.2100,
              e0.1950.2100=e0.1950.2100,
              asfr.1950.2100=asfr.1950.2100,
              repro.age.sex.dist.1950.2100=repro.age.sex.dist.1950.2100,
              births.1950.2100=births.1950.2100,
              cbr.1950.2100=cbr.1950.2100,
              yr.agespecificbirths.per.1000=fert,
              asdr.1950.2100.by5=asdr.1950.2100.by5,
              asdr.object=asdr.object,
              get.births.here=get.births.here,
              MCV1.coverage.1980to2100=MCV1.coverage,
              MCV2.coverage.1980to2100=MCV2.coverage,
              RCV1.coverage.1980to2100=RCV1.coverage,
              measlesSIA.coverage.1980to2100=measlesSIA.coverage,
              age.max.sia.measles=age.max.sia.measles,
              age.min.sia.measles=age.min.sia.measles,
              rubellaSIA.coverage.1980to2100=rubellaSIA.coverage,
              age.max.sia.rubella=age.max.sia.rubella,
              age.min.sia.rubella=age.min.sia.rubella,
              inaccessible.prop.1980.to.2100=(1-dtp.1980.to.2100)))
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


### Function to convert measles SIAs into something usable
###
##' @param - iso3code
##' @param - uncode
##' @param - upper.limit = upper coverage limit as proportion
###
### Returns age range, year, and coverage of each SIA
getWHO.measles.SIAs <- function(iso3code, uncode, upper.limit=0.97){
  
  df <- read_excel("./data/SIAs_WHO_DownloadedDec2021.xls", sheet="SIAs_2000_2022")
  df <- df %>% 
    filter(Intervention !="Rubella") %>% 
    filter(Activity=="Campaign" | Activity=="CatchUp" | Activity=="FollowUp" | Activity=="SIA") %>%
    filter(Extent!="unknown") %>%
    filter(`Implementation status` =="done") %>%
    filter(!is.na(Year),
           !is.na(`Age group`))
  
  mtch <- which(df$`ISO3 code`==iso3code,arr.ind=TRUE)
  df <- df[mtch,]
  
  if (nrow(df)!=0){
    
    age.range.new <- convertAgeSIA(df$`Age group`)
    years.sia.new <- df$Year
    
    coverage.new <- rep(NA, length(years.sia.new))
    for (s in 1:length(years.sia.new)){
      age.min.months <- age.range.new[s,1]
      age.min <- ceiling(age.range.new[s,1]/12)
      age.max <- ceiling(age.range.new[s,2]/12)
      y <- years.sia.new[s]
      
      if (df$Extent[s]=="national" | df$Extent[s]=="National"){
        if (!is.na(df$`Survey results (%)`[s])){
          coverage.new[s] <- (df$`Survey results (%)`[s])/100
        } else {
          coverage.new[s] <- (df$`% Reached`[s])/100
        }
      } else { #Extent is Rollover-National, Sub-national, or Sub-National
        #get national population of the target age
        out <- getDemography.wpp2019(uncode=uncode)
        if (age.min==1){
          tot.pop <- sum(out$pop.age.byageclasses.1950.2100[2:age.max,(y-1950)])*1000
          tot.pop <- as.numeric(tot.pop+(age.min.months/12*out$pop.age.byageclasses.1950.2100[1,(y-1950)]*1000))
        } else {
          tot.pop <- as.numeric(sum(out$pop.age.byageclasses.1950.2100[age.min:age.max,(y-1950)])*1000)
        }
        #multiply coverage by the targeted age then divide by national target age group to get national coverage
        if (!is.na(df$`Survey results (%)`[s])){
          coverage.new[s] <- (df$`Survey results (%)`[s])/100 * (df$`Target population`[s]) / tot.pop
        } else {
          coverage.new[s] <- (df$`Reached population`[s]) / tot.pop
        }
      }
    }
    
    coverage.new <- (pmin(coverage.new, upper.limit))
    age.range.new.min <- age.range.new[,1]
    age.range.new.max <- age.range.new[,2]
    
  } else {
    coverage.new <- 0
    age.range.new.min <- 12
    age.range.new.max <- 48
    years.sia.new <- 2000
  }
  
  return(list(age.min=age.range.new.min, age.max=age.range.new.max, year=years.sia.new, coverage=coverage.new))
}

### Function to convert rubella SIAs into something usable
###
##' @param - iso3code
##' @param - uncode
##' ##' @param - upper.limit = upper coverage limit as proportion
###
### Returns age range, year, and coverage of each SIA
getWHO.rubella.SIAs <- function(iso3code, uncode, upper.limit=0.97){
  
  df <- read_excel("./data/SIAs_WHO_DownloadedDec2021.xls", sheet="SIAs_2000_2022")
  df <- df %>% 
    filter(Intervention !="MCV", Intervention !="measles", Intervention !="Measles") %>% 
    filter(Activity=="Campaign" | Activity=="CatchUp" | Activity=="FollowUp" | Activity=="SIA") %>%
    filter(Extent!="unknown") %>%
    filter(`Implementation status` =="done") %>%
    filter(!is.na(Year),
           !is.na(`Age group`))
  
  mtch <- which(df$`ISO3 code`==iso3code,arr.ind=TRUE)
  df <- df[mtch,]
  
  if (nrow(df)!=0){
    
    age.range.new <- convertAgeSIA(df$`Age group`)
    years.sia.new <- df$Year
    
    coverage.new <- rep(NA, length(years.sia.new))
    for (s in 1:length(years.sia.new)){
      age.min.months <- age.range.new[s,1]
      age.min <- ceiling(age.range.new[s,1]/12)
      age.max <- ceiling(age.range.new[s,2]/12)
      y <- years.sia.new[s]
      
      if (df$Extent[s]=="national" | df$Extent[s]=="National"){
        if (!is.na(df$`Survey results (%)`[s])){
          coverage.new[s] <- (df$`Survey results (%)`[s])/100
        } else {
          coverage.new[s] <- (df$`% Reached`[s])/100
        }
      } else { #Extent is Rollover-National, Sub-national, or Sub-National
        #get national population of the target age
        out <- getDemography.wpp2019(uncode=uncode)
        if (age.min==1){
          tot.pop <- sum(out$pop.age.byageclasses.1950.2100[2:age.max,(y-1950)])*1000
          tot.pop <- as.numeric(tot.pop+(age.min.months/12*out$pop.age.byageclasses.1950.2100[1,(y-1950)]*1000))
        } else {
          tot.pop <- as.numeric(sum(out$pop.age.byageclasses.1950.2100[age.min:age.max,(y-1950)])*1000)
        }
        #multiply coverage by the targeted age then divide by national target age group to get national coverage
        if (!is.na(df$`Survey results (%)`[s])){
          coverage.new[s] <- (df$`Survey results (%)`[s])/100 * (df$`Target population`[s]) / tot.pop
        } else {
          coverage.new[s] <- (df$`Reached population`[s]) / tot.pop
        }
      }
    }
    
    coverage.new <- (pmin(coverage.new, upper.limit))
    age.range.new.min <- age.range.new[,1]
    age.range.new.max <- age.range.new[,2]
    
  } else {
    coverage.new <- 0
    age.range.new.min <- 12
    age.range.new.max <- 48
    years.sia.new <- 2000
  }
  
  return(list(age.min=age.range.new.min, age.max=age.range.new.max, year=years.sia.new, coverage=coverage.new))
}

### Function to convert the SIA age range into something usable
###
##' @param - sia.age.range - input as AgeGroup from the WHO SIA Data
###
### Returns age range
convertAgeSIA <- function(sia.age.range){
  
  age.lower <- rep(NA,length(sia.age.range))
  age.lower.unit <- rep(NA,length(sia.age.range))
  age.upper <- rep(NA,length(sia.age.range))
  age.upper.unit <- rep(NA,length(sia.age.range))
  
  sia.age.range[sia.age.range=="unknown"] <- "9-59 M"               # assume the unkowns are the usual 9mo-5yrs
  sia.age.range[sia.age.range=="Children at elementary"] <- "5-11 Y"
  sia.age.range[sia.age.range=="School-age"] <- "5-11 Y"
  sia.age.range[sia.age.range=="eligible children"] <- "9-59 M"
  sia.age.range[sia.age.range=="6 M+"] <- "6 M-15 Y"
  sia.age.range[sia.age.range=="<15 Y"] <- "9 M-14 Y"
  sia.age.range[sia.age.range=="<5 Y"] <- "9-59 M"
  sia.age.range[sia.age.range=="<4 Y"] <- "9-47 M"
  sia.age.range[sia.age.range=="5 Y"] <- "9-59 M"
  sia.age.range[sia.age.range=="14 Y"] <- "9 M-14 Y"
  sia.age.range[sia.age.range=="13 Y"] <- "9 M-13 Y"
  sia.age.range[sia.age.range=="<1 Y"] <- "6-12 M"
  sia.age.range[sia.age.range=="12 M"] <- "6-12 M"
  sia.age.range[sia.age.range=="1 Y"] <- "6-12 M"
  sia.age.range[sia.age.range=="4 Y"] <- "9-47 M"
  sia.age.range[sia.age.range=="9-5 Y"] <- "9-59 M"
  sia.age.range[sia.age.range=="1 Y school"] <- "5-6 Y"
  sia.age.range[sia.age.range=="1st year primary school"] <- "5-6 Y"
  
  sia.age.tmp <- rep(NA, length(sia.age.range))
  sia.age.tmp[grepl("-", sia.age.range)] <- as.vector(sia.age.range[grepl("-", sia.age.range)])
  
  for (i in 1:length(sia.age.tmp)){
    age.mat <- matrix(unlist(strsplit(sia.age.tmp[i], "-")), ncol=2, byrow=TRUE)
    # Define Lower Age
    if (grepl("M", age.mat[,1])){
      age.lower.unit[i] <- "M"
      age.lower[i] <- gsub("M", "", age.mat[,1])
      age.lower[i] <- gsub(" ", "", age.lower[i])
    } else if (grepl("Y", age.mat[,1])){
      age.lower.unit[i] <- "Y"
      age.lower[i] <- gsub("Y", "", age.mat[,1])
      age.lower[i] <- gsub(" ", "", age.lower[i])
    }
    # Define Upper Age
    if (grepl("M", age.mat[,2])){
      age.upper.unit[i] <- "M"
      age.upper[i] <- gsub("M", "", age.mat[,2])
      age.upper[i] <- gsub(" ", "", age.upper[i])
      if (grepl("<", age.upper[i])){
        age.upper[i] <- gsub("<", "", age.upper[i])
        age.upper[i] <- as.numeric(age.upper[i]) - 1
      }
    } else if (grepl("Y", age.mat[,2])){
      age.upper.unit[i] <- "Y"
      age.upper[i] <- gsub("Y", "", age.mat[,2])
      age.upper[i] <- gsub(" ", "", age.upper[i])
      if (grepl("<", age.upper[i])){
        age.upper[i] <- gsub("<", "", age.upper[i])
        age.upper[i] <- as.numeric(age.upper[i]) - 1
      }
    }
    
    if (grepl("M", age.mat[,1])==FALSE & grepl("Y", age.mat[,1])==FALSE){
      age.lower[i] <- gsub(" ", "", age.mat[,1])
      age.lower.unit[i] <- age.upper.unit[i]
    }
  }
  
  age.lower.res <- as.numeric(age.lower)
  age.upper.res <- as.numeric(age.upper)
  age.lower.res[age.lower.unit=="Y" & !is.na(age.lower.unit)] <- 12 * age.lower.res[age.lower.unit=="Y" & !is.na(age.lower.unit)]
  age.upper.res[age.upper.unit=="Y" & !is.na(age.upper.unit)] <- 12 * age.upper.res[age.upper.unit=="Y" & !is.na(age.upper.unit)]
  
  age.range <- as.matrix(cbind(age.lower.res, age.upper.res))
  
  return(age.range)
}


