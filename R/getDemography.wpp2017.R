#' Function to get population size, age structure, and fertility from 1950 to 2105
#'
#' @param uncode UN country code
#' @param age.classes age classes in years
#' @param if.missing.use.region logical; c(MHL, XK, TUV)
#'
#'
#' @importFrom stats smooth.spline predict
#' @importFrom zoo na.spline
#' @importFrom utils data
#' @import countrycode
#' @import wpp2017
#'
#' @return demography for each country
#' @export
#'

getDemography.wpp2017 <- function(uncode=NA, age.classes=c(1:101),
                                  if.missing.use.region=F){

  # deeply buried in model and won't change
  time.by5 <- seq(1950,2100,5)
  time <- seq(1950,2100,1)

  #UNPD data
  # data automatically loaded when install package
  #data("pop_wpp2017")
  #pop <- pop_wpp2017
  data("pop", package = "wpp2017", envir = environment())
  # wpp2017::pop, data from the wpp2017 R package, user will not change

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
    if (uncode == 212) { #Dominica
      region <- "LATIN AMERICA AND THE CARIBBEAN"
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
  #data(pop) # wpp2017
  #data("popproj_wpp2017")
  #popproj <- popproj_wpp2017
  data("popproj", package = "wpp2017", envir = environment()) # wpp2017

  pop.total.1950.2100.by5 <-  as.numeric(cbind(pop[pop$country_code==cc,3:ncol(pop)],
                                               popproj[popproj$country_code==cc,3:ncol(popproj)]))
  f <- smooth.spline(time.by5, pop.total.1950.2100.by5)
  pop.total.1950.2100 <- predict(f,time)$y

  #Population Age Structure by Age over time
  #data("popFprojMed_wpp2017")
  #popFprojMed <- popFprojMed_wpp2017

  #data("popMprojMed_wpp2017")
  #popMprojMed <- popMprojMed_wpp2017

  data("popF", package = "wpp2017", envir = environment()) # wpp2017
  data("popM", package = "wpp2017", envir = environment()) # wpp2017
  data("popFprojMed", package = "wpp2017", envir = environment()) # wpp2017
  data("popMprojMed", package = "wpp2017", envir = environment()) # wpp2017

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
  #data("tfr_wpp2017")
  #tfr <- tfr_wpp2017

  #data("tfrprojMed_wpp2017")
  #tfrprojMed <- tfrprojMed_wpp2017

  data("tfr", package = "wpp2017", envir = environment()) # wpp2017
  data("tfrprojMed", package = "wpp2017", envir = environment()) # wpp2017

  mid.time.by5 <- seq(1952, 2097, 5) #TFR is given over a range, therefore assume it is the mid-period
  tfr.1950.2100.by5 <- as.numeric(cbind(tfr[tfr$country_code==cc,3:15],
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
  #data("percentASFR_wpp2017")
  #percentASFR <- percentASFR_wpp2017

  data("percentASFR", package = "wpp2017", envir = environment()) # wpp2017

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
  data("mxF", package = "wpp2017", envir = environment()) # wpp2017

  #data("mxF_wpp2017")
  #mxF <- mxF_wpp2017

  asdr.1950.2100.by5 <- (mxF[mxF$country_code==cc,4:ncol(mxF)])
  asdr.1950.2100.by5 <- cbind(asdr.1950.2100.by5, "2100-2105" = asdr.1950.2100.by5[,ncol(asdr.1950.2100.by5)]) #revised to 2105 xxamy
  asdr.maxage <- 5*(length((mxF[mxF$country_code==cc,3]))-2)
  asdr.object <- new("nMx",
                     rate.years = seq(1950, 2095, 5), #xxamy - added rate.years
                     rates = asdr.1950.2100.by5,
                     mid.age = c(0.5, seq(2.5, (asdr.maxage+2.5), 5)))

  #Life Expectancy over time
  data("e0F", package = "wpp2017", envir = environment()) # wpp2017
  data("e0Fproj", package = "wpp2017", envir = environment()) # wpp2017

  #data("e0F_wpp2017")
  #e0F <- e0F_wpp2017

  #data("e0Fproj_wpp2017")
  #e0Fproj <- e0Fproj_wpp2017


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
