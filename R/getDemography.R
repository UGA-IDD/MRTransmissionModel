#' Function to get population size, age structure, fertility and death rates from 1950 to 2105
#'
#' @param uncode UN country code
#'
#' @importFrom utils read.csv
#'
#' @return demography for each country
#' @export
#'

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
