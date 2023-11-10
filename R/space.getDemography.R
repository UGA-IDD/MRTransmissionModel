#' Function to get spatial demography - this function is temporary and will need to be totally revamped
#'
#' @param uncode UN country code
#'
#' @return list, demographic data
#' @export
#'
space.getDemography <- function(uncode){

  # library(logspline)
  # setOldClass("oldlogspline")
  # library(KernSmooth)
  # library(scales) #for alpha function
  # library(zoo) #to fill in NAs from the dpt1 estimates
  # library(survival) #for the PDF of age of vaccination
  # library(countrycode) #to transfer b/w uncode and iso3codes and country names
  # library(readxl)
  # library(dplyr)
  #
  # dyn.load("./source/MRModel-funcs.so") #run "R CMD SHLIB source/MRModel-funcs.c" in terminal to compile
  # source("./source/build.R")
  # source("./source/base.R")
  # source("./source/user_interface.R")
  # source("./source/who_un_inputs.R")
  # source("./source/new_functions.R")
  #
  # library(wpp2019) #UN and WHO data model inputs require wpp2019
  setup <- setupCountry.Dec2021(country="Zambia")
  year <- 1980
  t.max <- length(year:2100)
  generation.time <- 0.5
  age.classes <- c(1:240, seq(252,1212,12))

  pop.total.1950.2100 <- rbind(setup$pop.total.1950.2100, setup$pop.total.1950.2100*0.5)
  pop.age.byageclasses.1950.2100 <- rbind(setup$pop.age.byageclasses.1950.2100, setup$pop.age.byageclasses.1950.2100*0.5)
  tfr.1950.2100 <- rbind(setup$tfr.1950.2100, setup$tfr.1950.2100)
  e0.1950.2100 <- rbind(setup$e0.1950.2100, setup$e0.1950.2100)
  asfr.1950.2100 <- rbind(setup$asfr.1950.2100, setup$asfr.1950.2100)
  repro.age.sex.dist.1950.2100 <- rbind(setup$repro.age.sex.dist.1950.2100, setup$repro.age.sex.dist.1950.2100)
  births.1950.2100 <- rbind(setup$births.1950.2100, setup$births.1950.2100)
  cbr.1950.2100 <- rbind(setup$cbr.1950.2100, setup$cbr.1950.2100)
  yr.agespecificbirths.per.1000 <- rbind(setup$yr.agespecificbirths.per.1000, setup$yr.agespecificbirths.per.1000)
  asdr.1950.2100.by5 <- rbind(setup$asdr.1950.2100.by5, setup$asdr.1950.2100.by5)
  asdr.object <- new("space.nMx",
                     rate.years = setup$asdr.object@rate.years,
                     rates = rbind(setup$asdr.object@rates, setup$asdr.object@rates),
                     mid.age = setup$asdr.object@mid.age,
                     n.subpops = 2)

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


