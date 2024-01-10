#' Function to get all the things you need to run setup coverage and demography for country (NOT Montagu)
#'
#' @param country character; country name
#'
#' @import countrycode
#' @import dplyr
#' @importFrom utils read.csv
#' @importFrom readxl read_xlsx
#'
#' @return demography and coverage for each country
#' @export
#'

setupCountry.Dec2021 <- function(country){

  # state, tran not defined for the get births function below
  uncode <- countrycode::countrycode(country, origin="country.name", destination="un")
  iso3code <- countrycode::countrycode(country, origin="country.name", destination="iso3c")

  demog <- getDemography.wpp2019(uncode=uncode)
  pop.total.1950.2100 <- demog$pop.total.1950.2100*1000
  pop.age.byageclasses.1950.2100 <- demog$pop.age.byageclasses.1950.2100*1000
  tfr.1950.2100 <- demog$tfr.1950.2100
  asfr.1950.2100 <- demog$asfr.1950.2100
  e0.1950.2100 <- demog$e0.1950.2100
  fert <- c(0,asfr.1950.2100[,(2015-1950+1)]*1000,0) # fert is defined here
  births.1950.2100 <- demog$births.1950.2100
  cbr.1950.2100 <- demog$cbr.1950.2100
  asdr.1950.2100.by5 <- demog$asdr.1950.2100.by5
  asdr.object <- demog$asdr.object
  repro.age.sex.dist.1950.2100 <- demog$repro.age.sex.dist.1950.2100

  # MCV routine coverage
  filep <- system.file("/extdata/MCV_WUENIC_coverage_estimates_DownloadedDec2021.xlsx",
                       package = "MRTransmissionModel")
  df <- read_xlsx(path = filep, sheet="Sheet1")

  #df <- read_xlsx("./data/MCV_WUENIC_coverage_estimates_DownloadedDec2021.xlsx", sheet="Sheet1")
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
  filep2 <- system.file("/extdata/RCV_WUENIC_coverage_estimates_DownloadedDec2021.xlsx",
                       package = "MRTransmissionModel")
  df <- read_xlsx(path = filep2, sheet="Sheet1")

  #df <- read_xlsx("./data/RCV_WUENIC_coverage_estimates_DownloadedDec2021.xlsx", sheet="Sheet1")
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
                              lower.age.boundary = c(0,15,20,25,30,35,40,45,50,101)) {
    #xxamy - age.class change from lower.age.boundary = c(0,15,20,25,30,35,40,45,49,100)) {
    return(get.births.fun(state=state,tran=tran,fert.curve=fert.curve,lower.age.boundary=lower.age.boundary))

  }

  #Inaccessible from DPT1 numbers over time
  dtp <- MRTransmissionModel::dtp1_estimates

  #dtp <- read.csv("./data/dtp1_estimates.csv")
  dtp.c <- dtp[dtp$ISO_code==iso3code,]
  dtp.c <- dtp.c[order(dtp.c$Year),]
  dtp.1980.to.2100 <- dtp.c$vacc_coverage[1:length(1980:2017)]
  #filling in zeros using zoo package
  dtp.1980.to.2100 <- zoo::na.locf(dtp.1980.to.2100, fromLast=T, na.rm=F) #fills in using value after NA
  dtp.1980.to.2100 <- zoo::na.locf(dtp.1980.to.2100) #fills in using previous value before NA
  dtp.1980.to.2100 <- c(dtp.1980.to.2100, rep(dtp.1980.to.2100[length(dtp.1980.to.2100)], length(2018:2101))) #add one more year 1980to2101

  if (length(dtp.1980.to.2100)==0){ #Palestine and Kosovo dont have DTP1 estimates
    # should this be RCV1.coverage?
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
