#' Function to get all the things you need to run setup coverage and demography for country (Montagu for touchstone 202110gavi_v3)
#'
#' @param country country iso3code
#'
#' @return demography and coverage for each country
#' @export
#'

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

  #Coverage estimates from Montagu, and inaccessible population
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
