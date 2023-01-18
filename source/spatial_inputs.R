## Functions to create data inputs based on subnational vaccination data and subnational demography data
## Authors: Amy Winter


### Function to get all the things you need to run
### setup coverage and demography for country (NOT Montagu)
###
##'  @param - country - character - country name 
###
### Returns demography and coverage for each country
space.setupCountry <- function(country){ 
  
  uncode <- countrycode::countrycode(country, origin="country.name", destination="un")
  iso3code <- countrycode::countrycode(country, origin="country.name", destination="iso3c")
  
  demog <- space.getDemography(uncode=uncode)
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
  
  # MCV1
  library(tidyverse)
  df <- read.csv("data/ihme_mcv1_coverage_ad2.csv")
  MCV1.coverage <- df %>%
    filter(ISO3==iso3code) %>%
    select(Admin2, Year, "MCV1...Mean") %>%
    pivot_wider(names_from=Year, values_from="MCV1...Mean") %>%
    rename(district = Admin2) %>%
    arrange(district)
  
  ## MCV2 
  load("/Users/winter/Dropbox/MeaslesPostCOVIDZambia/analyses/DHSDistrictEstimates/output/DHS2018_115/params_MCV2.RData") #loading prop.access
  mcv2.df.tmp <- prop.access
  dist.names <- names(mcv2.df.tmp)
  mcv2 <- as.vector(mcv2.df.tmp)
  dist.names[is.na(match(dist.names, mcv1.df$district))]
  mcv1.df$district[is.na(match(mcv1.df$district, dist.names))]
  dist.names <- recode(dist.names, 
                            `Milenge`="Milengi",
                            `Shiwang'Andu`="Shiwan'gandu",
                            `Lavushi Manda`="Lavushimanda",
                            `Senga`="Senga Hill")
  #add the missing districts and give them the mean mcv2
  mcv2 <- c(mcv2, rep(mean(mcv2), length(mcv1.df$district[is.na(match(mcv1.df$district, dist.names))])))
  dist.names <- c(dist.names, mcv1.df$district[is.na(match(mcv1.df$district, dist.names))])
  mcv2.df <- tibble(district = dist.names, '2018' = mcv2) %>% arrange(district)
  MCV2.coverage <- MCV1.coverage %>% naniar::replace_with_na_all(condition = ~ .x < 1)
  MCV2.coverage[,which(colnames(MCV1.coverage)=="2018")] <- mcv2.df$`2018`
  
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
              MCV1.coverage.2000to2021=MCV1.coverage,
              MCV2.coverage.2000to2021=MCV2.coverage))
}