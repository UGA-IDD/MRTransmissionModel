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
#' @import wpp2022
#' @import tidyr
#' @import dplyr
#' @import tibble
#'
#' @return demography for each country
#' @export
#'

getDemography.wpp2022 <- function(uncode=NA, age.classes=c(1:101),
                                  if.missing.use.region=F){
  # nice to have accessible to user
  cc <- uncode

  time.by5 <- seq(1950,2100,5)
  time <- seq(1950,2100,1)

  #UNPD data
  data("pop", package = "wpp2022", envir = environment()) # wpp2022

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
  data("pop1dt", package = "wpp2022", envir = environment())
  data("popproj1dt", package = "wpp2022", envir = environment())

  pop.total.1950.2100 <- bind_rows(pop1dt %>% dplyr::filter(country_code==uncode) %>% dplyr::select(year, pop),
                       popproj1dt %>% dplyr::filter(country_code==uncode) %>% dplyr::select(year, pop)) %>%
    dplyr::filter(year %in% time) %>%
    dplyr::pull(pop)

  #Population Age Structure by Age over time
  data("popAge1dt", package = "wpp2022", envir = environment())
  data("popprojAge1dt", package = "wpp2022", envir = environment())

  pop.age.byageclasses.1950.2100 <- bind_rows(popAge1dt %>% dplyr::filter(country_code==uncode) %>% dplyr::select(year, age, pop) %>% dplyr::mutate(year=as.integer(year)),
                                              popprojAge1dt %>% dplyr::filter(country_code==uncode) %>% dplyr::select(year, age, pop) %>% dplyr::mutate(year=as.integer(year))) %>%
    dplyr::mutate(age = age+1) %>%
    dplyr::filter(year %in% time,
                  age %in% age.classes) %>%
    tidyr::pivot_wider(names_from = year, values_from = pop) %>%
    tibble::column_to_rownames(var="age")

  #Sex Distribution grouped by 5 year reproductive ages over time
  repro.age.sex.dist.1950.2100 <- bind_rows(popAge1dt %>% dplyr::filter(country_code==uncode) %>% dplyr::select(year, age, pop, popF) %>% dplyr::mutate(year=as.integer(year)),
                                              popprojAge1dt %>% dplyr::filter(country_code==uncode) %>% dplyr::select(year, age, pop, popF) %>% dplyr::mutate(year=as.integer(year))) %>%
    dplyr::mutate(age = age) %>%
    dplyr::filter(year %in% time,
                  age %in% 15:49) %>%
    dplyr::mutate(prop = popF/pop) %>%
    dplyr::select(age, year, prop) %>%
    tidyr::pivot_wider(names_from = year, values_from = prop) %>%
    tibble::column_to_rownames(var="age")

  #TFR over time
  data("tfr1dt", package = "wpp2022", envir = environment())
  data("tfrproj1dt", package = "wpp2022", envir = environment())

  tfr.1950.2100 <- bind_rows(tfr1dt %>% dplyr::filter(country_code==uncode) %>% dplyr::select(year, tfr),
                             tfrproj1dt %>% dplyr::filter(country_code==uncode) %>% dplyr::select(year, tfr)) %>%
    dplyr::filter(year %in% time) %>%
    dplyr::pull(tfr)
  names(tfr.1950.2100) <- time

  #ASFR overtime
  data("percentASFR1dt", package = "wpp2022", envir = environment())
  asfr.1950.2100 <- percentASFR1dt %>%
    dplyr::filter(country_code==uncode) %>%
    dplyr::select(year, age, pasfr) %>%
    dplyr::filter(age %in% 15:49) %>%
    tidyr::pivot_wider(names_from = year, values_from = pasfr) %>%
    tibble::column_to_rownames(var="age")

  #Births over time
  data("misc1dt", package = "wpp2022", envir = environment())
  data("miscproj1dt", package = "wpp2022", envir = environment())
  births.1950.2100 <- bind_rows(misc1dt %>% dplyr::filter(country_code==uncode) %>% dplyr::select(year, births),
                                miscproj1dt %>% dplyr::filter(country_code==uncode) %>% dplyr::select(year, births)) %>%
    dplyr::filter(year %in% time) %>%
    dplyr::pull(births)
  names(births.1950.2100) <- time

  #Crude Birth Rates over time
  cbr.1950.2100 <- bind_rows(misc1dt %>% dplyr::filter(country_code==uncode) %>% dplyr::select(year, cbr),
                             miscproj1dt %>% dplyr::filter(country_code==uncode) %>% dplyr::select(year, cbr)) %>%
    dplyr::filter(year %in% time) %>%
    dplyr::pull(cbr)
  names(cbr.1950.2100) <- time

  #Age Specific Death Rate over time
  data("mx1dt", package = "wpp2022", envir = environment())
  asdr.1950.2100.long <- mx1dt %>%
    dplyr::filter(country_code==uncode) %>%
    dplyr::select(year, age, mxB) %>%
    dplyr::mutate(age = age+1) %>%
    dplyr::filter(year %in% time,
                  age %in% age.classes)
  asdr.1950.2100 <- asdr.1950.2100.long %>%
    tidyr::pivot_wider(names_from = year, values_from = mxB) %>%
    tibble::column_to_rownames(var="age")
  asdr.object <- new("nMx",
                     rate.years = unique(asdr.1950.2100.long$year),
                     rates = asdr.1950.2100,
                     mid.age = unique(asdr.1950.2100.long$age)-0.5)


  #Life Expectancy over time
  data("e01dt", package = "wpp2022", envir = environment())
  data("e0proj1dt", package = "wpp2022", envir = environment())
  e0.1950.2100 <- bind_rows(e01dt %>% dplyr::filter(country_code==uncode) %>% dplyr::select(year, e0B),
                            e0proj1dt %>% dplyr::filter(country_code==uncode) %>% dplyr::select(year, e0B)) %>%
    dplyr::filter(year %in% time) %>%
    dplyr::pull(e0B)
  names(e0.1950.2100) <- time

  return(list(pop.total.1950.2100=pop.total.1950.2100,
              pop.age.byageclasses.1950.2100=pop.age.byageclasses.1950.2100,
              tfr.1950.2100=tfr.1950.2100,
              e0.1950.2100=e0.1950.2100,
              asfr.1950.2100=asfr.1950.2100, #different from wpp2019 b/c 35x151 instead of 7x151
              repro.age.sex.dist.1950.2100=repro.age.sex.dist.1950.2100, #different from wpp2019 b/c 35x151 instead of 7x151
              births.1950.2100=births.1950.2100,
              cbr.1950.2100=cbr.1950.2100,
              asdr.1950.2100=asdr.1950.2100, #different from wpp2019 b/c 101x151 instead of 22x31
              asdr.object=asdr.object)) #different from wpp2019 b/c 101x151 instead of 22x31
}
