
#' Function nested inside getDemography.wpp2022() to only pull pop and age.structure
#'
#' @param uncode UN country code
#' @param age.classes age classes in years
#' @param if.missing.use.region logical; c(MHL, XK, TUV)
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

getDemography.wpp2022.mini <- function(uncode=NA, age.classes=c(1:101),
                                       if.missing.use.region=F){

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
    stop(paste("No UNPD data for uncode ", uncode, "; NAs produced", sep=""))
    #print(paste("No UNPD data for uncode ", uncode, "; NAs produced", sep=""))
    #pop.total.1950.2100=NA
    #pop.age.byageclasses.1950.2100=NA
    #tfr.1950.2100=NA
    #births.1950.2100=NA
    #cbr.1950.2100=NA
    #e0.1950.2100=NA
    #asfr.1950.2100=NA
    #repro.age.sex.dist.1950.2100=NA
    #asdr.1950.2100.by5=NA
    #asdr.object=NA
  }

  #Population total in years 1950 to 2100
  data("pop1dt", package = "wpp2022", envir = environment())
  data("popproj1dt", package = "wpp2022", envir = environment())

  pop.total.1950.2100 <- bind_rows(pop1dt %>% dplyr::filter(country_code==cc) %>% dplyr::select(year, pop),
                                   popproj1dt %>% dplyr::filter(country_code==cc) %>% dplyr::select(year, pop)) %>%
    dplyr::filter(year %in% time) %>%
    dplyr::pull(pop)

  #Population Age Structure by Age over time
  data("popAge1dt", package = "wpp2022", envir = environment())
  data("popprojAge1dt", package = "wpp2022", envir = environment())

  pop.age.byageclasses.1950.2100 <- bind_rows(popAge1dt %>% dplyr::filter(country_code==cc) %>% dplyr::select(year, age, pop) %>% dplyr::mutate(year=as.integer(year)),
                                              popprojAge1dt %>% dplyr::filter(country_code==cc) %>% dplyr::select(year, age, pop) %>% dplyr::mutate(year=as.integer(year))) %>%
    dplyr::mutate(age = age+1) %>%
    dplyr::filter(year %in% time,
                  age %in% age.classes) %>%
    tidyr::pivot_wider(names_from = year, values_from = pop) %>%
    tibble::column_to_rownames(var="age")

  return(list(cc=cc,
              time=time,
              time.by5=time.by5,
              pop.total.1950.2100=pop.total.1950.2100,
              pop.age.byageclasses.1950.2100=pop.age.byageclasses.1950.2100))

}
