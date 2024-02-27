#' Function to get spatial demography
#'
#' @param age.classes numeric - in months
#' @param iso3code character
#'
#' @import dplyr
#' @import schoolentrydataZambia
#'
#' @return list, demographic data
#'
space.getDemography <- function(age.classes, iso3code){

  library(schoolentrydataZambia)
  zamb_data <- schoolentrydataZambia::getDemography(iso3code)
  subpop.names <- zamb_data$pop.total.1950.2100 %>% arrange(district) %>% pull(district)

  # arranging all by district name so that will be consistent with subpop.names above
  pop.total.1950.2100 <- zamb_data$pop.total.1950.2100 %>% arrange(district) %>% select(-c(1:3)) %>% as.data.frame()
  cbr.1950.2100 <- zamb_data$cbr.1950.2100 %>% arrange(district) %>% select(-c(1:3)) %>% as.data.frame()
  asdr.1950.2100.by5 <- zamb_data$asdr.1950.2100 %>% arrange(district) %>% select(-c(1:4)) %>% as.data.frame()
  asdr.rate.years <- zamb_data$asdr.rate.years
  asdr.rate.mid.age <- zamb_data$asdr.rate.mid.age

  # getting pop by by age.classes - and make sure consistent with pop total
  pop.age.1950.2100 <- zamb_data$pop.age.byageclasses.1950.2100 %>%
    pivot_longer(cols=-c(ISO3, district, province, age), names_to="year", values_to="pop")
  years <- unique(pop.age.1950.2100$year)
  ages <- unique(pop.age.1950.2100$age)
  #setup empty matrix
  pop.age.byageclasses.1950.2100 <- matrix(NA, nrow=(length(age.classes)*length(subpop.names)), ncol=length(years))
  #fill matrix by district and year
  for (d in 1:length(subpop.names)){
    for (y in 1:length(years)){
      tmp.pop.ages <- pop.age.1950.2100 %>% filter(district==subpop.names[d]) %>% filter(year==years[y]) %>% pull(pop)
      f <- smooth.spline(ages, tmp.pop.ages)
      tmp.pop.age.classes <- predict(f,age.classes/12)$y
      #reduce to 1 if negative
      tmp.pop.age.classes[tmp.pop.age.classes<0] <- 1
      #adjust for varying bin width!
      tmp.pop.age.classes <- tmp.pop.age.classes*diff(c(0,age.classes))
      #rescale the age distribution to the estimated population each year
      true.pop <- zamb_data$pop.total.1950.2100  %>%
        pivot_longer(cols=-c(ISO3, district, province), names_to="year", values_to="pop") %>%
        filter(district==subpop.names[d]) %>% filter(year==years[y]) %>% pull(pop)
      tmp.pop.age.classes <- true.pop*(tmp.pop.age.classes/sum(tmp.pop.age.classes))
      #fill in the matrix
      pop.age.byageclasses.1950.2100[(d*length(age.classes)-(length(age.classes)-1)):(d*length(age.classes)),y] <- tmp.pop.age.classes
    }
  }

  return(list(pop.total.1950.2100=pop.total.1950.2100,
              pop.age.byageclasses.1950.2100=pop.age.byageclasses.1950.2100,
              cbr.1950.2100=cbr.1950.2100,
              asdr.1950.2100.by5=asdr.1950.2100.by5,
              asdr.rate.years=asdr.rate.years,
              asdr.rate.mid.age=asdr.rate.mid.age,
              subpop.names=subpop.names))
}


