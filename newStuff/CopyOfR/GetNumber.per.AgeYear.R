#' Function to get number of individuals per age in years
#'
#' @param vec vector of length age.classes
#' @param age.classes age.classes
#'
#' @return number of individuals per age in years
#' @export
#'

GetNumber.per.AgeYear <- function(vec, age.classes) {

  upper.age.year <- (age.classes[length(age.classes)]-12)/12  #xxamy - age.class change from upper.age.year <- (age.classes[length(age.classes)]-1)/12 #using only to age 58, because by 59 all mostly dead
  diff <- diff(age.classes)
  top.one.month.age <- age.classes[which(diff==unique(diff)[2])[1]]
  top.one.month.age.inyears <- (top.one.month.age)/12 #xxamy - age.class change from top.one.month.age.inyears <- (top.one.month.age-1)/12

  pop.per.youngage.year <- rep(0,top.one.month.age.inyears)
  for (y in 1:top.one.month.age.inyears) {
    for (u in 0:11){
      pop.per.youngage.year[y] <- pop.per.youngage.year[y]+vec[(12*y)-u] #xxamy - age.class change from pop.per.youngage.year[y] <- pop.per.youngage.year[y]+vec[(12*y+1)-u] #xxamy - leaves out age group 1
    }
  }

  pop.per.age.year <- c(pop.per.youngage.year, vec[(top.one.month.age+1):length(vec)])
  return(pop.per.age.year)
}
