#' Aggregate a model state vector from monthly age classes to single-year age groups
#'
#' Sums the fine-grained (1-month) age classes into single-year bins and
#' appends the coarser (annual) age classes for older ages unchanged.
#'
#' @param vec numeric vector. State vector for one epi class and one time point,
#'   with one element per age class in \code{age.classes}.
#' @param age.classes numeric vector. Upper bounds of age classes in months,
#'   as used throughout the model (e.g. \code{c(1:240, seq(252, 1212, 12))}).
#'
#' @return numeric vector of length equal to the number of single-year age groups
#'   spanned by \code{age.classes}, with values summed across the 12 monthly
#'   sub-classes for each year.
#' @export
#'

GetNumber.per.AgeYear <- function(vec, age.classes) {

  upper.age.year <- (age.classes[length(age.classes)]-12)/12  #xxamy - age.class change from upper.age.year <- (age.classes[length(age.classes)]-1)/12 #using only to age 58, because by 59 all mostly dead
  diff <- diff(age.classes)
  top.one.month.age <- age.classes[which(diff==unique(diff)[2])[1]]
  top.one.month.age.inyears <- (top.one.month.age)/12 #xxamy - age.class change from top.one.month.age.inyears <- (top.one.month.age-1)/12

  #pop.per.youngage.year <- rep(0,top.one.month.age.inyears)
  #for (y in 1:top.one.month.age.inyears) {
  #  for (u in 0:11){
  #    pop.per.youngage.year[y] <- pop.per.youngage.year[y]+vec[(12*y)-u] #xxamy - age.class change from pop.per.youngage.year[y] <- pop.per.youngage.year[y]+vec[(12*y+1)-u] #xxamy - leaves out age group 1
  #  }
  #}
  # the following replaced for former to efficiency
  pop.per.youngage.year <- colSums(matrix(vec[seq_len(top.one.month.age)], nrow = 12))

  pop.per.age.year <- c(pop.per.youngage.year, vec[(top.one.month.age+1):length(vec)])
  return(pop.per.age.year)
}
