 #' Function to pull in Age Specific Death Rates and return survivorship for running out transients
#'
#' @param age.classes a set of age classes for which tran is being built
#' @param generation.time generation time
#' @param check xxx
#' @param nMx object of age specific death rates over time and mid-age per rate
#' @param year year to pull age specific death rates that will be used to run out the transients (should coincide with DFE year)
#'
#' @importFrom stats smooth.spline predict
#' @importFrom graphics points
#'
#' @return returns the age profile of survivorship in units of the generation time
#' @export
#'

create.surv.prob.over.age <- function(age.classes, generation.time, check=F, nMx=NULL, year=1990){

  rate.years <- seq(1950,2100,5)
  index <- min(which(findInterval(rate.years, year)==1))
  rates.overtime <- nMx@rates
  mid.age <- nMx@mid.age
  generation.time.year <- generation.time/12   # put generation.time in terms of years rather than months (0.5 in months) now (0.041667 in years)

  yearly.mortality <- rates.overtime[,index]

  # Fit a smooth spline to this
  fit <- smooth.spline(mid.age, log(yearly.mortality), df=length(mid.age))

  # Convert the fit to our ages
  survs <- 1-exp(predict(fit, age.classes/12)$y)
  # gives me only the $y or estiamted m.x.n for age.classes/12 (puts age.classes in years rather than months)
  # for each age predict the survival rate based on the smooth.spline
  # exp to get back to m.x.n from it's previous log state and take 1-m.x.n to get survival rate
  #plot(survs)

  # Convert yearly rates to generation time rates
  survs <- 1+(log(survs)/(1/generation.time.year))

  # Need to adjust last age class to deplete
  survs[length(survs)] <-  0.5^(12*1/generation.time)

  if (check) {
    plot(mid.age, yearly.mortality, type="b", xlim=c(0,100) ,ylim=range(yearly.mortality[mid.age<100]),pch=19)
    points(age.classes/12,(1-survs^(1/generation.time.year)), col=2,type="b")
    points(fit$x,exp(fit$y),col=4, type="l")
  }

  return(survs)
}
