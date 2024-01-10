#' Function to pull Age Specific Death Rates and return survivorship over time
#'
#' @param age.classes a set of age classes for which tran is being built
#' @param generation.time generation time
#' @param nMx object of age specific death rates over time and mid-age per rate
#' @param nMx.years years to interpolate age specific death rates
#' @param check logical;
#'
#' @return returns the age profile (rows) of survivorship (in units of the generation time) by time (cols) as a matrix
#' @export
#'
create.surv.prob.over.age.time <- function(age.classes, generation.time, nMx=NULL, nMx.years, check=F){

  mid.age <- nMx@mid.age
  generation.time.year <- generation.time/12   # put generation.time in terms of years rather than months (0.5 in months) now (0.041667 in years)
  no.gens.in.year <- 12*1/generation.time

  rate.years <- nMx@rate.years #seq(1950,2100,5) #xxamy - revised what feeds into rate.years here
  indexes <- findInterval(nMx.years, rate.years)
  rates.years <- data.frame(nMx@rates[,indexes])
  rates.generations <- rates.years[,rep(seq_len(ncol(rates.years)), each=no.gens.in.year)]
  rates.generations <- cbind(rates.generations[,1], rates.generations) #adding additional column to match the length of the experiment time steps

  survs <- matrix(NA, nrow=length(age.classes), ncol=ncol(rates.generations))

  for (t in 1:ncol(rates.generations)){
    yearly.mortality <- rates.generations[,t]

    # Fit a smooth spline to this
    fit <- smooth.spline(mid.age, log(yearly.mortality), df=length(mid.age))

    # Convert the fit to our ages
    survs[,t] <- 1-exp(predict(fit, age.classes/12)$y)
    # gives me only the $y or estimated m.x.n for age.classes/12 (puts age.classes in years rather than months)
    # for each age predict the survival rate based on the smooth.spline
    # exp to get back to m.x.n from it's previous log state and take 1-m.x.n to get survival rate
    #plot(survs)

    if (any(survs[,t]<0)){
      print(paste0("age-specific survival model fit producing NAs, see create.surv.prob.over.age.time() at time:", t))
      if (t!=1) survs[,t] <-  survs[,(t-1)] #hack, use the previous year, return to later
      if (t==1) { #hack if year 1 to use first non-negative
        survs[which(survs[,t]<0),t] <-  survs[min(which(survs[,t]>0)),t]
      }
    }
    # Convert yearly rates to generation time rates
    survs[,t] <- 1+(log(survs[,t])/(1/generation.time.year))

    # Need to adjust last age class to deplete
    survs[,t][length(survs[,t])] <-  0.5^no.gens.in.year
  }

  if (check) {
    plot(mid.age, yearly.mortality, type="b", xlim=c(0,50) ,ylim=range(yearly.mortality[mid.age<60]),pch=19)
    points(age.classes/12,(1-survs[,t]^(1/generation.time.year)), col=2,type="b")
    points(fit$x,exp(fit$y),col=4, type="l")
  }
  return(survs)

}
