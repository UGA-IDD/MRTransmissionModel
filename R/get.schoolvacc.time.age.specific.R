#' Get School Vaccination time age-specific
#'
#' @param time.step numeric - in months
#' @param age.classes vector - in months
#' @param time.specific.cov vector
#' @param age.min vector - in months
#' @param age.max vector - in months
#' @param obj.vcdf object
#' @param obj.prob.vsucc object
#'
#' @return list for school enrollment vaccination by age and time
#' @export
#'

get.schoolvacc.time.age.specific <- function(time.step=0.5, age.classes = c(1:240, seq(252,1212,12)),
                                             time.specific.cov = rep(0.6, 50), age.min=rep(12, 50), age.max=rep(23,50),
                                             obj.vcdf=get.vcdf.uniform(36, 12*15),
                                             obj.prob.vsucc = pvacsuccess(1:(14*12), new("vsucc.constant", success.rate=0.8))){

  time.specific.routine <- prop.fail.byage <- matrix(0, length(time.specific.cov), length(age.classes))
  prop.fail <- one.minus.ve <- rep(0, length(time.specific.cov))

  #get probability of successful school vaccination
  for (j in 1:length(time.specific.cov)) {
    if (time.specific.cov[j]!=0){

      #old
      #cdf.vaccination <- obj.vcdf@cdf[obj.vcdf@ages %in% c(age.min[j]:age.max[j])]
      #prob.success <- obj.prob.vsucc@prob.vsucc[obj.prob.vsucc@ages %in% c(age.min[j]:age.max[j])]
      #if (length(cdf.vaccination)!=length(prob.success)) stop("problem matching age in vcdf and vsucc, school vaccination")
      #cdf.scaled.vaccination <- time.specific.cov[j]*(cdf.vaccination) #xxamy - different from get.routine.time.age.specific() b/c time.specific.cov[j]*(cdf.vaccination/cdf.vaccination[length(cdf.vaccination)])
      #cdf.scaled.vaccination <- c(0,cdf.scaled.vaccination)
      #pdf <- diff(cdf.scaled.vaccination)
      #h <- pdf/(1-cdf.scaled.vaccination[2:length(cdf.scaled.vaccination)-1]) #hazard
      #low.age <- (age.min[j]-1):(age.max[j]-1)
      #high.age <- age.min[j]:age.max[j]
      #age.sz <- high.age-low.age
      #ts.per.class <- age.sz/time.step
      #final.pdf <- 1-(1-h)^(1/ts.per.class)
      #pdf.scaled.vaccination <- pmax(final.pdf,0)
      #time.specific.routine[j,age.classes %in% c(age.min[j]:age.max[j])] <-  pdf.scaled.vaccination*prob.success #pdf*prob.success.MR1 = rep(0.04, 12); sum(rep(0.04, 12)) = 0.48 = ve1*mcv1, check
      #one.minus.ve1[j] <- sum(pdf/sum(pdf)*(1-prob.success)) #using pdf/sum(pdf) to get weighted average 1-VE1 = 0.2, check
      #prop.fail.MR1[j] <- sum(pdf*(1-prob.success)) # = 0.12 = (1-ve1)*mcv1, proportion of total population with primary failure from MR1, check
      #prop.fail.MR1.byage[j,age.classes %in% c(age.min[j]:age.max[j])] <- pdf*(1-prob.success)

      #New - School vaccination should be more akin to SIA vaccination, rather than routine (which b/c happening at each discrete time interval, "everyone" is exposed to all probabilities in the pdf)
      cdf.vaccination <- obj.vcdf@cdf[obj.vcdf@ages %in% c(age.min[j]:age.max[j])]
      prob.success <- obj.prob.vsucc@prob.vsucc[obj.prob.vsucc@ages %in% c(age.min[j]:age.max[j])]
      if (length(cdf.vaccination)!=length(prob.success)) stop("problem matching age in vcdf and vsucc, school vaccination")
      cdf.scaled.vaccination <- time.specific.cov[j]*(cdf.vaccination) #xxamy - different from get.routine.time.age.specific() b/c time.specific.cov[j]*(cdf.vaccination/cdf.vaccination[length(cdf.vaccination)])
      cdf.scaled.vaccination <- c(0,cdf.scaled.vaccination)
      pdf <- diff(cdf.scaled.vaccination)
      pdf <- c(pdf[1], rep(rollsum(pdf[2:length(pdf)], 12)[seq(1, length(rollsum(pdf[2:length(pdf)], 12)), 12)], each=12)) #each prob is used 12 times
      #pdf <- pdf*12 #xxamy hack - generally works excpet when rate is super high then get prob>1
      ##xxamy, this decision was made b/c individuals are only eligible once a year, and the pdf intervals are in months, therefore multiply by 12
      ##xxamy, thought experiment, if only 7yo were enrolling in school at 95%, then we would create a vector pdf <- rep(0.95, 12), where sum(pdf)=12*0.95
      ##xxamy, thought experiment continued, cdf above plateaus at 0.95 across more ages, therefore sum(pdf)=0.95, but we want sum to be 12*0.95, therefore multiply pdf by 12
      ##xxamy, thought experiment continued, if 7yo and 8yo were enrolling in school at 95%, then we would create a vector pdf <- rep(0.95/2, 24), b/c a birth cohort will have two opportunities to enroll once at 7yo and then again one year later at 8yo, therfore the proportion enrolled is 0.475 in each year
      time.specific.routine[j,age.classes %in% c(age.min[j]:age.max[j])] <-  pdf*prob.success
      one.minus.ve[j] <- sum(pdf/sum(pdf)*(1-prob.success)) #using pdf/sum(pdf) to get weighted average 1-VE1 = 0.2, check
      prop.fail[j] <- sum(pdf*(1-prob.success)) # = 0.12 = (1-ve1)*mcv1, proportion of total population with primary failure from MR1, check
      prop.fail.byage[j,age.classes %in% c(age.min[j]:age.max[j])] <- pdf*(1-prob.success)

    }
  }

  return(list(
    age.time.specific.schoolvacc=time.specific.routine,
    one.minus.ve.schoolvacc=one.minus.ve,
    prop.fail.schoolvacc=prop.fail,
    prop.fail.schoolvacc.byage = prop.fail.byage))

}
