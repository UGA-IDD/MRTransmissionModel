#' Get School Vaccination time age-specific
#'
#' @param time.step xxx
#' @param age.classes xxx
#' @param time.specific.cov xxx
#' @param age.min xxx
#' @param age.max xxx
#' @param obj.vcdf xxx
#' @param obj.prob.vsucc xxx
#'
#' @return list for routine time
#' @export
#'

get.schoolvacc.time.age.specific <- function(time.step=0.5, age.classes = c(1:240, seq(252,1212,12)),
                                             time.specific.cov = rep(0.6, 50), age.min=rep(12, 50), age.max=rep(23,50),
                                             obj.vcdf=get.vcdf.uniform(12, 23),
                                             obj.prob.vsucc = pvacsuccess(1:(14*12), new("vsucc.constant", success.rate=0.8))){

  time.specific.routine <- prop.fail.MR1.byage <- matrix(0, length(time.specific.cov), length(age.classes))
  prop.fail.MR1 <- one.minus.ve1 <- rep(0, length(time.specific.cov))

  #get probability of successful school vaccination
  for (j in 1:length(time.specific.cov)) {
    if (time.specific.cov[j]!=0){

      cdf.vaccination <- obj.vcdf@cdf[obj.vcdf@ages %in% c(age.min[j]:age.max[j])]
      prob.success <- obj.prob.vsucc@prob.vsucc[obj.prob.vsucc@ages %in% c(age.min[j]:age.max[j])]
      if (length(cdf.vaccination)!=length(prob.success)) stop("problem matching age in vcdf and vsucc, school vaccination")
      cdf.scaled.vaccination <- time.specific.cov[j]*(cdf.vaccination) #xxamy - different from get.routine.time.age.specific() b/c time.specific.cov[j]*(cdf.vaccination/cdf.vaccination[length(cdf.vaccination)])
      cdf.scaled.vaccination <- c(0,cdf.scaled.vaccination)
      pdf <- diff(cdf.scaled.vaccination)
      h <- pdf/(1-cdf.scaled.vaccination[2:length(cdf.scaled.vaccination)-1]) #hazard
      low.age <- (age.min[j]-1):(age.max[j]-1)
      high.age <- age.min[j]:age.max[j]
      age.sz <- high.age-low.age
      ts.per.class <- age.sz/time.step
      final.pdf <- 1-(1-h)^(1/ts.per.class)
      pdf.scaled.vaccination <- pmax(final.pdf,0)
      time.specific.routine[j,age.classes %in% c(age.min[j]:age.max[j])] <-  pdf.scaled.vaccination*prob.success #pdf*prob.success.MR1 = rep(0.04, 12); sum(rep(0.04, 12)) = 0.48 = ve1*mcv1, check
      #time.specific.routine[j,age.classes %in% c(age.min.MR1[j]:age.max.MR1[j])] <-  pdf*prob.success.MR1
      one.minus.ve1[j] <- sum(pdf/sum(pdf)*(1-prob.success)) #using pdf/sum(pdf) to get weighted average 1-VE1 = 0.2, check
      prop.fail.MR1[j] <- sum(pdf*(1-prob.success)) # = 0.12 = (1-ve1)*mcv1, proportion of total population with primary failure from MR1, check
      prop.fail.MR1.byage[j,age.classes %in% c(age.min[j]:age.max[j])] <- pdf*(1-prob.success)

    }
  }

  return(list(
    age.time.specific.schoolvacc=time.specific.routine,
    one.minus.ve.schoolvacc=one.minus.ve1,
    prop.fail.schoolvacc=prop.fail.MR1,
    prop.fail.schoolvacc.byage = prop.fail.MR1.byage))

}
