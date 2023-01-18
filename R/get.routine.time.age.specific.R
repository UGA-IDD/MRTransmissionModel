#' Get Routine time age-specific
#'
#' @param time.step xxx
#' @param age.classes xxx
#' @param time.specific.MR1cov xxx
#' @param age.min.MR1 xxx
#' @param age.max.MR1 xxx
#' @param time.specific.MR2cov xxx
#' @param age.min.MR2 xxx
#' @param age.max.MR2 xxx
#' @param obj.vcdf.MR1 xxx
#' @param obj.vcdf.MR2 xxx
#' @param obj.prob.vsucc xxx
#' @param MR1MR2correlation xxx
#'
#' @return list for routine time
#' @export
#'

get.routine.time.age.specific <- function(time.step=0.5, age.classes = c(1:240, seq(252,1212,12)),
                                          time.specific.MR1cov = rep(0.6, 50), age.min.MR1=rep(12, 50), age.max.MR1=rep(23,50),
                                          time.specific.MR2cov=rep(0.5, 50), age.min.MR2=rep(24, 50), age.max.MR2=rep(35, 50),
                                          obj.vcdf.MR1=get.vcdf.uniform(12, 23), obj.vcdf.MR2=get.vcdf.uniform(24, 35),
                                          obj.prob.vsucc = pvacsuccess(1:(14*12), new("vsucc.constant", success.rate=0.8)),
                                          MR1MR2correlation=F){

  time.specific.routine <- prop.fail.MR1.byage <- prop.fail.MR2.byage <- matrix(0, length(time.specific.MR1cov), length(age.classes))
  prop.fail.MR1 <- one.minus.ve1 <-  prop.fail.MR2 <- one.minus.ve2 <- rep(0, length(time.specific.MR1cov))

  #get probability of successful MR1
  for (j in 1:length(time.specific.MR1cov)) {
    if (time.specific.MR1cov[j]!=0){

      cdf.vaccination.MR1 <- obj.vcdf.MR1@cdf[obj.vcdf.MR1@ages %in% c(age.min.MR1[j]:age.max.MR1[j])]
      prob.success.MR1 <- obj.prob.vsucc@prob.vsucc[obj.prob.vsucc@ages %in% c(age.min.MR1[j]:age.max.MR1[j])]
      if (length(cdf.vaccination.MR1)!=length(prob.success.MR1)) stop("problem matching age in vcdf and vsucc, MR1")
      cdf.scaled.vaccination <- time.specific.MR1cov[j]*(cdf.vaccination.MR1/cdf.vaccination.MR1[length(cdf.vaccination.MR1)])
      cdf.scaled.vaccination <- c(0,cdf.scaled.vaccination)
      pdf <- diff(cdf.scaled.vaccination)
      h <- pdf/(1-cdf.scaled.vaccination[2:length(cdf.scaled.vaccination)-1]) #hazard
      low.age <- (age.min.MR1[j]-1):(age.max.MR1[j]-1)
      high.age <- age.min.MR1[j]:age.max.MR1[j]
      age.sz <- high.age-low.age
      ts.per.class <- age.sz/time.step
      final.pdf <- 1-(1-h)^(1/ts.per.class)
      pdf.scaled.vaccination <- pmax(final.pdf,0)
      time.specific.routine[j,age.classes %in% c(age.min.MR1[j]:age.max.MR1[j])] <-  pdf.scaled.vaccination*prob.success.MR1 #pdf*prob.success.MR1 = rep(0.04, 12); sum(rep(0.04, 12)) = 0.48 = ve1*mcv1, check
      #time.specific.routine[j,age.classes %in% c(age.min.MR1[j]:age.max.MR1[j])] <-  pdf*prob.success.MR1
      one.minus.ve1[j] <- sum(pdf/sum(pdf)*(1-prob.success.MR1)) #using pdf/sum(pdf) to get weighted average 1-VE1 = 0.2, check
      prop.fail.MR1[j] <- sum(pdf*(1-prob.success.MR1)) # = 0.12 = (1-ve1)*mcv1, proportion of total population with primary failure from MR1, check
      prop.fail.MR1.byage[j,age.classes %in% c(age.min.MR1[j]:age.max.MR1[j])] <- pdf*(1-prob.success.MR1)

    }
  }

  one.minus.ve1.lagged <- c(one.minus.ve1[1], one.minus.ve1[-length(one.minus.ve1)]) #lag one year so that it applies to the correct birth cohort
  prop.MR1failANDnoMR1.laggged <- (1-c(time.specific.MR1cov[1], time.specific.MR1cov[-length(time.specific.MR1cov)])) + c(prop.fail.MR1[1], prop.fail.MR1[-length(prop.fail.MR1)]) #lag one year so that it applies to the correct birth cohort = 0.52, check


  #get probability of successful MR2
  for (j in 1:length(time.specific.MR2cov)) {
    if (time.specific.MR2cov[j]!=0){
      cdf.vaccination.MR2 <- obj.vcdf.MR2@cdf[obj.vcdf.MR2@ages %in% c(age.min.MR2[j]:age.max.MR2[j])]
      prob.success.MR2 <- obj.prob.vsucc@prob.vsucc[obj.prob.vsucc@ages %in% c(age.min.MR2[j]:age.max.MR2[j])]
      if (length(cdf.vaccination.MR2)!=length(prob.success.MR2)) stop("prablem matching age in vcdf and vsucc, MR2")
      age.range <- length(age.min.MR2[j]:age.max.MR2[j])
      #if (MR1MR2correlation) cdf.scaled.vaccination <- one.minus.ve1.lagged[j]*time.specific.MR2cov[j]*(cdf.vaccination.MR2/cdf.vaccination.MR2[length(cdf.vaccination.MR2)])/prop.MR1failANDnoMR1.laggged[j] #last index is 0.19230769 = 0.1/0.52 = (mcv2*(1-ve1)) / (1-ve1*mcv1), proportion of susceptibles vaccinated with MCV2, assuming correlation, gtg
      if (MR1MR2correlation) cdf.scaled.vaccination <- one.minus.ve1.lagged[j]*time.specific.MR2cov[j]*(cdf.vaccination.MR2/cdf.vaccination.MR2[length(cdf.vaccination.MR2)]) #last index is 0.1 = (mcv2*(1-ve1)), proportion of total population vaccinated with MCV2, assuming correlation, gtg
      if (!MR1MR2correlation) cdf.scaled.vaccination <- time.specific.MR2cov[j]*(cdf.vaccination.MR2/cdf.vaccination.MR2[length(cdf.vaccination.MR2)]) #last index is 0.5 = mcv2, proportion of total population (or susceptibles) vaccinated with MCV2, no correlation, gtg
      cdf.scaled.vaccination <- c(0,cdf.scaled.vaccination)
      pdf <- diff(cdf.scaled.vaccination)
      h <- pdf/(1-cdf.scaled.vaccination[2:length(cdf.scaled.vaccination)-1]) #hazard
      low.age <- (age.min.MR2[j]-1):(age.max.MR2[j]-1)
      high.age <- age.min.MR2[j]:age.max.MR2[j]
      age.sz <- high.age-low.age
      ts.per.class <- age.sz/time.step
      final.pdf <- 1-(1-h)^(1/ts.per.class)
      pdf.scaled.vaccination <- pmax(final.pdf,0)
      time.specific.routine[j,age.classes %in% c(age.min.MR2[j]:age.max.MR2[j])] <- pdf.scaled.vaccination*prob.success.MR2
      #time.specific.routine[j,age.classes %in% c(age.min.MR2[j]:age.max.MR2[j])] <- pdf*prob.success.MR2
      #IF assuming correlation, pdf*prob.success.MR2 = rep(0.1282051, 12); sum(rep(0.1282051, 12)) = 0.1538462 = 0.08/0.52 = (mcv2*ve2*(1-ve1)) / (1-ve1*mcv1), proportion of susceptibles successfully vaccinated with MCV2, gtg
      #IF assuming no correlation pdf*prob.success.MR2 = rep(0.333, 12); sum(rep(0.333, 12)) = 0.4 = 0.5*0.8 = mcv2*ve2, check
      one.minus.ve2[j] <- sum(pdf/sum(pdf)*(1-prob.success.MR2)) # = 0.2 = 1-ve2, check
      #prop.fail.MR2[j] <- sum(pdf*(1-prob.success.MR2))*prop.MR1failANDnoMR1.laggged[j] # 0.02, proportion of total population with primary failure from MR2, gtg
      prop.fail.MR2[j] <- sum(pdf*(1-prob.success.MR2)) # 0.02, proportion of total population with primary failure from MR2, gtg
      prop.fail.MR2.byage[j,age.classes %in% c(age.min.MR2[j]:age.max.MR2[j])] <- pdf*(1-prob.success.MR2)
    }
  }


  return(list(
    age.time.specific.routine=time.specific.routine,
    one.minus.ve1=one.minus.ve1,
    one.minus.ve2=one.minus.ve2,
    prop.fail.MR1=prop.fail.MR1,
    prop.fail.MR1.byage = prop.fail.MR1.byage,
    prop.fail.MR2=prop.fail.MR2,
    prop.fail.MR2.byage = prop.fail.MR2.byage))

}
