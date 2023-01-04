#' Wrapper function to add space to get.routine.time.age.specific
#'
#' @param time.step numeric, time step in months
#' @param age.classes vector, age classes in months, same as tran object
#' @param space.time.specific.MR1cov matrix, space (rows) and year (columns) specific routine MR1 coverage as a proportion
#' @param age.min.MR1 vector, time specific minimum age eligible for MR1
#' @param age.max.MR1 vector, time specific maximum age eligible for MR1
#' @param space.time.specific.MR2cov matrix, space (rows) and year (columns) specific routine MR2 coverage as a proportion
#' @param age.min.MR2 vector, time specific minimum age eligible for MR2
#' @param age.max.MR2 vector, time specific maximum age eligible for MR2
#' @param obj.vcdf.MR1 vaccine.cdf.byage object, age distribution of routine MR1
#' @param obj.vcdf.MR2 vaccine.cdf.byage object, age distribution of routine MR2
#' @param obj.prob.vsucc prob.vsucc.byage object, VE by age relevant to both routine vacccine doses
#' @param MR1MR2correlation boolean, if MR1 and MR2 are dependent (TRUE) or independent (FALSE)
#'
#' @return list, features specific to routine vaccination, including age.time.specific.routine matrix which is age and space (rows) by year (cols) of proportion routine vaccination
#' @export
#'
space.wrapper.get.routine.time.age.specific <- function(time.step=0.5, age.classes = c(1:240, seq(252,1212,12)),
                                                        space.time.specific.MR1cov = rbind(rep(0.7, 50),rep(0.6, 50)), age.min.MR1=rep(12, 50), age.max.MR1=rep(23,50),
                                                        space.time.specific.MR2cov = rbind(rep(0.4, 50),rep(0.3, 50)), age.min.MR2=rep(24, 50), age.max.MR2=rep(35, 50),
                                                        obj.vcdf.MR1=get.vcdf.uniform(12, 23), obj.vcdf.MR2=get.vcdf.uniform(24, 35),
                                                        obj.prob.vsucc = pvacsuccess(1:(14*12), new("vsucc.constant", success.rate=0.8)),
                                                        MR1MR2correlation=F){
  n.subpops <- nrow(space.time.specific.MR1cov)
  for (s in 1:n.subpops){
    tmp <- get.routine.time.age.specific(time.step, age.classes,
                                         time.specific.MR1cov=space.time.specific.MR1cov[s,], age.min.MR1, age.max.MR1,
                                         time.specific.MR2cov=space.time.specific.MR2cov[s,], age.min.MR2, age.max.MR2,
                                         obj.vcdf.MR1, obj.vcdf.MR2,
                                         obj.prob.vsucc, MR1MR2correlation)

    #we need to transpose these matrices so that age (rows) and time (cols)
    tmp$age.time.specific.routine <- t(tmp$age.time.specific.routine)
    tmp$prop.fail.MR1.byage <- t(tmp$prop.fail.MR1.byage)
    tmp$prop.fail.MR2.byage <- t(tmp$prop.fail.MR2.byage)

    if (s==1) {
      out <- tmp
    }

    # row bind so that each new sub-population adds length(age.classes) number of rows to the bottom of the matrix
    if (s!=1){
      out$age.time.specific.routine <- rbind(out$age.time.specific.routine, tmp$age.time.specific.routine)
      out$prop.fail.MR1.byage <- rbind(out$prop.fail.MR1.byage, tmp$prop.fail.MR1.byage)
      out$prop.fail.MR2.byage <- rbind(out$prop.fail.MR2.byage, tmp$prop.fail.MR2.byage)
      out$one.minus.ve1 <- rbind(out$one.minus.ve1, tmp$one.minus.ve1)
      out$one.minus.ve2 <- rbind(out$one.minus.ve2, tmp$one.minus.ve2)
      out$prop.fail.MR1 <- rbind(out$prop.fail.MR1, tmp$prop.fail.MR1)
      out$prop.fail.MR2 <- rbind(out$prop.fail.MR2, tmp$prop.fail.MR2)
    }
  }
  return(out)
}
