#' Wrapper function to add space to get.schoolvacc.time.age.specific
#'
#' @param time.step numeric, time step in months
#' @param age.classes vector, age classes in months, same as tran object
#' @param space.time.specific.cov matrix, space (rows) and year (columns) specific school vaccination coverage as a proportion
#' @param age.min vector, time specific minimum age eligible for school vaccination
#' @param age.max vector, time specific maximum age eligible for school vaccination
#' @param list.obj.vcdf list of vaccine.cdf.byage objects, each item is a subpop, age distribution of school vaccination
#' @param obj.prob.vsucc prob.vsucc.byage object, VE by age
#'
#' @return list, features specific to school vaccination, including age.time.specific.school.vaccination matrix which is age and space (rows) by year (cols) of proportion school vaccination
#' @export
#'
space.wrapper.get.schoolvacc.time.age.specific <- function(time.step=0.5, age.classes = c(1:240, seq(252,1212,12)),
                                                                   space.time.specific.cov = rbind(rep(0.7, 50),rep(0.6, 50)), age.min=rep(12, 50), age.max=rep(23,50),
                                                                   list.obj.vcdf = replicate(2,list(get.vcdf.uniform(72, 83))),
                                                                   obj.prob.vsucc = pvacsuccess(1:(14*12), new("vsucc.constant", success.rate=0.8))){
  n.subpops <- nrow(space.time.specific.cov)
  for (s in 1:n.subpops){
    tmp <- get.schoolvacc.time.age.specific(time.step, age.classes,
                                            time.specific.cov=space.time.specific.cov[s,],
                                            age.min, age.max,
                                            obj.vcdf = list.obj.vcdf[[s]],
                                            obj.prob.vsucc)

    #we need to transpose these matrices so that age (rows) and time (cols)
    tmp$age.time.specific.schoolvacc <- t(tmp$age.time.specific.schoolvacc)
    tmp$prop.fail.schoolvacc.byage <- t(tmp$prop.fail.schoolvacc.byage)

    if (s==1) {
      out <- tmp
    }

    # row bind so that each new sub-population adds length(age.classes) number of rows to the bottom of the matrix
    if (s!=1){
      out$age.time.specific.schoolvacc <- rbind(out$age.time.specific.schoolvacc, tmp$age.time.specific.schoolvacc)
      out$prop.fail.schoolvacc.byage <- rbind(out$prop.fail.schoolvacc.byage, tmp$prop.fail.school.vaccination.byage)
      out$one.minus.ve.schoolvacc <- rbind(out$one.minus.ve.schoolvacc, tmp$prop.fail.schoolvacc.byage)
      out$prop.fail.schoolvacc <- rbind(out$prop.fail.schoolvacc, tmp$prop.fail.schoolvacc)
    }
  }
  return(out)
}
