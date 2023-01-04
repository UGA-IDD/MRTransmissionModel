#' Wrapper function to add space to get.sia.time.age.specific
#'
#' @param age.classes vector, age classes in months, same as tran object
#' @param space.time.specific.SIAcov matrix, space (rows) and year (columns) specific campaign MR coverage as a proportion
#' @param age.min.sia vector, time specific minimum age eligible for MR campaign dose
#' @param age.max.sia vector, time specific maximum age eligible for MR campaign dose
#' @param obj.prob.vsucc prob.vsucc.byage object, VE by age relevant to campaign vacccine doses
#'
#' @return list, features specific to campaign vaccination, including age.time.specific.SIA matrix which is age and space (rows) by year (cols) of proportion campaign vaccination
#' @export
#'
space.wrapper.get.sia.time.age.specific <- function(age.classes=c(1:240, seq(252,1212,12)),
                                                    space.time.specific.SIAcov=rbind(c(rep(0,29),0.5,rep(0,19),0.7),c(rep(0,29),0.7,rep(0,19),0.8)),
                                                    age.min.sia=rep(12,50),
                                                    age.max.sia=rep(60,50),
                                                    obj.prob.vsucc=pvacsuccess(1:240, new("vsucc.constant", success.rate=0.85))){
  n.subpops <- nrow(space.time.specific.SIAcov)
  for (s in 1:n.subpops){
    tmp <- get.sia.time.age.specific(age.classes=age.classes,
                                     time.specific.SIAcov=space.time.specific.SIAcov[s,],
                                     age.min.sia=age.min.sia,
                                     age.max.sia=age.max.sia,
                                     obj.prob.vsucc=obj.prob.vsucc)

    #we need to transpose these matrices so that age (rows) and time (cols)
    tmp$age.time.specific.SIA <- t(tmp$age.time.specific.SIA)

    if (s==1) {
      out <- tmp
    }

    # row bind so that each new sub-population adds length(age.classes) number of rows to the bottom of the matrix
    if (s!=1){
      out$age.time.specific.SIA <- rbind(out$age.time.specific.SIA, tmp$age.time.specific.SIA)
      out$prop.fail.SIA <- rbind(out$prop.fail.SIA, tmp$prop.fail.SIA)
    }
  }
  return(out)
}
