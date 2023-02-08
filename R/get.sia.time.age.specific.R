#' Get SIA time age-specific
#'
#' @param age.classes xxx
#' @param time.specific.SIAcov xxx
#' @param age.min.sia xxx
#' @param age.max.sia xxx
#' @param obj.prob.vsucc xxx
#'
#' @return xxx
#' @export
#'

get.sia.time.age.specific <- function(age.classes=c(1:240, seq(252,1212,12)),
                                      time.specific.SIAcov=c(rep(0,29),0.5,rep(0,19),0.7),
                                      age.min.sia=rep(12,50),
                                      age.max.sia=rep(60,50),
                                      obj.prob.vsucc=pvacsuccess(1:240, new("vsucc.constant", success.rate=0.85))){

  ## Get SIA age specific estimates for each year
  age.time.specific.SIA <-  matrix(0, length(time.specific.SIAcov), length(age.classes))
  prop.fail.SIA <- rep(0, length(time.specific.SIAcov))

  for (j in 1:length(time.specific.SIAcov)){
    if (time.specific.SIAcov[j]!=0){
      prob.success.SIA <- obj.prob.vsucc@prob.vsucc[obj.prob.vsucc@ages %in% c(age.min.sia[j]:age.max.sia[j])]
      age.min.index <- which(age.classes==age.min.sia[j])
      age.max.index <-  which(age.classes==age.max.sia[j])
      vacc.SIA <- rep(time.specific.SIAcov[j], length(age.min.index:age.max.index))
      age.time.specific.SIA[j,(age.min.index:age.max.index)] <- vacc.SIA*prob.success.SIA
      prop.fail.SIA[j] <- (vacc.SIA*(1-prob.success.SIA))[1]
    }
  }

  return(list(age.time.specific.SIA=age.time.specific.SIA,
              prop.fail.SIA=prop.fail.SIA))

}
