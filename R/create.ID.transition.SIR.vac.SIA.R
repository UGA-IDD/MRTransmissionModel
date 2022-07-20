#' Method to create an ID.transtion.SIR.vac object
#'
#' @param n.age.class xxx
#' @param aging.rate xxx
#' @param survival.rate xxx
#' @param waifw xxx
#' @param routine.vac xxx
#' @param routine.vac.age.index xxx
#' @param time.step xxx
#' @param age.class xxx
#' @param n.epi.class xxx
#' @param birth.rate xxx
#' @param sia.vac xxx
#' @param sia.vsucc xxx
#' @param introduction.rate xxx
#' @param exponent xxx
#' @param frequency.dep xxx
#' @param is.stochastic xxx
#' @param get.births xxx
#'
#' @importFrom methods new
#'
#' @return transition matrix
#' @export
#'

create.ID.transition.SIR.vac.SIA <- function(n.age.class,
                                             aging.rate,
                                             survival.rate,
                                             waifw,
                                             routine.vac = 0, #xxamy
                                             routine.vac.age.index = 12, #xxamy age index based on age classes
                                             time.step = 1/2, #time step in months == generation time
                                             age.class = numeric(),
                                             n.epi.class = 3,
                                             birth.rate = 0,
                                             sia.vac = 0,
                                             sia.vsucc = new("vsucc.constant", success.rate=1),
                                             introduction.rate = 0,
                                             exponent = 1,
                                             frequency.dep=F,
                                             is.stochastic=F,
                                             get.births = function(pop,tran){return(c(tran@birth.rate,rep(0,length(pop)-1)))}) {

  tmp <- create.ID.transition.SIR.vac(n.age.class = n.age.class,
                                      aging.rate = aging.rate,
                                      survival.rate = survival.rate,
                                      waifw = waifw,
                                      age.class = age.class,
                                      n.epi.class = n.epi.class,
                                      birth.rate = birth.rate,
                                      introduction.rate = introduction.rate,
                                      exponent = exponent,
                                      frequency.dep=frequency.dep,
                                      is.stochastic=is.stochastic,
                                      get.births = get.births)

  #get the %vaccinated in each age group per time step
  #time step must be in years!
  #vac.per <- create.vacc.per.time.step.by.age(vpdf, age.class, time.step) #xxamy
  vac.per <- new("vacc.per.time.step.by.age", #xxamy
                 pvacc.in.age.class = numeric(length=n.age.class))
  vac.per[routine.vac.age.index] <- routine.vac #xxamy

  rc <- new("ID.transition.SIR.vac.SIA",
            n.epi.class = tmp@n.epi.class,
            epi.class = tmp@epi.class,
            s.inds = tmp@s.inds,
            i.inds = tmp@i.inds,
            r.inds = tmp@r.inds,
            epi.class.label = c("S","I","R"),
            n.age.class = tmp@n.age.class,
            aging.rate = tmp@aging.rate,
            survival.rate = tmp@survival.rate,
            birth.rate = tmp@birth.rate,
            waifw = tmp@waifw,
            age.surv.matrix = tmp@age.surv.matrix,
            vac.per = vac.per,
            sia.vac = sia.vac,
            sia.vsucc = sia.vsucc,
            introduction.rate = introduction.rate,
            exponent = exponent,
            age.class = tmp@age.class,
            frequency.dep=frequency.dep,
            is.stochastic=is.stochastic,
            get.births = get.births
  )

  return(rc)
}
