#' Function to Run Outbreak Response Simulations/Experiments
#'
#' @param uncode UN country code
#' @param generation.time generation time in months
#' @param pop.rescale numeric; population by which you want to rescale at pop.time
#' @param pop.time numeric; time in YEARS that you want to rescale the population
#' @param is.stochastic logical
#' @param t.max numeric; time in years to run the experiment
#' @param rescale.WAIFW logical
#' @param yr.births.per.1000.acrossyears vector of length t.max; crude birth rate per 1000 per year
#' @param asdr.object nMx object; age specific death rates
#' @param year numeric; simulation start year
#' @param EXt0 object from EX.Country.part1()
#' @param time.specific.MR1cov numeric vector; MR1 coverage by year
#' @param time.specific.MR2cov numeric vector; MR2 coverage by year
#' @param time.specific.SIAcov numeric vector; SIA coverage by year
#' @param time.specific.min.age.MR1 numeric vector; minimum age for MR1 by year (months)
#' @param time.specific.max.age.MR1 numeric vector; maximum age for MR1 by year (months)
#' @param time.specific.min.age.MR2 numeric vector; minimum age for MR2 by year (months)
#' @param time.specific.max.age.MR2 numeric vector; maximum age for MR2 by year (months)
#' @param time.specific.min.age.SIA numeric vector; minimum age for SIA by year (months)
#' @param time.specific.max.age.SIA numeric vector; maximum age for SIA by year (months)
#' @param obj.vcdf.MR1 vaccine.cdf.byage object for MR1
#' @param obj.vcdf.MR2 vaccine.cdf.byage object for MR2
#' @param obj.prob.vsucc prob.vsucc.byage object
#' @param sia.timing.in.year numeric; fraction of year at which SIA occurs
#' @param MR1MR2correlation logical; whether MR1 and MR2 doses are correlated
#' @param MR1SIAcorrelation numeric in [0,1]; Pearson rho between MR1 and SIA doses
#' @param MR2SIAcorrelation numeric in [0,1]; Pearson rho between MR2 and SIA doses
#' @param intro.rate numeric vector; introduction rate per age class per time step
#' @param or.reporting.rate numeric; reporting rate (care-seeking x clinical recognition, pre-testing)
#'   in (0,1]; scalar or vector of length n.age. Default 1.
#' @param or.non.meas.cases.by.age.month matrix; non-measles suspected cases with rows = model
#'   age classes and columns = 12 calendar months.
#' @param or.Se numeric; diagnostic test sensitivity in [0,1].
#' @param or.Sp numeric; diagnostic test specificity in [0,1].
#' @param or.n.confirmations.target numeric; trigger threshold (sum of obs.TP over window must reach
#'   this) and per-step testing stop rule (default 5).
#' @param or.vacc.agedist.percentile numeric in (0,1]; percentile of cumulative case age CDF used as
#'   the campaign upper age bound.
#' @param or.confirmation.delay numeric; time steps from end of surveillance window to case confirmation
#' @param or.response.delay numeric; time steps from confirmation to campaign delivery
#' @param or.trigger.window numeric; number of time steps to sum obs.TP over for trigger metric
#' @param or.trigger.window.age.lower numeric; lower age bound for surveillance trigger (months); NA = all ages
#' @param or.trigger.window.age.upper numeric; upper age bound for surveillance trigger (months); NA = all ages
#' @param or.vacc.coverage numeric in [0,1]; OBR campaign coverage
#' @param or.min.interval numeric; minimum time steps between successive OBR triggers
#' @param or.start.timestep numeric; earliest time step at which OBR can trigger (1-indexed); default 1
#' @param stop_testing_at_trigger logical; if TRUE, campaign age targeting uses suspected cases
#'   (true.reported.all + non-measles background; reporting process only); if FALSE, uses obs.TP
#'   accumulated through the response window (diagnostic confirmation continues)
#'
#' @include setClasses.R
#' @importFrom methods new
#'
#' @return experiment results
#' @export

EX.Country.part2.OR <- function(uncode,
                                generation.time          = 0.5,
                                pop.rescale              = NULL,
                                pop.time                 = 4,
                                is.stochastic            = FALSE,
                                t.max                    = 40,
                                rescale.WAIFW            = FALSE,
                                yr.births.per.1000.acrossyears,
                                asdr.object,
                                year                     = 1990,
                                EXt0,
                                time.specific.MR1cov,
                                time.specific.MR2cov,
                                time.specific.SIAcov,
                                time.specific.min.age.MR1,
                                time.specific.max.age.MR1,
                                time.specific.min.age.MR2,
                                time.specific.max.age.MR2,
                                time.specific.min.age.SIA,
                                time.specific.max.age.SIA,
                                obj.vcdf.MR1             = get.vcdf.normal(6, 12),
                                obj.vcdf.MR2             = get.vcdf.normal(15, 21),
                                obj.prob.vsucc           = pvacsuccess(1:(20*12), get.boulianne.vsucc()),
                                sia.timing.in.year       = (3/12),
                                MR1MR2correlation        = TRUE,
                                MR1SIAcorrelation        = 1,
                                MR2SIAcorrelation        = 1,
                                intro.rate,
                                or.reporting.rate              = 1,
                                or.non.meas.cases.by.age.month,
                                or.Se,
                                or.Sp,
                                or.n.confirmations.target      = 5,
                                or.vacc.agedist.percentile,
                                or.confirmation.delay,
                                or.response.delay,
                                or.trigger.window,
                                or.trigger.window.age.lower    = NA_real_,
                                or.trigger.window.age.upper    = NA_real_,
                                or.vacc.coverage,
                                or.min.interval,
                                or.start.timestep              = 1,
                                stop_testing_at_trigger        = FALSE) {

  EX <- new("experiment.updatedemog.vaccinationchange.vaccinationcorrelation.outbreakresponse")

  name <- countrycode::countrycode(uncode, origin = "un", destination = "country.name")

  EX@trans        <- EXt0@trans
  EX@state.t0     <- EXt0@state.t0
  EX@maternal.obj <- EXt0@maternal.obj
  EX@name         <- paste("OBR Experiment:", name, year, "to", (year + t.max), sep = " ")
  EX@description  <- "Frequency Dependent Stochastic"
  EX@t.min        <- 0
  EX@t0.doy       <- 0
  EX@step.size    <- EXt0@step.size
  EX@season.obj   <- EXt0@season.obj
  EX@t.max        <- t.max

  no.time.steps.in.experiment <- round((EX@t.max - EX@t.min) / EX@step.size) + 1
  no.gens.in.year <- (no.time.steps.in.experiment - 1) / t.max

  EX@trans@is.stochastic <- is.stochastic
  if (is.stochastic) EX@state.t0[, 1] <- round(EX@state.t0[, 1])

  EX@pop.rescale.each.timestep <- rep(NaN, no.time.steps.in.experiment)
  if (!is.null(pop.rescale)) {
    for (p in 1:length(pop.rescale)) {
      EX@pop.rescale.each.timestep[pop.time[p] * no.gens.in.year - 1] <- pop.rescale[p]
    }
  }

  EX@births.per.1000.each.timestep <- c(yr.births.per.1000.acrossyears[1],
                                         rep(yr.births.per.1000.acrossyears,
                                             each = no.gens.in.year)) * generation.time / 12

  EX@surv.each.timestep <- create.surv.prob.over.age.time(
    EX@trans@age.class, generation.time,
    nMx       = asdr.object,
    nMx.years = seq(year, (year + t.max - 1), 1))

  EX@time.specific.MR1cov      <- time.specific.MR1cov
  EX@time.specific.MR2cov      <- time.specific.MR2cov
  EX@time.specific.SIAcov      <- time.specific.SIAcov
  EX@time.specific.min.age.MR1 <- time.specific.min.age.MR1
  EX@time.specific.max.age.MR1 <- time.specific.max.age.MR1
  EX@time.specific.min.age.MR2 <- time.specific.min.age.MR2
  EX@time.specific.max.age.MR2 <- time.specific.max.age.MR2
  EX@time.specific.min.age.SIA <- time.specific.min.age.SIA
  EX@time.specific.max.age.SIA <- time.specific.max.age.SIA
  EX@obj.vcdf.MR1              <- obj.vcdf.MR1
  EX@obj.vcdf.MR2              <- obj.vcdf.MR2
  EX@obj.prob.vsucc            <- obj.prob.vsucc
  EX@sia.timing.in.year        <- sia.timing.in.year
  EX@MR1MR2correlation         <- MR1MR2correlation
  EX@MR1SIAcorrelation         <- as.numeric(MR1SIAcorrelation)
  EX@MR2SIAcorrelation         <- as.numeric(MR2SIAcorrelation)

  if (length(intro.rate) > 1) {
    EX@intro.rate <- intro.rate
  } else {
    EX@intro.rate <- EXt0@trans@introduction.rate
  }

  or.total.delay                     <- or.confirmation.delay + or.response.delay
  EX@or.reporting.rate               <- or.reporting.rate
  EX@or.non.meas.cases.by.age.month  <- or.non.meas.cases.by.age.month
  EX@or.Se                           <- or.Se
  EX@or.Sp                           <- or.Sp
  EX@or.n.confirmations.target       <- or.n.confirmations.target
  EX@or.vacc.agedist.percentile      <- or.vacc.agedist.percentile
  EX@or.total.delay                  <- or.total.delay
  EX@or.response.delay               <- or.response.delay
  EX@or.trigger.window               <- or.trigger.window
  EX@or.trigger.window.age.lower     <- or.trigger.window.age.lower
  EX@or.trigger.window.age.upper     <- or.trigger.window.age.upper
  EX@or.vacc.coverage                <- or.vacc.coverage
  EX@or.min.interval                 <- or.min.interval
  EX@or.start.timestep               <- or.start.timestep
  EX@stop_testing_at_trigger         <- stop_testing_at_trigger

  return(run(EX, rescale.WAIFW = rescale.WAIFW))
}
