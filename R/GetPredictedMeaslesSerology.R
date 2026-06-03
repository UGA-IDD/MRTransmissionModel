#' Get model-predicted measles serology for observed serology rows
#'
#' Runs an MSIRV simulation (Part 1 transient burn-in, then Part 2 from
#' \code{year} to \code{year + t.max}) and returns predicted measles
#' seroprevalence aligned with every row of a serology data frame.
#'
#' @param serodata data.frame in long format with columns:
#'   \code{survey.time.point} (model time-step index, not calendar year),
#'   \code{age.bin.lower}, \code{age.bin.upper} (single-year age group bounds,
#'   inclusive), \code{n_tested}, and \code{n_positive}.
#' @param setup country setup object returned by a \code{setupCountry_*} helper.
#' @param year numeric. Simulation start year. Should remain 1980 because
#'   vaccination inputs are indexed from 1980.
#' @param t.max numeric. Number of years to simulate in Part 2.
#' @param R0 numeric. Basic reproduction number used to scale the WAIFW matrix.
#' @param rho numeric. Pearson correlation between RI and SIA doses,
#'   in \code{[-1, 1]} (default 1 = maximum correlation).
#' @param scale.sia numeric in \code{(0, 1]}. Multiplicative scaling factor
#'   applied to all SIA coverage values before simulation (default 1 = no
#'   scaling). Values below 1 reduce effective SIA coverage proportionally.
#' @param age.classes numeric vector. Upper bounds of age classes in months
#'   (default \code{c(1:240, seq(252, 1212, 12))}).
#' @param generation.time numeric. Generation time in months (default 0.5,
#'   giving 24 time steps per year).
#' @param seasonal.amp numeric. Amplitude of the seasonal cosine forcing
#'   (default 0.15).
#' @param age0is6to11monly logical. If \code{TRUE}, the age-0 (< 1 year)
#'   seroprevalence cell is computed from age classes 7-12 (months 6-11,
#'   0-indexed) only.
#'
#' @return The input \code{serodata} data.frame with three columns appended:
#'   \describe{
#'     \item{pred.seroprev}{model-predicted seroprevalence for that row}
#'     \item{pred.imm.pop}{predicted immune count (M + R + V) in that age bin
#'       at that time point}
#'     \item{pred.pop}{predicted total population in that age bin at that time
#'       point}
#'   }
#'
#' @details
#' \code{serodata$survey.time.point} must be expressed as a model time-step
#' index, not a calendar year. For a simulation starting in \code{year} with
#' \code{generation.time = 0.5} (24 steps/year), calendar year \code{y}
#' corresponds to time-step index \code{(y - year) * 24}.
#'
#' Seroprevalence is extracted only at the unique time points present in
#' \code{serodata$survey.time.point}, so post-processing cost scales with the
#' number of distinct survey years, not total simulation length.
#'
#' @seealso \code{\link{LogLikMeaslesSerology}}, \code{\link{FitMeaslesSerology}}
#'
#' @export
GetPredictedMeaslesSerology <- function(
    serodata,
    setup,
    year             = 1980,
    t.max,
    R0,
    rho              = 1,
    scale.sia        = 1,
    age.classes      = c(1:240, seq(252, 1212, 12)),
    generation.time  = 0.5,
    seasonal.amp     = 0.15,
    age0is6to11monly = FALSE) {

  req.cols <- c("survey.time.point", "age.bin.lower", "age.bin.upper",
                "n_tested", "n_positive")
  miss.cols <- setdiff(req.cols, names(serodata))
  if (length(miss.cols) > 0)
    stop("serodata is missing required columns: ", paste(miss.cols, collapse = ", "))

  if (any(serodata$age.bin.lower > serodata$age.bin.upper))
    stop("All age.bin.lower values must be <= age.bin.upper")

  if (any(serodata$survey.time.point < 1))
    stop("survey.time.point values must be >= 1")

  sia.cov <- setup$measlesSIA.coverage.1980to2100[(year - 1980 + 1):((year - 1980 + 1) + t.max)] * scale.sia

  med.pop.2020to2100 <- median(setup$pop.total.1950.2100[71:151])
  if (med.pop.2020to2100 > 70000000) {
    intro.rate.1950.2100 <- med.pop.2020to2100 * 0.015 / 1000000
  } else if (med.pop.2020to2100 < 5000000) {
    intro.rate.1950.2100 <- med.pop.2020to2100 * 0.05 / 1000000
  } else {
    intro.rate.1950.2100 <- 1
  }

  EXt0 <- EX.Country.part1(
    uncode = setup$uncode,
    generation.time = generation.time,
    age.classes = age.classes,
    maternal.decay.rt = 0.45,
    exponent = 0.97,
    yr.births.per.1000.acrossyears = rep(setup$cbr.1950.2100[year - 1950 + 1], 20),
    intro.rate = 1 / 24 / 320,
    targeted.intro = FALSE,
    R0 = R0,
    t.max = 20,
    get.births = setup$get.births.here,
    seasonal.amp = seasonal.amp,
    flat.WAIFW = FALSE,
    country.specific.WAIFW = TRUE,
    asdr.object = setup$asdr.object,
    year = year,
    use_montagu_demog = FALSE
  )

  vacc_succ_obj <- pvacsuccess(age.classes, get.boulianne.vsucc())
  vacc_succ_obj@prob.vsucc <- rep(1, length(vacc_succ_obj@ages))

  tmp.res <- EX.Country.part2(
    uncode = setup$uncode,
    pop.rescale = setup$pop.total.1950.2100[seq(1990, (year + t.max - 10), 10) - 1950 + 1],
    pop.time = seq(1990, (year + t.max - 10), 10) - year,
    is.stochastic = FALSE,
    get.births = setup$get.births.here,
    t.max = t.max,
    rescale.WAIFW = FALSE,
    yr.births.per.1000.acrossyears = setup$cbr.1950.2100[(year - 1950 + 1):((year - 1950 + 1) + t.max)],
    asdr.object = setup$asdr.object,
    year = year,
    EXt0 = EXt0,
    time.specific.MR1cov = setup$MCV1.coverage.1980to2100[(year - 1980 + 1):((year - 1980 + 1) + t.max)] * 0.85,
    time.specific.MR2cov = setup$MCV2.coverage.1980to2100[(year - 1980 + 1):((year - 1980 + 1) + t.max)] * 0.90,
    time.specific.SIAcov = sia.cov,
    time.specific.min.age.MR1 = rep(10, length(year:(year + t.max))),
    time.specific.max.age.MR1 = rep(12, length(year:(year + t.max))),
    time.specific.min.age.MR2 = rep(19, length(year:(year + t.max))),
    time.specific.max.age.MR2 = rep(25, length(year:(year + t.max))),
    time.specific.min.age.SIA = setup$age.min.sia.measles[(year - 1980 + 1):((year - 1980 + 1) + t.max)],
    time.specific.max.age.SIA = setup$age.max.sia.measles[(year - 1980 + 1):((year - 1980 + 1) + t.max)],
    obj.vcdf.MR1 = get.vcdf.uniform(10, 12),
    obj.vcdf.MR2 = get.vcdf.uniform(19, 25),
    obj.prob.vsucc = vacc_succ_obj,
    sia.timing.in.year = 1 / 12,
    MR1MR2correlation = TRUE,
    MR1SIAcorrelation = rho,
    MR2SIAcorrelation = rho,
    SIAinacc = FALSE,
    SIAinefficient = FALSE,
    intro.rate = intro.rate.1950.2100 / 24 / 320
  )

  tp.unique <- sort(unique(serodata$survey.time.point))

  pred.age <- GetMeaslesSeroprevalence.per.TimePoint(
    res              = tmp.res@result,
    trans            = tmp.res@experiment.def@trans,
    epi.state        = tmp.res@experiment.def@trans@epi.class,
    no.gens.in.year  = 24,
    time.point.index = tp.unique,
    age0is6to11monly = age0is6to11monly
  )

  bins.unique <- unique(serodata[, c("age.bin.lower", "age.bin.upper")])
  rownames(bins.unique) <- NULL

  pred.bin <- AggregateMeaslesSeroprevalenceByAgeBins(
    imm.pop       = pred.age$imm.pop,
    pop.age       = pred.age$pop.age,
    age.bin.lower = bins.unique$age.bin.lower,
    age.bin.upper = bins.unique$age.bin.upper
  )

  time.lookup   <- match(serodata$survey.time.point, tp.unique)
  bin.key.all   <- paste(serodata$age.bin.lower, serodata$age.bin.upper, sep = "_")
  bin.key.unique <- paste(bins.unique$age.bin.lower, bins.unique$age.bin.upper, sep = "_")
  bin.lookup    <- match(bin.key.all, bin.key.unique)

  serodata$pred.seroprev <- pred.bin$seroprev.bin[cbind(bin.lookup, time.lookup)]
  serodata$pred.imm.pop  <- pred.bin$imm.pop.bin[cbind(bin.lookup, time.lookup)]
  serodata$pred.pop      <- pred.bin$pop.bin[cbind(bin.lookup, time.lookup)]

  return(serodata)
}
