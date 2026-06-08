#' Combined log-likelihood for measles serology and age-stratified cases
#'
#' Runs a single MSIRV simulation and computes the total log-likelihood as the
#' sum of the binomial log-likelihood for serological data and the
#' negative-binomial log-likelihood for age-stratified case counts.
#'
#' @param R0 numeric. Basic reproduction number.
#' @param rho numeric. Pearson correlation between RI and SIA doses,
#'   in \code{[-1, 1]} (default 1 = maximum correlation).
#' @param scale.sia numeric in \code{(0, 1]}. Multiplicative scaling factor
#'   applied to all SIA coverage values (default 1 = no scaling).
#' @param phi numeric in \code{(0, 1)}. Case reporting fraction.
#' @param kappa numeric \code{> 0}. Negative-binomial dispersion parameter.
#' @param serodata data.frame with columns \code{survey.time.point},
#'   \code{age.bin.lower}, \code{age.bin.upper}, \code{n_tested},
#'   \code{n_positive}.
#' @param casedata data.frame with columns \code{age.lower} (integer years,
#'   inclusive lower bound of age band), \code{age.upper} (integer years,
#'   inclusive upper bound), \code{year} (calendar year), and \code{cases}
#'   (observed count). Pre-aggregate to 5-year bands before passing.
#' @param setup country setup object returned by a \code{setupCountry_*} helper.
#' @param year numeric. Simulation start year (default 1980).
#' @param t.max numeric. Number of years to simulate.
#' @param age.classes numeric vector. Upper bounds of age classes in months
#'   (default \code{c(1:240, seq(252, 1212, 12))}).
#' @param generation.time numeric. Generation time in months (default 0.5).
#' @param seasonal.amp numeric. Seasonal forcing amplitude (default 0.15).
#' @param age0is6to11monly logical. If \code{TRUE}, the age-0 seroprevalence
#'   cell uses months 7-12 only (default \code{FALSE}).
#' @param eps numeric. Small value used to bound serology predicted
#'   probabilities away from 0 and 1, and to floor expected case counts
#'   (default \code{1e-10}).
#'
#' @return numeric scalar. Total log-likelihood (serology + cases).
#'
#' @seealso \code{\link{GetPredictedMeaslesSerologyAndCases}},
#'   \code{\link{LogLikMeaslesCases}},
#'   \code{\link{FitMeaslesSerologyAndCases}}
#'
#' @export
LogLikMeaslesSerologyAndCases <- function(
    R0,
    rho              = 1,
    scale.sia        = 1,
    phi,
    kappa,
    serodata,
    casedata,
    setup,
    year             = 1980,
    t.max,
    age.classes      = c(1:240, seq(252, 1212, 12)),
    generation.time  = 0.5,
    seasonal.amp     = 0.15,
    age0is6to11monly = FALSE,
    eps              = 1e-10) {

  pred <- GetPredictedMeaslesSerologyAndCases(
    serodata         = serodata,
    casedata         = casedata,
    setup            = setup,
    year             = year,
    t.max            = t.max,
    R0               = R0,
    rho              = rho,
    scale.sia        = scale.sia,
    age.classes      = age.classes,
    generation.time  = generation.time,
    seasonal.amp     = seasonal.amp,
    age0is6to11monly = age0is6to11monly
  )

  p <- pred$serodata$pred.seroprev
  p <- pmax(eps, pmin(1 - eps, p))

  ll.sero <- sum(dbinom(
    x    = pred$serodata$n_positive,
    size = pred$serodata$n_tested,
    prob = p,
    log  = TRUE
  ))

  ll.case <- LogLikMeaslesCases(
    casedata = pred$casedata,
    phi      = phi,
    kappa    = kappa,
    eps      = eps
  )

  return(ll.sero + ll.case)
}
