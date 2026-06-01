#' Binomial log-likelihood for measles serology calibration
#'
#' Runs a full MSIRV simulation for the given parameters and computes the
#' binomial log-likelihood of the observed serological data.
#'
#' @param R0 numeric. Basic reproduction number used to scale the WAIFW matrix.
#' @param rho numeric. Pearson correlation between RI and SIA doses,
#'   in \code{[-1, 1]} (default 0 = independence).
#' @param serodata data.frame with columns \code{survey.time.point},
#'   \code{age.bin.lower}, \code{age.bin.upper}, \code{n_tested},
#'   \code{n_positive}.
#' @param setup country setup object returned by a \code{setupCountry_*} helper.
#' @param year numeric. Simulation start year (default 1980).
#' @param t.max numeric. Number of years to simulate.
#' @param age.classes numeric vector. Upper bounds of age classes in months.
#' @param generation.time numeric. Generation time in months (default 0.5).
#' @param seasonal.amp numeric. Seasonal forcing amplitude (default 0.15).
#' @param age0is6to11monly logical. Passed as \code{age0is6to11monly} to
#'   \code{GetPredictedMeaslesSerology()}.
#' @param eps numeric. Small value used to bound predicted probabilities away
#'   from 0 and 1 before evaluating the log-likelihood (default 1e-10).
#' @param .cache environment or \code{NULL}. Passed to
#'   \code{GetPredictedMeaslesSerology()} to enable Part 1 result caching by
#'   \code{R0}. See \code{\link{GetPredictedMeaslesSerology}} for details.
#'
#' @return numeric scalar. Total binomial log-likelihood across all rows of
#'   \code{serodata}.
#' @export
LogLikMeaslesSerology <- function(
    R0,
    rho,
    serodata,
    setup,
    year = 1980,
    t.max,
    age.classes = c(1:240, seq(252, 1212, 12)),
    generation.time = 0.5,
    seasonal.amp = 0.15,
    age0is6to11monly = FALSE,
    eps = 1e-10,
    .cache = NULL){

  pred <- GetPredictedMeaslesSerology(
    serodata = serodata,
    setup = setup,
    year = year,
    t.max = t.max,
    R0 = R0,
    rho = rho,
    age.classes = age.classes,
    generation.time = generation.time,
    seasonal.amp = seasonal.amp,
    age0is6to11monly = age0is6to11monly,
    .cache = .cache
  )

  p <- pred$pred.seroprev
  p <- pmax(eps, pmin(1 - eps, p))

  ll <- sum(dbinom(
    x = pred$n_positive,
    size = pred$n_tested,
    prob = p,
    log = TRUE
  ))

  return(ll)
}
