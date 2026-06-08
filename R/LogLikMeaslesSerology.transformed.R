#' Binomial log-likelihood for measles serology calibration on transformed scale
#'
#' Convenience wrapper that accepts unconstrained transformed parameters,
#' back-transforms them, and returns the binomial log-likelihood.
#'
#' @param theta numeric vector of transformed parameters. Length 2 when
#'   \code{fix.rho} is a number (\code{log(R0)}, \code{qlogis(scale.sia)}),
#'   or length 3 when \code{fix.rho = NA} (adds \code{atanh(rho)} as third
#'   element).
#' @param serodata data.frame containing observed serological data with columns:
#'   \describe{
#'     \item{survey.time.point}{model time point index for each observation}
#'     \item{age.bin.lower}{lower bound of age bin (inclusive)}
#'     \item{age.bin.upper}{upper bound of age bin (inclusive)}
#'     \item{n_tested}{number of individuals tested}
#'     \item{n_positive}{number of seropositive individuals}
#'   }
#' @param setup country setup object returned by a setupCountry_* helper.
#' @param year numeric. Simulation start year (default 1980).
#' @param t.max numeric. Number of years to simulate.
#' @param fix.rho numeric scalar or \code{NA}. When a number, \code{rho} is
#'   fixed at that value and \code{theta} must have length 2. When \code{NA},
#'   \code{rho} is a free parameter in \code{theta[3]}. Default \code{1}.
#' @param age.classes numeric vector. Upper bounds of age classes in months
#'   (default \code{c(1:240, seq(252, 1212, 12))}).
#' @param generation.time numeric. Generation time in months (default 0.5).
#' @param seasonal.amp numeric. Seasonal forcing amplitude (default 0.15).
#' @param age0is6to11monly logical. If \code{TRUE}, the age-0 seroprevalence
#'   cell uses months 7-12 only (default \code{FALSE}).
#' @param eps numeric. Small value to bound predicted probabilities away from
#'   0 and 1 (default 1e-10).
#'
#' @return numeric scalar. The total log-likelihood across all observations.
#'
#' @seealso \code{\link{TransformMeaslesSerologyParameters}},
#'   \code{\link{LogLikMeaslesSerology}},
#'   \code{\link{NegLogLikMeaslesSerology.transformed}}
#'
#' @export
LogLikMeaslesSerology.transformed <- function(
    theta,
    serodata,
    setup,
    year             = 1980,
    t.max,
    fix.rho          = 1,
    age.classes      = c(1:240, seq(252, 1212, 12)),
    generation.time  = 0.5,
    seasonal.amp     = 0.15,
    age0is6to11monly = FALSE,
    eps              = 1e-10) {

  par.nat <- TransformMeaslesSerologyParameters(theta = theta, fix.rho = fix.rho)

  ll <- LogLikMeaslesSerology(
    R0               = par.nat$R0,
    rho              = par.nat$rho,
    scale.sia        = par.nat$scale.sia,
    serodata         = serodata,
    setup            = setup,
    year             = year,
    t.max            = t.max,
    age.classes      = age.classes,
    generation.time  = generation.time,
    seasonal.amp     = seasonal.amp,
    age0is6to11monly = age0is6to11monly,
    eps              = eps
  )

  return(ll)
}
