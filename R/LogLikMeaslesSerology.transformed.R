#' Binomial log-likelihood for measles serology calibration on transformed scale
#'
#' This function computes the binomial log-likelihood of observed measles
#' serological data using unconstrained transformed parameters. It is mainly
#' useful as a convenience wrapper when optimizing or sampling on the
#' transformed scale.
#'
#' @param theta numeric vector of transformed parameters:
#'   \describe{
#'     \item{theta[1]}{log-scale parameter for \code{R0}}
#'     \item{theta[2]}{atanh-scale parameter for \code{rho}}
#'   }
#' @param serodata data.frame containing observed serological data with columns:
#'   \describe{
#'     \item{survey.time.point}{model time point index for each observation}
#'     \item{age.bin.lower}{lower bound of age bin (inclusive)}
#'     \item{age.bin.upper}{upper bound of age bin (inclusive)}
#'     \item{n_tested}{number of individuals tested}
#'     \item{n_positive}{number of seropositive individuals}
#'   }
#' @param setup country setup object returned by a setupCountry_* helper.
#' @param year numeric. Simulation start year (typically 1980).
#' @param t.max numeric. Number of years to simulate.
#' @param ... additional arguments passed to \code{LogLikMeaslesSerology()}.
#'
#' @return numeric scalar. The total log-likelihood across all observations.
#'
#' @details
#' The transformed parameters are mapped to the natural scale using
#' \code{TransformMeaslesSerologyParameters()} and then passed to
#' \code{LogLikMeaslesSerology()}.
#'
#' @seealso \code{\link{TransformMeaslesSerologyParameters}},
#'   \code{\link{LogLikMeaslesSerology}},
#'   \code{\link{NegLogLikMeaslesSerology.transformed}}
#'
#' @export
LogLikMeaslesSerology.transformed <- function(theta, serodata, setup, year = 1980, t.max, ...){

  par.nat <- TransformMeaslesSerologyParameters(theta = theta)

  ll <- LogLikMeaslesSerology(
    R0 = par.nat$R0,
    rho = par.nat$rho,
    serodata = serodata,
    setup = setup,
    year = year,
    t.max = t.max,
    ...
  )

  return(ll)
}
