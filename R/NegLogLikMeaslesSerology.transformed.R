#' Negative log-likelihood for measles serology calibration on transformed scale
#'
#' This function computes the negative binomial log-likelihood of observed
#' measles serological data using unconstrained transformed parameters. It is
#' intended for use with optimization routines such as \code{optim()}.
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
#' @return numeric scalar. The negative log-likelihood.
#'
#' @details
#' This function transforms unconstrained parameters to the natural scale using
#' \code{TransformMeaslesSerologyParameters()}, evaluates the log-likelihood via
#' \code{LogLikMeaslesSerology()}, and returns its negative. This allows
#' optimization to proceed on an unconstrained scale while ensuring that:
#' \itemize{
#'   \item \code{R0 > 0}
#'   \item \code{-1 < rho < 1}
#' }
#'
#' @seealso \code{\link{TransformMeaslesSerologyParameters}},
#'   \code{\link{LogLikMeaslesSerology.transformed}},
#'   \code{\link{LogLikMeaslesSerology}}
#'
#' @export
NegLogLikMeaslesSerology.transformed <- local({

  iter <- 0  # persists across calls

  function(theta, serodata, setup, year = 1980, t.max, ...){

    iter <<- iter + 1

    R0  <- exp(theta[1])
    rho <- tanh(theta[2])

    cat(sprintf("Iter %d | R0 = %.3f | rho = %.3f\n",
                iter, R0, rho))

    ll <- LogLikMeaslesSerology(
      R0 = R0,
      rho = rho,
      serodata = serodata,
      setup = setup,
      year = year,
      t.max = t.max,
      ...
    )

    return(-ll)
  }
})
