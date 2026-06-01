#' Transform unconstrained calibration parameters to natural scale
#'
#' This helper maps unconstrained optimization or MCMC parameters to the
#' natural parameter space used by the measles serology calibration model.
#'
#' @param theta numeric vector of transformed parameters:
#'   \describe{
#'     \item{theta[1]}{log-scale parameter for \code{R0}}
#'     \item{theta[2]}{atanh-scale parameter for \code{rho}}
#'   }
#'
#' @return A named list containing:
#' \describe{
#'   \item{R0}{basic reproduction number, constrained to \code{R0 > 0}}
#'   \item{rho}{RI-SIA dose correlation, constrained to \code{-1 < rho < 1}}
#' }
#'
#' @details
#' The transformations are:
#' \deqn{R0 = \exp(\theta_1)}
#' \deqn{\rho = \tanh(\theta_2)}
#'
#' @examples
#' TransformMeaslesSerologyParameters(c(log(12), atanh(0.3)))
#'
#' @export
TransformMeaslesSerologyParameters <- function(theta){

  if(length(theta) != 2){
    stop("theta must be a numeric vector of length 2")
  }

  return(list(
    R0  = exp(theta[1]),
    rho = tanh(theta[2])
  ))
}
