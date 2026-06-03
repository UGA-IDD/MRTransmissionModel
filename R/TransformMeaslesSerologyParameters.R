#' Transform unconstrained calibration parameters to natural scale
#'
#' Maps unconstrained optimization or MCMC parameters to the natural parameter
#' space used by the measles serology calibration model.
#'
#' @param theta numeric vector of length 2 (when \code{fix.rho} is a number) or
#'   length 3 (when \code{fix.rho = NA}):
#'   \describe{
#'     \item{theta[1]}{log-scale parameter for \code{R0}}
#'     \item{theta[2]}{logit-scale parameter for \code{scale.sia}}
#'     \item{theta[3]}{atanh-scale parameter for \code{rho} (only when
#'       \code{fix.rho = NA})}
#'   }
#' @param fix.rho numeric scalar or \code{NA}. When a number, \code{rho} is
#'   fixed at that value and \code{theta} must have length 2. When \code{NA},
#'   \code{rho} is treated as a free parameter encoded in \code{theta[3]} on the
#'   atanh scale, and \code{theta} must have length 3. Default \code{1}.
#'
#' @return A named list containing:
#' \describe{
#'   \item{R0}{basic reproduction number, constrained to \code{R0 > 0}}
#'   \item{scale.sia}{SIA coverage scaling factor, constrained to \code{(0, 1)}}
#'   \item{rho}{RI-SIA dose correlation, in \code{(-1, 1)}}
#' }
#'
#' @details
#' The transformations are:
#' \deqn{R0 = \exp(\theta_1)}
#' \deqn{\text{scale.sia} = \text{plogis}(\theta_2)}
#' \deqn{\rho = \tanh(\theta_3) \quad \text{(only when fix.rho = NA)}}
#'
#' @examples
#' # Fixed rho (2-parameter model)
#' TransformMeaslesSerologyParameters(c(log(12), qlogis(0.7)), fix.rho = 1)
#'
#' # Estimated rho (3-parameter model)
#' TransformMeaslesSerologyParameters(c(log(12), qlogis(0.7), atanh(0.5)), fix.rho = NA)
#'
#' @export
TransformMeaslesSerologyParameters <- function(theta, fix.rho = 1) {

  if (!is.na(fix.rho)) {
    if (length(theta) != 2)
      stop("theta must have length 2 when fix.rho is a number")
    return(list(
      R0        = exp(theta[1]),
      scale.sia = plogis(theta[2]),
      rho       = fix.rho
    ))
  } else {
    if (length(theta) != 3)
      stop("theta must have length 3 when fix.rho = NA (rho is estimated)")
    return(list(
      R0        = exp(theta[1]),
      scale.sia = plogis(theta[2]),
      rho       = tanh(theta[3])
    ))
  }
}
