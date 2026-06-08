#' Transform unconstrained calibration parameters to natural scale
#'
#' Maps unconstrained optimization or MCMC parameters to the natural parameter
#' space used by the measles serology (and optionally case) calibration model.
#'
#' @param theta numeric vector. Elements depend on which parameters are free
#'   (see Details). All free parameters appear in the order: \code{log(R0)},
#'   \code{qlogis(scale.sia)}, \code{atanh(rho)}, \code{qlogis(phi)},
#'   \code{log(kappa)} — each element is omitted when the corresponding
#'   parameter is fixed via \code{fix.R0}, \code{fix.rho}, or
#'   \code{include.case.params}.
#' @param fix.rho numeric scalar or \code{NA}. When a number, \code{rho} is
#'   fixed at that value and \code{theta} needs no atanh element. When
#'   \code{NA}, \code{rho} is a free parameter encoded on the atanh scale.
#'   Default \code{1}.
#' @param fix.R0 numeric scalar or \code{NA}. When a number, \code{R0} is
#'   fixed at that value and \code{theta} needs no log element. When \code{NA}
#'   (default), \code{R0} is a free parameter encoded on the log scale.
#' @param include.case.params logical. If \code{TRUE}, \code{theta} must
#'   also include \code{qlogis(phi)} and \code{log(kappa)} as the last two
#'   elements. Default \code{FALSE}.
#'
#' @return A named list containing:
#' \describe{
#'   \item{R0}{basic reproduction number, \code{> 0}}
#'   \item{scale.sia}{SIA coverage scaling factor, in \code{(0, 1)}}
#'   \item{rho}{RI-SIA dose correlation, in \code{(-1, 1)}}
#'   \item{phi}{case reporting fraction in \code{(0, 1)} (only when
#'     \code{include.case.params = TRUE})}
#'   \item{kappa}{negative-binomial dispersion, \code{>= 0.1} (only when
#'     \code{include.case.params = TRUE}; floored at 0.1 to prevent NegBin
#'     collapse)}
#' }
#'
#' @details
#' Free parameters appear in theta in this fixed order, with fixed parameters
#' simply omitted:
#' \enumerate{
#'   \item \code{log(R0)} — present when \code{fix.R0 = NA}
#'   \item \code{qlogis(scale.sia)} — always present
#'   \item \code{atanh(rho)} — present when \code{fix.rho = NA}
#'   \item \code{qlogis(phi)} — present when \code{include.case.params = TRUE}
#'   \item \code{log(kappa)} — present when \code{include.case.params = TRUE}
#' }
#'
#' @examples
#' # Fixed rho, serology only (2-parameter model)
#' TransformMeaslesSerologyParameters(c(log(12), qlogis(0.7)), fix.rho = 1)
#'
#' # Estimated rho, serology only (3-parameter model)
#' TransformMeaslesSerologyParameters(c(log(12), qlogis(0.7), atanh(0.5)), fix.rho = NA)
#'
#' # Fixed rho, with case params (4-parameter model)
#' TransformMeaslesSerologyParameters(
#'   c(log(12), qlogis(0.7), qlogis(0.1), log(1)),
#'   fix.rho = 1, include.case.params = TRUE)
#'
#' # Fixed R0 and rho, case params only (3-parameter profile fit)
#' TransformMeaslesSerologyParameters(
#'   c(qlogis(0.5), qlogis(0.01), log(1)),
#'   fix.R0 = 14, fix.rho = 1, include.case.params = TRUE)
#'
#' @export
TransformMeaslesSerologyParameters <- function(theta, fix.rho = 1, fix.R0 = NA,
                                               include.case.params = FALSE) {

  n.free <- as.integer(is.na(fix.R0)) + 1L + as.integer(is.na(fix.rho)) +
            if (include.case.params) 2L else 0L

  if (length(theta) != n.free)
    stop(sprintf(
      "theta must have length %d (fix.R0 = %s, fix.rho = %s, include.case.params = %s)",
      n.free,
      if (is.na(fix.R0)) "NA" else fix.R0,
      if (is.na(fix.rho)) "NA" else fix.rho,
      include.case.params
    ))

  idx <- 1L

  if (is.na(fix.R0)) {
    R0  <- exp(theta[idx]); idx <- idx + 1L
  } else {
    R0  <- fix.R0
  }

  scale.sia <- plogis(theta[idx]); idx <- idx + 1L

  if (is.na(fix.rho)) {
    rho <- tanh(theta[idx]); idx <- idx + 1L
  } else {
    rho <- fix.rho
  }

  pars <- list(R0 = R0, scale.sia = scale.sia, rho = rho)

  if (include.case.params) {
    pars$phi   <- plogis(theta[idx])
    pars$kappa <- max(exp(theta[idx + 1L]), 0.1)
  }

  return(pars)
}
