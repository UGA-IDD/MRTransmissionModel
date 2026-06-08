#' Fit measles serology calibration model by maximum likelihood
#'
#' Maximizes the binomial likelihood of observed serological data over
#' \code{R0} and \code{scale.sia} (and optionally \code{rho}) using
#' \code{optim()} on the transformed (unconstrained) parameter scale. When
#' \code{fix.R0} is a number, \code{R0} is held fixed and only
#' \code{scale.sia} (and optionally \code{rho}) are optimized — useful for
#' profile likelihood over an R0 grid.
#'
#' @param serodata data.frame containing observed serological data with columns:
#'   \describe{
#'     \item{survey.time.point}{model time-step index for each observation}
#'     \item{age.bin.lower}{lower bound of age bin in single-year units (inclusive)}
#'     \item{age.bin.upper}{upper bound of age bin in single-year units (inclusive)}
#'     \item{n_tested}{number of individuals tested}
#'     \item{n_positive}{number of seropositive individuals}
#'   }
#' @param setup country setup object returned by a \code{setupCountry_*} helper.
#' @param year numeric. Simulation start year (typically 1980).
#' @param t.max numeric. Number of years to simulate.
#' @param par.init numeric vector. Starting values on the transformed scale.
#'   If \code{NULL} (default), sensible defaults are used. Length and layout
#'   match the free parameters determined by \code{fix.R0} and \code{fix.rho}
#'   (see \code{\link{TransformMeaslesSerologyParameters}}).
#' @param fix.rho numeric scalar or \code{NA}. When a number, \code{rho} is
#'   fixed at that value. When \code{NA}, \code{rho} is also estimated.
#'   Default \code{1}.
#' @param fix.R0 numeric scalar or \code{NA}. When a number, \code{R0} is
#'   fixed at that value and the R0 prior is not applied — useful for profile
#'   likelihood. When \code{NA} (default), \code{R0} is estimated.
#' @param prior.R0.meanlog numeric. Mean of the log-normal prior on \code{R0}
#'   (log scale). Ignored when \code{fix.R0} is a number. Default \code{log(14)}.
#' @param prior.R0.sdlog numeric. SD of the log-normal prior on \code{R0}
#'   (log scale). Ignored when \code{fix.R0} is a number. Default \code{0.4}.
#' @param method character. Optimization method passed to \code{optim()}.
#'   Default \code{"Nelder-Mead"}.
#' @param hessian logical. Should \code{optim()} return the Hessian at the
#'   optimum? Default \code{FALSE}.
#' @param age.classes numeric vector. Upper bounds of age classes in months
#'   (default \code{c(1:240, seq(252, 1212, 12))}).
#' @param generation.time numeric. Generation time in months (default 0.5).
#' @param seasonal.amp numeric. Seasonal forcing amplitude (default 0.15).
#' @param age0is6to11monly logical. If \code{TRUE}, the age-0 seroprevalence
#'   cell uses months 7-12 only (default \code{FALSE}).
#' @param eps numeric. Small value to bound predicted probabilities away from
#'   0 and 1 (default 1e-10).
#'
#' @return A named list containing:
#' \describe{
#'   \item{optim}{raw output object from \code{optim()}}
#'   \item{par.transformed}{named numeric vector of fitted parameters on the
#'     transformed scale}
#'   \item{par.natural}{named list with fitted \code{R0}, \code{scale.sia},
#'     and \code{rho} on the natural scale}
#'   \item{logLik}{maximized log-likelihood (scalar)}
#'   \item{predictions}{\code{serodata} with columns \code{pred.seroprev},
#'     \code{pred.imm.pop}, and \code{pred.pop} appended at the fitted
#'     parameters}
#'   \item{convergence}{integer convergence code from \code{optim()} (0 = success)}
#'   \item{message}{character convergence message, or \code{NA} if none}
#'   \item{cov.transformed}{covariance matrix on the transformed scale (from
#'     Hessian inversion), or \code{NULL} if \code{hessian = FALSE} or the
#'     Hessian is singular}
#' }
#'
#' @seealso \code{\link{NegLogLikMeaslesSerology.transformed}},
#'   \code{\link{TransformMeaslesSerologyParameters}},
#'   \code{\link{GetPredictedMeaslesSerology}}
#'
#' @export
FitMeaslesSerology <- function(
    serodata,
    setup,
    year             = 1980,
    t.max,
    par.init         = NULL,
    fix.rho          = 1,
    fix.R0           = NA,
    prior.R0.meanlog = log(14),
    prior.R0.sdlog   = 0.4,
    method           = "Nelder-Mead",
    hessian          = FALSE,
    age.classes      = c(1:240, seq(252, 1212, 12)),
    generation.time  = 0.5,
    seasonal.amp     = 0.15,
    age0is6to11monly = FALSE,
    eps              = 1e-10
) {

  r0.free  <- is.na(fix.R0)
  rho.free <- is.na(fix.rho)
  n.par    <- as.integer(r0.free) + 1L + as.integer(rho.free)

  if (is.null(par.init)) {
    par.init <- c(
      if (r0.free)  log(12)  else NULL,
                    qlogis(0.5),
      if (rho.free) atanh(0) else NULL
    )
  }

  if (length(par.init) != n.par) {
    stop(sprintf(
      "par.init must have length %d (fix.R0 = %s, fix.rho = %s)",
      n.par,
      if (is.na(fix.R0)) "NA" else fix.R0,
      if (is.na(fix.rho)) "NA" else fix.rho
    ))
  }

  fit <- optim(
    par              = par.init,
    fn               = NegLogLikMeaslesSerology.transformed,
    serodata         = serodata,
    setup            = setup,
    year             = year,
    t.max            = t.max,
    fix.rho          = fix.rho,
    fix.R0           = fix.R0,
    prior.R0.meanlog = prior.R0.meanlog,
    prior.R0.sdlog   = prior.R0.sdlog,
    method           = method,
    hessian          = hessian,
    age.classes      = age.classes,
    generation.time  = generation.time,
    seasonal.amp     = seasonal.amp,
    age0is6to11monly = age0is6to11monly,
    eps              = eps
  )

  par.nat <- TransformMeaslesSerologyParameters(fit$par, fix.rho = fix.rho,
                                                fix.R0 = fix.R0)

  par.names <- c(
    if (r0.free)  "log_R0"          else NULL,
                  "logit_scale.sia",
    if (rho.free) "atanh_rho"       else NULL
  )

  pred <- GetPredictedMeaslesSerology(
    serodata         = serodata,
    setup            = setup,
    year             = year,
    t.max            = t.max,
    R0               = par.nat$R0,
    rho              = par.nat$rho,
    scale.sia        = par.nat$scale.sia,
    age.classes      = age.classes,
    generation.time  = generation.time,
    seasonal.amp     = seasonal.amp,
    age0is6to11monly = age0is6to11monly
  )

  cov.transformed <- if (hessian && !is.null(fit$hessian)) {
    tryCatch(
      solve(fit$hessian),
      error = function(e) {
        warning("Hessian is singular; cov.transformed set to NULL")
        NULL
      }
    )
  } else NULL

  list(
    optim           = fit,
    par.transformed = stats::setNames(fit$par, par.names),
    par.natural     = par.nat,
    logLik          = -fit$value,
    predictions     = pred,
    convergence     = fit$convergence,
    message         = if (!is.null(fit$message)) fit$message else NA_character_,
    cov.transformed = cov.transformed
  )
}
