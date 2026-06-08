#' Negative log-likelihood for joint serology+cases calibration on transformed scale
#'
#' Computes the negative combined log-likelihood (serology + case counts) using
#' unconstrained transformed parameters. Intended for use with \code{optim()}.
#'
#' @param theta numeric vector of transformed parameters. Length depends on
#'   \code{fix.R0} and \code{fix.rho}: see
#'   \code{\link{TransformMeaslesSerologyParameters}} for the full layout.
#'   With defaults (\code{fix.R0 = NA}, \code{fix.rho = 1}): length 4 —
#'   \code{(log(R0), qlogis(scale.sia), qlogis(phi), log(kappa))}.
#'   When \code{fix.R0} is a number: length 3 —
#'   \code{(qlogis(scale.sia), qlogis(phi), log(kappa))}.
#' @param serodata data.frame with columns \code{survey.time.point},
#'   \code{age.bin.lower}, \code{age.bin.upper}, \code{n_tested},
#'   \code{n_positive}.
#' @param casedata data.frame with columns \code{age.lower}, \code{age.upper},
#'   \code{year}, \code{cases}.
#' @param setup country setup object returned by a \code{setupCountry_*} helper.
#' @param year numeric. Simulation start year (default 1980).
#' @param t.max numeric. Number of years to simulate.
#' @param fix.rho numeric scalar or \code{NA}. When a number, \code{rho} is
#'   fixed and \code{theta} has no atanh element. Default \code{1}.
#' @param fix.R0 numeric scalar or \code{NA}. When a number, \code{R0} is
#'   fixed at that value and \code{theta} has no log element. The R0 prior is
#'   not applied when \code{fix.R0} is a number. Default \code{NA}.
#' @param prior.R0.meanlog numeric. Mean of log-normal prior on \code{R0} (log
#'   scale). Ignored when \code{fix.R0} is a number. Default \code{log(14)}.
#' @param prior.R0.sdlog numeric. SD of log-normal prior on \code{R0} (log
#'   scale). Ignored when \code{fix.R0} is a number. Default \code{0.4}.
#' @param age.classes numeric vector. Default \code{c(1:240, seq(252, 1212, 12))}.
#' @param generation.time numeric. Default \code{0.5}.
#' @param seasonal.amp numeric. Default \code{0.15}.
#' @param age0is6to11monly logical. Default \code{FALSE}.
#' @param eps numeric. Default \code{1e-10}.
#'
#' @return numeric scalar. The negative combined log-likelihood (plus negative
#'   log prior on R0 when \code{fix.R0 = NA}).
#'
#' @details
#' An internal iteration counter is maintained across calls within the same R
#' session and printed alongside the current parameter values. Re-source the
#' file to reset the counter.
#'
#' @seealso \code{\link{TransformMeaslesSerologyParameters}},
#'   \code{\link{LogLikMeaslesSerologyAndCases}},
#'   \code{\link{FitMeaslesSerologyAndCases}}
#'
#' @export
NegLogLikMeaslesSerologyAndCases.transformed <- local({

  iter <- 0

  function(theta,
           serodata,
           casedata,
           setup,
           year             = 1980,
           t.max,
           fix.rho          = 1,
           fix.R0           = NA,
           prior.R0.meanlog = log(14),
           prior.R0.sdlog   = 0.4,
           age.classes      = c(1:240, seq(252, 1212, 12)),
           generation.time  = 0.5,
           seasonal.amp     = 0.15,
           age0is6to11monly = FALSE,
           eps              = 1e-10) {

    iter <<- iter + 1

    pars <- TransformMeaslesSerologyParameters(theta,
                                               fix.rho             = fix.rho,
                                               fix.R0              = fix.R0,
                                               include.case.params = TRUE)

    cat(sprintf(
      "Iter %d | R0 = %.3f | scale.sia = %.4f | rho = %.4f | phi = %.4f | kappa = %.4f\n",
      iter, pars$R0, pars$scale.sia, pars$rho, pars$phi, pars$kappa
    ))

    ll <- LogLikMeaslesSerologyAndCases(
      R0               = pars$R0,
      rho              = pars$rho,
      scale.sia        = pars$scale.sia,
      phi              = pars$phi,
      kappa            = pars$kappa,
      serodata         = serodata,
      casedata         = casedata,
      setup            = setup,
      year             = year,
      t.max            = t.max,
      age.classes      = age.classes,
      generation.time  = generation.time,
      seasonal.amp     = seasonal.amp,
      age0is6to11monly = age0is6to11monly,
      eps              = eps
    )

    ll.prior <- if (is.na(fix.R0)) {
      dlnorm(pars$R0, meanlog = prior.R0.meanlog,
             sdlog = prior.R0.sdlog, log = TRUE)
    } else {
      0
    }

    return(-(ll + ll.prior))
  }
})
