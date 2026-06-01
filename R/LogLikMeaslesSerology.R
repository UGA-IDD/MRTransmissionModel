#' Binomial log-likelihood for measles serology calibration
#'
#' Runs a full MSIRV simulation for the given parameters and computes the
#' binomial log-likelihood of the observed serological data.
#'
#' @param R0 numeric. Basic reproduction number used to scale the WAIFW matrix.
#' @param sia.scale numeric. Multiplicative adjustment applied to SIA coverage;
#'   clamped to \code{[0, 1]}.
#' @param serodata data.frame with columns \code{survey.time.point},
#'   \code{age.bin.lower}, \code{age.bin.upper}, \code{n_tested},
#'   \code{n_positive}.
#' @param setup country setup object returned by a \code{setupCountry_*} helper.
#' @param year numeric. Simulation start year (default 1980).
#' @param t.max numeric. Number of years to simulate.
#' @param age.classes numeric vector. Upper bounds of age classes in months.
#' @param generation.time numeric. Generation time in months (default 0.5).
#' @param seasonal.amp numeric. Seasonal forcing amplitude (default 0.15).
#' @param age0is9to12monly logical. Passed as \code{age0is9to11monly} to
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
    sia.scale,
    serodata,
    setup,
    year = 1980,
    t.max,
    age.classes = c(1:240, seq(252, 1212, 12)),
    generation.time = 0.5,
    seasonal.amp = 0.15,
    age0is9to12monly = FALSE,
    eps = 1e-10,
    .cache = NULL){

  pred <- GetPredictedMeaslesSerology(
    serodata = serodata,
    setup = setup,
    year = year,
    t.max = t.max,
    R0 = R0,
    sia.scale = sia.scale,
    age.classes = age.classes,
    generation.time = generation.time,
    seasonal.amp = seasonal.amp,
    age0is9to11monly = age0is9to12monly,
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



#' Transform unconstrained calibration parameters to natural scale
#'
#' This helper maps unconstrained optimization or MCMC parameters to the
#' natural parameter space used by the measles serology calibration model.
#'
#' @param theta numeric vector of transformed parameters:
#'   \describe{
#'     \item{theta[1]}{log-scale parameter for \code{R0}}
#'     \item{theta[2]}{logit-scale parameter for \code{sia.scale}}
#'   }
#'
#' @return A named list containing:
#' \describe{
#'   \item{R0}{basic reproduction number, constrained to \code{R0 > 0}}
#'   \item{sia.scale}{SIA coverage adjustment factor, constrained to
#'   \code{0 < sia.scale < 1}}
#' }
#'
#' @details
#' The transformations are:
#' \deqn{R0 = \exp(\theta_1)}
#' \deqn{sia.scale = \frac{1}{1+\exp(-\theta_2)}}
#'
#' @examples
#' TransformMeaslesSerologyParameters(c(log(12), qlogis(0.8)))
#'
#' @export
TransformMeaslesSerologyParameters <- function(theta){

  if(length(theta) != 2){
    stop("theta must be a numeric vector of length 2")
  }

  return(list(
    R0 = exp(theta[1]),
    sia.scale = plogis(theta[2])
  ))
}



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
#'     \item{theta[2]}{logit-scale parameter for \code{sia.scale}}
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
#' @examples
#' LogLikMeaslesSerology.transformed(
#'   theta = c(log(12), qlogis(0.8)),
#'   serodata = serodata,
#'   setup = setup,
#'   year = 1980,
#'   t.max = t.max
#' )
#'
#' @export
LogLikMeaslesSerology.transformed <- function(theta, serodata, setup, year = 1980, t.max, ...){

  par.nat <- TransformMeaslesSerologyParameters(theta = theta)

  ll <- LogLikMeaslesSerology(
    R0 = par.nat$R0,
    sia.scale = par.nat$sia.scale,
    serodata = serodata,
    setup = setup,
    year = year,
    t.max = t.max,
    ...
  )

  return(ll)
}


#' Negative log-likelihood for measles serology calibration on transformed scale
#'
#' This function computes the negative binomial log-likelihood of observed
#' measles serological data using unconstrained transformed parameters. It is
#' intended for use with optimization routines such as \code{optim()}.
#'
#' @param theta numeric vector of transformed parameters:
#'   \describe{
#'     \item{theta[1]}{log-scale parameter for \code{R0}}
#'     \item{theta[2]}{logit-scale parameter for \code{sia.scale}}
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
#'   \item \code{0 < sia.scale < 1}
#' }
#'
#' @seealso \code{\link{TransformMeaslesSerologyParameters}},
#'   \code{\link{LogLikMeaslesSerology.transformed}},
#'   \code{\link{LogLikMeaslesSerology}}
#'
#' @examples
#' fit <- optim(
#'   par = c(log(12), qlogis(0.8)),
#'   fn = NegLogLikMeaslesSerology.transformed,
#'   serodata = serodata,
#'   setup = setup,
#'   year = 1980,
#'   t.max = t.max,
#'   method = "Nelder-Mead"
#' )
#'
#' @export
NegLogLikMeaslesSerology.transformed <- local({

  iter <- 0  # persists across calls

  function(theta, serodata, setup, year = 1980, t.max, ...){

    iter <<- iter + 1

    R0 <- exp(theta[1])
    sia.scale <- plogis(theta[2])

    cat(sprintf("Iter %d | R0 = %.3f | sia.scale = %.3f\n",
                iter, R0, sia.scale))

    ll <- LogLikMeaslesSerology(
      R0 = R0,
      sia.scale = sia.scale,
      serodata = serodata,
      setup = setup,
      year = year,
      t.max = t.max,
      ...
    )

    return(-ll)
  }
})

#' Fit measles serology calibration model by maximum likelihood
#'
#' Maximizes the binomial likelihood of observed serological data over
#' \code{R0} and \code{sia.scale} using \code{optim()} on the transformed
#' (unconstrained) parameter scale. Optionally runs a coarse grid search first
#' to find a good starting point. A Part 1 result cache is maintained
#' internally across all likelihood evaluations so that the expensive transient
#' burn-in is not repeated when only \code{sia.scale} changes.
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
#' @param par.init numeric vector of length 2. Starting values on the
#'   transformed scale: \code{c(log(R0), qlogis(sia.scale))}. Overridden by
#'   the grid search best point when \code{n.grid > 1}.
#' @param method character. Optimization method passed to \code{optim()}.
#'   Default is \code{"Nelder-Mead"} (derivative-free; one evaluation per
#'   step). Use \code{"BFGS"} for faster convergence when the surface is
#'   smooth; with the Part 1 cache active, the extra gradient evaluations are
#'   cheaper than without caching.
#' @param hessian logical. Should \code{optim()} return the Hessian at the
#'   optimum? Default \code{FALSE}.
#' @param n.grid integer. Number of grid points along each parameter axis for
#'   an initial coarse grid search. A value of \code{n.grid = k} evaluates
#'   \code{k^2} points. Set to \code{0} (default) to skip the grid search and
#'   start \code{optim()} directly from \code{par.init}. Values of 4-6 are
#'   usually sufficient to identify the basin of attraction.
#' @param R0.range numeric vector of length 2. Lower and upper bounds for the
#'   \code{R0} grid search on the natural scale (default \code{c(4, 30)}).
#'   Ignored when \code{n.grid <= 1}.
#' @param sia.scale.range numeric vector of length 2. Lower and upper bounds
#'   for the \code{sia.scale} grid search (default \code{c(0.2, 1.0)}).
#'   Ignored when \code{n.grid <= 1}.
#' @param age.classes numeric vector. Upper bounds of age classes in months
#'   (default \code{c(1:240, seq(252, 1212, 12))}). Pass a coarser vector
#'   (e.g. \code{c(1:60, seq(72, 1212, 12))}) to speed up the simulation.
#' @param ... additional arguments passed through to
#'   \code{NegLogLikMeaslesSerology.transformed()} and downstream functions
#'   (e.g. \code{generation.time}, \code{seasonal.amp},
#'   \code{age0is9to12monly}).
#'
#' @return A named list containing:
#' \describe{
#'   \item{optim}{raw output object from \code{optim()}}
#'   \item{par.transformed}{named numeric vector of fitted parameters on the
#'     transformed scale (\code{log_R0}, \code{logit_sia.scale})}
#'   \item{par.natural}{named list with fitted \code{R0} and \code{sia.scale}
#'     on the natural scale}
#'   \item{logLik}{maximized log-likelihood (scalar)}
#'   \item{predictions}{\code{serodata} with columns \code{pred.seroprev},
#'     \code{pred.imm.pop}, and \code{pred.pop} appended at the fitted
#'     parameters}
#'   \item{convergence}{integer convergence code from \code{optim()} (0 = success)}
#'   \item{message}{character convergence message from \code{optim()}, or
#'     \code{NA} if none}
#' }
#'
#' @details
#' Parameters are optimized on the transformed scale:
#' \deqn{R0 = \exp(\theta_1), \quad sia.scale = \mathrm{logistic}(\theta_2)}
#' ensuring \code{R0 > 0} and \code{0 < sia.scale < 1} throughout.
#'
#' A single \code{.cache} environment is created at the start of each
#' \code{FitMeaslesSerology} call and shared across all likelihood evaluations
#' (grid search and \code{optim()}). The Part 1 transient burn-in result is
#' stored by \code{R0} and reused whenever \code{R0} has not changed since the
#' previous call — most beneficial during grid search (same \code{R0} across
#' all \code{sia.scale} values in a row) and during derivative-free
#' optimization steps that hold \code{R0} roughly fixed.
#'
#' When \code{n.grid > 1} the grid is ordered so that \code{sia.scale} varies
#' fastest, maximising Part 1 cache hits.
#'
#' @seealso \code{\link{NegLogLikMeaslesSerology.transformed}},
#'   \code{\link{TransformMeaslesSerologyParameters}},
#'   \code{\link{GetPredictedMeaslesSerology}}
#'
#' @examples
#' # Basic fit from a manual starting point
#' fit <- FitMeaslesSerology(
#'   serodata = serodata,
#'   setup = setup,
#'   year = 1980,
#'   t.max = t.max,
#'   par.init = c(log(16), qlogis(0.8))
#' )
#' fit$par.natural
#' head(fit$predictions)
#'
#' # With a 5x5 grid search to find a good starting point automatically
#' fit <- FitMeaslesSerology(
#'   serodata = serodata,
#'   setup = setup,
#'   year = 1980,
#'   t.max = t.max,
#'   n.grid = 5,
#'   R0.range = c(5, 25),
#'   sia.scale.range = c(0.3, 1.0)
#' )
#'
#' @export
FitMeaslesSerology <- function(
    serodata,
    setup,
    year = 1980,
    t.max,
    par.init = c(log(12), qlogis(0.8)),
    method = "Nelder-Mead",
    hessian = FALSE,
    n.grid = 0,
    R0.range = c(4, 30),
    sia.scale.range = c(0.2, 1.0),
    age.classes = c(1:240, seq(252, 1212, 12)),
    ...
){

  if(length(par.init) != 2){
    stop("par.init must be a numeric vector of length 2")
  }

  .cache <- new.env(parent = emptyenv())

  if (n.grid > 1) {
    R0.grid  <- exp(seq(log(R0.range[1]),       log(R0.range[2]),       length.out = n.grid))
    sia.grid <- seq(sia.scale.range[1], sia.scale.range[2], length.out = n.grid)
    # sia varies fastest so consecutive evals share R0 and hit the Part 1 cache
    grid.vals <- expand.grid(sia = sia.grid, R0 = R0.grid)

    cat(sprintf(
      "Grid search: %d points, R0 in [%.1f, %.1f], sia.scale in [%.2f, %.2f]...\n",
      nrow(grid.vals), R0.range[1], R0.range[2], sia.scale.range[1], sia.scale.range[2]
    ))

    grid.negll <- mapply(function(R0, sia) {
      theta <- c(log(R0), qlogis(sia))
      tryCatch(
        NegLogLikMeaslesSerology.transformed(
          theta       = theta,
          serodata    = serodata,
          setup       = setup,
          year        = year,
          t.max       = t.max,
          .cache      = .cache,
          age.classes = age.classes,
          ...
        ),
        error = function(e) Inf
      )
    }, grid.vals$R0, grid.vals$sia)

    best.idx <- which.min(grid.negll)
    par.init <- c(log(grid.vals$R0[best.idx]), qlogis(grid.vals$sia[best.idx]))
    cat(sprintf(
      "Grid best: R0=%.2f, sia.scale=%.3f (negLL=%.2f)\n",
      grid.vals$R0[best.idx], grid.vals$sia[best.idx], grid.negll[best.idx]
    ))
  }

  fit <- optim(
    par         = par.init,
    fn          = NegLogLikMeaslesSerology.transformed,
    serodata    = serodata,
    setup       = setup,
    year        = year,
    t.max       = t.max,
    method      = method,
    hessian     = hessian,
    .cache      = .cache,
    age.classes = age.classes,
    ...
  )

  par.nat <- TransformMeaslesSerologyParameters(fit$par)

  pred <- GetPredictedMeaslesSerology(
    serodata    = serodata,
    setup       = setup,
    year        = year,
    t.max       = t.max,
    R0          = par.nat$R0,
    sia.scale   = par.nat$sia.scale,
    .cache      = .cache,
    age.classes = age.classes,
    ...
  )

  out <- list(
    optim           = fit,
    par.transformed = stats::setNames(fit$par, c("log_R0", "logit_sia.scale")),
    par.natural     = par.nat,
    logLik          = -fit$value,
    predictions     = pred,
    convergence     = fit$convergence,
    message         = if(!is.null(fit$message)) fit$message else NA_character_
  )

  return(out)
}
