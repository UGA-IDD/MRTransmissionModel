#' Fit measles serology calibration model by maximum likelihood
#'
#' Maximizes the binomial likelihood of observed serological data over
#' \code{R0} and \code{rho} using \code{optim()} on the transformed
#' (unconstrained) parameter scale. Optionally runs a coarse grid search first
#' to find a good starting point. A Part 1 result cache is maintained
#' internally across all likelihood evaluations so that the expensive transient
#' burn-in is not repeated when only \code{rho} changes.
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
#'   transformed scale: \code{c(log(R0), atanh(rho))}. Overridden by
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
#' @param rho.range numeric vector of length 2. Lower and upper bounds
#'   for the \code{rho} grid search (default \code{c(-0.9, 0.9)}).
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
#'     transformed scale (\code{log_R0}, \code{atanh_rho})}
#'   \item{par.natural}{named list with fitted \code{R0} and \code{rho}
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
#' \deqn{R0 = \exp(\theta_1), \quad \rho = \tanh(\theta_2)}
#' ensuring \code{R0 > 0} and \code{-1 < rho < 1} throughout.
#'
#' A single \code{.cache} environment is created at the start of each
#' \code{FitMeaslesSerology} call and shared across all likelihood evaluations
#' (grid search and \code{optim()}). The Part 1 transient burn-in result is
#' stored by \code{R0} and reused whenever \code{R0} has not changed since the
#' previous call — most beneficial during grid search (same \code{R0} across
#' all \code{rho} values in a row) and during derivative-free
#' optimization steps that hold \code{R0} roughly fixed.
#'
#' When \code{n.grid > 1} the grid is ordered so that \code{rho} varies
#' fastest, maximising Part 1 cache hits.
#'
#' @seealso \code{\link{NegLogLikMeaslesSerology.transformed}},
#'   \code{\link{TransformMeaslesSerologyParameters}},
#'   \code{\link{GetPredictedMeaslesSerology}}
#'
#' @export
FitMeaslesSerology <- function(
    serodata,
    setup,
    year = 1980,
    t.max,
    par.init = c(log(12), atanh(0)),
    method = "Nelder-Mead",
    hessian = FALSE,
    n.grid = 0,
    R0.range = c(4, 30),
    rho.range = c(-0.9, 0.9),
    age.classes = c(1:240, seq(252, 1212, 12)),
    ...
){

  if(length(par.init) != 2){
    stop("par.init must be a numeric vector of length 2")
  }

  .cache <- new.env(parent = emptyenv())

  if (n.grid > 1) {
    R0.grid  <- exp(seq(log(R0.range[1]), log(R0.range[2]), length.out = n.grid))
    rho.grid <- seq(rho.range[1], rho.range[2], length.out = n.grid)
    # rho varies fastest so consecutive evals share R0 and hit the Part 1 cache
    grid.vals <- expand.grid(rho = rho.grid, R0 = R0.grid)

    cat(sprintf(
      "Grid search: %d points, R0 in [%.1f, %.1f], rho in [%.2f, %.2f]...\n",
      nrow(grid.vals), R0.range[1], R0.range[2], rho.range[1], rho.range[2]
    ))

    grid.negll <- mapply(function(R0, rho) {
      theta <- c(log(R0), atanh(rho))
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
    }, grid.vals$R0, grid.vals$rho)

    best.idx <- which.min(grid.negll)
    par.init <- c(log(grid.vals$R0[best.idx]), atanh(grid.vals$rho[best.idx]))
    cat(sprintf(
      "Grid best: R0=%.2f, rho=%.3f (negLL=%.2f)\n",
      grid.vals$R0[best.idx], grid.vals$rho[best.idx], grid.negll[best.idx]
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
    rho         = par.nat$rho,
    .cache      = .cache,
    age.classes = age.classes,
    ...
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

  out <- list(
    optim           = fit,
    par.transformed = stats::setNames(fit$par, c("log_R0", "atanh_rho")),
    par.natural     = par.nat,
    logLik          = -fit$value,
    predictions     = pred,
    convergence     = fit$convergence,
    message         = if(!is.null(fit$message)) fit$message else NA_character_,
    cov.transformed = cov.transformed
  )

  return(out)
}
