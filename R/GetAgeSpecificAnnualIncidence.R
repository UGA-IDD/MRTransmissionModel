#' Extract age-band-specific annual incidence from an MSIRV simulation result
#'
#' Computes model-predicted annual incidence (new infections) by age band and
#' calendar year by taking the change in the recovered (R) compartment over
#' each year. Accepts pre-aggregated age bands (e.g. 5-year groups) via
#' \code{age.lower} and \code{age.upper} columns, avoiding the zero-incidence
#' problem that arises when many single-year age cells have near-zero predicted
#' infections but small non-zero observed counts (surveillance noise).
#'
#' @param res result matrix from an MSIRV simulation (\code{tmp.res@result}).
#' @param trans transition object (\code{tmp.res@experiment.def@trans}).
#'   Used to identify which age classes fall within each age band.
#' @param epi.state integer vector (\code{tmp.res@experiment.def@trans@epi.class}).
#'   Used to select the R compartment rows (epi.state == 4).
#' @param year.sim numeric. Simulation start year (typically 1980).
#' @param casedata data.frame with columns \code{age.lower} (integer years,
#'   inclusive lower bound of age band), \code{age.upper} (integer years,
#'   inclusive upper bound of age band), \code{year} (calendar year), and
#'   \code{cases} (observed case count). All years must fall within the
#'   simulation range and all ages must be covered by the model age classes.
#'
#' @return \code{casedata} with one column appended:
#' \describe{
#'   \item{pred.incidence}{model-predicted new infections summed across all
#'     model age classes within \code{[age.lower, age.upper]} for that year}
#' }
#'
#' @details
#' Annual incidence for age band \code{[a_low, a_high]} in year \code{y}:
#' \deqn{\Delta R = \sum_{months \in [a_{low}, a_{high}]} R_{month, t_{end}} -
#'       \sum_{months \in [a_{low}, a_{high}]} R_{month, t_{start}}}
#' where \code{t_start = (y - year.sim) * 24 + 1} and
#' \code{t_end = t_start + 23}. Negative values are floored at 0.
#'
#' Pre-aggregate \code{casedata} to 5-year bands before passing to avoid
#' zero-incidence cells driven by surveillance noise in single-year age bins:
#' \preformatted{
#' casedata.5yr <- casedata |>
#'   dplyr::mutate(age.lower = floor(age / 5) * 5,
#'                 age.upper = age.lower + 4) |>
#'   dplyr::group_by(age.lower, age.upper, year) |>
#'   dplyr::summarise(cases = sum(cases), .groups = "drop")
#' }
#'
#' @export
GetAgeSpecificAnnualIncidence <- function(res, trans, epi.state, year.sim, casedata) {

  req.cols <- c("age.lower", "age.upper", "year", "cases")
  miss.cols <- setdiff(req.cols, names(casedata))
  if (length(miss.cols) > 0)
    stop("casedata is missing required columns: ", paste(miss.cols, collapse = ", "))

  if (any(casedata$age.lower > casedata$age.upper))
    stop("All age.lower values must be <= age.upper")

  r.mat       <- res[epi.state == 4, , drop = FALSE]
  age.classes <- trans@age.class

  n              <- nrow(casedata)
  pred.incidence <- numeric(n)

  for (i in seq_len(n)) {
    lower.month <- casedata$age.lower[i] * 12
    upper.month <- (casedata$age.upper[i] + 1) * 12
    age.rows    <- which(age.classes > lower.month & age.classes <= upper.month)

    if (length(age.rows) == 0)
      stop(sprintf(
        "No age classes found for age band %d-%d years — check age.classes covers this range",
        casedata$age.lower[i], casedata$age.upper[i]
      ))

    t.start <- (casedata$year[i] - year.sim) * 24 + 1
    t.end   <- t.start + 23

    pred.incidence[i] <- max(0,
      sum(r.mat[age.rows, t.end]) - sum(r.mat[age.rows, t.start])
    )
  }

  casedata$pred.incidence <- pred.incidence
  return(casedata)
}
