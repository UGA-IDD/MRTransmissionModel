#' Aggregate age-specific measles seroprevalence into user-defined age bins
#'
#' This function aggregates model-predicted measles seroprevalence from
#' single-year age groups into broader age bins defined by lower and upper
#' age bounds. Aggregation is performed using population-weighted sums of
#' immune individuals and total population.
#'
#' @param imm.pop numeric matrix. Immune population counts (M + R + V) by
#'   single-year age (rows) and time point (columns).
#' @param pop.age numeric matrix. Total population counts by single-year age
#'   (rows) and time point (columns). Must have the same dimensions as
#'   \code{imm.pop}.
#' @param age.bin.lower integer vector. Lower bounds of age bins (inclusive).
#' @param age.bin.upper integer vector. Upper bounds of age bins (inclusive).
#'
#' @return A list containing:
#' \describe{
#'   \item{imm.pop.bin}{matrix of aggregated immune population counts by age bin
#'   (rows) and time point (columns)}
#'   \item{pop.bin}{matrix of aggregated total population counts by age bin
#'   (rows) and time point (columns)}
#'   \item{seroprev.bin}{matrix of aggregated seroprevalence by age bin
#'   (rows) and time point (columns)}
#' }
#'
#' @details
#' Seroprevalence is calculated as:
#' \deqn{seroprev = \frac{\sum_a (M_a + R_a + V_a)}{\sum_a N_a}}
#'
#' The vectors \code{age.bin.lower} and \code{age.bin.upper} must be of equal
#' length and define valid index ranges corresponding to rows of
#' \code{imm.pop} and \code{pop.age}.
#'
#' If the total population in a bin is zero, seroprevalence is set to 1 by
#' convention.
#'
#' @examples
#' lower <- c(0, 5, 10)
#' upper <- c(4, 9, 14)
#'
#' out <- AggregateMeaslesSeroprevalenceByAgeBins(
#'   imm.pop = imm.pop,
#'   pop.age = pop.age,
#'   age.bin.lower = lower,
#'   age.bin.upper = upper
#' )
#'
#' @export
AggregateMeaslesSeroprevalenceByAgeBins <- function(imm.pop, pop.age,
                                                age.bin.lower, age.bin.upper){

  if(length(age.bin.lower) != length(age.bin.upper)){
    stop("age.bin.lower and age.bin.upper must have the same length")
  }

  n_bins <- length(age.bin.lower)
  n_time <- ncol(pop.age)

  bin.imm <- matrix(NA_real_, nrow = n_bins, ncol = n_time)
  bin.pop <- matrix(NA_real_, nrow = n_bins, ncol = n_time)
  bin.seroprev <- matrix(NA_real_, nrow = n_bins, ncol = n_time)

  for(b in seq_len(n_bins)){
    idx <- age.bin.lower[b]:age.bin.upper[b]

    bin.imm[b, ] <- colSums(imm.pop[idx, , drop = FALSE], na.rm = TRUE)
    bin.pop[b, ] <- colSums(pop.age[idx, , drop = FALSE], na.rm = TRUE)
    bin.seroprev[b, ] <- bin.imm[b, ] / bin.pop[b, ]
  }

  bin.seroprev[bin.pop == 0] <- 1

  rownames(bin.imm) <- paste0(age.bin.lower, "-", age.bin.upper)
  rownames(bin.pop) <- paste0(age.bin.lower, "-", age.bin.upper)
  rownames(bin.seroprev) <- paste0(age.bin.lower, "-", age.bin.upper)

  colnames(bin.imm) <- colnames(imm.pop)
  colnames(bin.pop) <- colnames(pop.age)
  colnames(bin.seroprev) <- colnames(pop.age)

  return(list(
    imm.pop.bin = bin.imm,
    pop.bin = bin.pop,
    seroprev.bin = bin.seroprev
  ))
}
