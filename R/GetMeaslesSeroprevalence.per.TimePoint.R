#' Age-specific measles seroprevalence at requested simulation time points
#'
#' Extracts immune (M + R + V) and total population counts from an MSIRV
#' simulation result matrix, aggregates them from monthly to annual age groups,
#' and returns seroprevalence only at the time points requested. Only
#' \code{time.point.index} columns are processed; the rest of the result matrix
#' is never touched.
#'
#' @param res numeric matrix. Full result matrix from an MSIRV simulation
#'   (\code{n.states x n.timesteps}).
#' @param trans transition object containing \code{age.class} and
#'   \code{epi.class} slots.
#' @param epi.state integer vector of length \code{nrow(res)} giving the
#'   epidemiological state index (1 = M, 2 = S, 3 = I, 4 = R, 5 = V) for
#'   each row of \code{res}.
#' @param no.gens.in.year numeric. Number of model generations per year
#'   (not currently used internally but kept for interface compatibility).
#' @param time.point.index integer vector. Column indices of \code{res} at
#'   which seroprevalence should be extracted (values in \code{1:ncol(res)}).
#'   Only these columns are processed.
#' @param age0is6to11monly logical. If \code{TRUE}, the age-0 row (< 1 year)
#'   is replaced by the sum of months 10-12 only, approximating the window
#'   when maternal antibodies have waned.
#'
#' @return A named list with three matrices, each of dimension
#'   \code{n.age.years x length(time.point.index)}, with columns named by
#'   \code{time.point.index}:
#' \describe{
#'   \item{imm.pop}{immune population (M + R + V) by single-year age group}
#'   \item{pop.age}{total population by single-year age group}
#'   \item{seroprev.age}{seroprevalence (imm.pop / pop.age); set to 1 where
#'     pop.age is zero}
#' }
GetMeaslesSeroprevalence.per.TimePoint <- function(res, trans, epi.state, no.gens.in.year,
                                                   time.point.index, age0is6to11monly = FALSE){

  m <- res[epi.state == 1, , drop = FALSE]
  r <- res[epi.state == 4, , drop = FALSE]
  v <- res[epi.state == 5, , drop = FALSE]

  n.age.yrs <- floor(max(trans@age.class) / 12)

  n_tp <- length(time.point.index)
  m.age   <- matrix(NA_real_, nrow = n.age.yrs, ncol = n_tp)
  r.age   <- matrix(NA_real_, nrow = n.age.yrs, ncol = n_tp)
  v.age   <- matrix(NA_real_, nrow = n.age.yrs, ncol = n_tp)
  pop.age <- matrix(NA_real_, nrow = n.age.yrs, ncol = n_tp)

  for(i in seq_along(time.point.index)){
    t <- time.point.index[i]
    m.age[, i] <- GetNumber.per.AgeYear(vec = m[, t], age.classes = trans@age.class)
    r.age[, i] <- GetNumber.per.AgeYear(vec = r[, t], age.classes = trans@age.class)
    v.age[, i] <- GetNumber.per.AgeYear(vec = v[, t], age.classes = trans@age.class)

    if(age0is6to11monly){
      m.age[1, i] <- sum(m[7:12, t])
      r.age[1, i] <- sum(r[7:12, t])
      v.age[1, i] <- sum(v[7:12, t])
    }

    pop.all <- GetNumber.per.AgeGroup(state = res[, t], trans = trans)
    pop.age[, i] <- GetNumber.per.AgeYear(vec = pop.all, age.classes = trans@age.class)
    if(age0is6to11monly) pop.age[1, i] <- sum(pop.all[7:12])
  }

  imm.pop <- m.age + r.age + v.age

  seroprev.age <- imm.pop / pop.age
  seroprev.age[pop.age == 0] <- 1

  colnames(imm.pop) <- time.point.index
  colnames(pop.age) <- time.point.index
  colnames(seroprev.age) <- time.point.index

  return(list(
    imm.pop = imm.pop,
    pop.age = pop.age,
    seroprev.age = seroprev.age
  ))
}
