#' Negative-binomial log-likelihood for age-stratified measles case counts
#'
#' Computes the negative-binomial log-likelihood of observed case counts given
#' model-predicted incidence and a reporting fraction. Expected observed cases
#' are \eqn{E[\text{cases}] = \phi \times \text{pred.incidence}}, with
#' variance \eqn{\mu + \mu^2 / \kappa}.
#'
#' @param casedata data.frame with columns \code{cases} (observed count) and
#'   \code{pred.incidence} (model-predicted new infections, as appended by
#'   \code{\link{GetAgeSpecificAnnualIncidence}} or
#'   \code{\link{GetPredictedMeaslesSerologyAndCases}}).
#' @param phi numeric in \code{(0, 1)}. Case reporting fraction. Scales
#'   predicted incidence to expected observed counts.
#' @param kappa numeric \code{> 0}. Negative-binomial size (dispersion)
#'   parameter. Larger values → lower overdispersion; as
#'   \code{kappa → Inf} the distribution converges to Poisson.
#' @param eps numeric. Small floor applied to expected case counts to keep
#'   the NegBin mean strictly positive (default \code{1e-6}).
#'
#' @return numeric scalar. Total negative-binomial log-likelihood across all
#'   rows of \code{casedata}.
#'
#' @details
#' \deqn{\mu_i = \max(\epsilon,\; \phi \times \text{pred.incidence}_i)}
#' \deqn{\ell = \sum_i \log \text{NegBin}(\text{cases}_i \mid \mu_i, \kappa)}
#'
#' @seealso \code{\link{LogLikMeaslesSerologyAndCases}},
#'   \code{\link{GetAgeSpecificAnnualIncidence}}
#'
#' @export
LogLikMeaslesCases <- function(casedata, phi, kappa, eps = 1e-6) {

  req.cols <- c("cases", "pred.incidence")
  miss.cols <- setdiff(req.cols, names(casedata))
  if (length(miss.cols) > 0)
    stop("casedata is missing required columns: ", paste(miss.cols, collapse = ", "))

  mu <- pmax(eps, phi * casedata$pred.incidence)

  ll <- sum(dnbinom(
    x    = casedata$cases,
    size = kappa,
    mu   = mu,
    log  = TRUE
  ))

  return(ll)
}
