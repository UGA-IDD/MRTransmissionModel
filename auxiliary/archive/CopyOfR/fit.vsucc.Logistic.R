#' Fits a vsucc.Logistic object to some data using logistic regression
#'
#' @param data a data frame of age (in months), number of success and total tests
#' @param full.efficacy the maximum efficacy we wish to allow
#'
#' @importFrom stats glm binomial
#' @importFrom methods new
#'
#' @return a vsucc.Logistic object
#' @export
#'

fit.vsucc.Logistic <- function(data, full.efficacy) {
  succ.fail <- c()
  age <- c()
  for (i in 1:nrow(data)) {
    succ.fail <- c(succ.fail, rep (1,data$successes[i]))
    succ.fail <- c(succ.fail, rep (0,data$N[i]-data$successes[i]))
    age <- c(age, rep(data$age[i], data$N[i]))
  }

  mdl<-glm(succ.fail~age, family=binomial(link="logit"))

  res <- new("vsucc.Logistic",
             intercept=as.numeric(mdl$coef[1]),
             mo.eff=as.numeric(mdl$coef[2]),
             full.efficacy=as.numeric(full.efficacy))

  return(res)
}
