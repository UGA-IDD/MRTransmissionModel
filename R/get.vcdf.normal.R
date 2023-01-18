#' Function to find the non-scaled vaccination CDF (assuming normal pdf)
#'
#' @param min.age numeric; age in months
#' @param max.age numeric; age in months
#'
#' @importFrom stats pnorm
#' @return CDF vaccinated across age in months
#' @export
#'

get.vcdf.normal <- function(min.age=25, max.age=36){

  ages <- min.age:max.age
  age.range <- length(ages)
  cdf <- c(pnorm(seq(-2,2,4/age.range),0,1)[-1])

  vcdf <- new("vaccine.cdf.byage",
              cdf=cdf, ages=ages)

  return(vcdf)

}
