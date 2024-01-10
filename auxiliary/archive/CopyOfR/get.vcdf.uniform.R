#' Function to find the non-scaled vaccination CDF (assuming uniform pdf)
#'
#' @param min.age numeric; age in months
#' @param max.age numeric; age in months
#'
#' @return CDF vaccinated across age in months
#' @export
#'

get.vcdf.uniform <- function(min.age=25, max.age=36){

  ages <- min.age:max.age
  cdf <- cumsum(rep(1/length(ages), length(ages)))

  vcdf <- new("vaccine.cdf.byage",
              cdf=cdf, ages=ages)

  return(vcdf)

}
