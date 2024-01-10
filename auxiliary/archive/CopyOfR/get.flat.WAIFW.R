#' Make a flat WAIFW matrix
#'
#' @param age.class.boundries the upper age limit for each age class in years
#'
#' @return a flat WAIFW matrix
#' @export
#'

get.flat.WAIFW <- function (age.class.boundries = (1:120/12)) {
  n.age.cats <- length(age.class.boundries)
  ages.to.use <- (age.class.boundries +
                    c(0,age.class.boundries[2:n.age.cats-1]))/2

  foi.matrix <- matrix(1,nrow= n.age.cats,
                       ncol=n.age.cats)

  colnames(foi.matrix) <- age.class.boundries
  rownames(foi.matrix) <- age.class.boundries

  return(foi.matrix)
}
