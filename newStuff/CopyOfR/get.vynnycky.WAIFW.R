#' Make a WAIFW matrix to mimic Emilia's - currently data-free
#'
#' @param age.class.boundries the upper age limit for each age class in years
#' @param beta_young beta on young individuls (<=13 years old)
#' @param beta_old beta on old individuals (>13 years old)
#'
#' @return a WAIFW matrix that mimics Emilia's
#' @export
#'

get.vynnycky.WAIFW <- function (age.class.boundries = (1:120/12), beta_young, beta_old) {
  #get the force of infection for each age class.
  #make sure we use the mid point of each category
  n.age.cats <- length(age.class.boundries)
  ages.to.use <- (age.class.boundries +
                    c(0,age.class.boundries[2:n.age.cats-1]))/2

  #make the unscaled matrix version of the transmission over age
  foi.matrix <- matrix(1,nrow= n.age.cats,
                       ncol=n.age.cats)

  #Emilia's WAIFW has old and young betas, and contact between then is beta_old*0.7
  foi.matrix[ages.to.use<=13,ages.to.use<=13] <- beta_young
  foi.matrix[ages.to.use>13,ages.to.use>13] <- beta_old
  foi.matrix[ages.to.use<=13,ages.to.use>13] <- beta_old*0.7
  foi.matrix[ages.to.use>13,ages.to.use<=13] <- beta_old*0.7

  colnames(foi.matrix) <- age.class.boundries
  rownames(foi.matrix) <- age.class.boundries

  return(foi.matrix)
}
