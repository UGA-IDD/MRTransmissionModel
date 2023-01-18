#' Make a diagonal WAIFW matrix
#'
#' @param age.class.boundries the upper age limit for each age class in years
#'
#' @return a diagonal WAIFW matrix
#' @export
#'

get.diagonal.WAIFW <- function (age.class.boundries = (1:120/12)) {
  n.age.cats <- length(age.class.boundries)
  ages.to.use <- (age.class.boundries +
                    c(0,age.class.boundries[2:n.age.cats-1]))/2

  foi.matrix <- matrix(0,nrow= n.age.cats,
                       ncol=n.age.cats)
  diag(foi.matrix) <- 1
  colnames(foi.matrix) <- age.class.boundries
  rownames(foi.matrix) <- age.class.boundries

  # fully diagonal perhaps a bit extreme?
  #  here can infect 6months above and below (at least with current age config)
  #foi.matrix[cbind(2:n.age.cats,1:(n.age.cats-1))] <- 1
  #foi.matrix[cbind(1:(n.age.cats-1),2:n.age.cats)] <- 1

  ## Updated:
  # age structure for < 15y is monthly, so let them infect if they are the same year
  # ie: 1-12m with 1-12m, 13-24m with 13-24m, etc.

  for(i in 1:15){
    for(j in 1:12){

      foi.matrix[cbind(rep((12*(i-1)+1):(12*(i-1)+12),12),c(rep((12*(i-1)+j),12)))] <- 1

    }}
  return(foi.matrix)
}
