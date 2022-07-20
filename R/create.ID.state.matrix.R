#' Function creates an ID state matrix with the given number of epidemiologic categories, age classes and labels
#'
#' @param n.age.class the number of age classes
#' @param n.epi.class the number of states in the model
#' @param epi.class.label the names of the states
#' @param value starting values, defauls to all 0
#'
#' @importFrom methods new
#' @return an ID.state.matrix object
#' @export
#'

create.ID.state.matrix <- function(n.age.class,
                                   n.epi.class = 3,
                                   epi.class.label = c("S","I","R"),
                                   value = rep(0,n.age.class*n.epi.class)) {
  rc <- new("ID.state.matrix",
            nrow = n.age.class*n.epi.class,
            ncol = 1,
            n.epi.class = n.epi.class,
            #repeat epi classes for each age
            epi.class = rep(1:n.epi.class, n.age.class),
            epi.class.label = epi.class.label,
            n.age.class = n.age.class
  )

  rc[,1] <- value
  return(rc)
}
