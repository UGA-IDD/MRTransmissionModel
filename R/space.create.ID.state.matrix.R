#' Function creates an ID state matrix with the given number of epidemiological categories, age classes, and subpopulations
#'
#' @param n.age.class the number of age classes
#' @param n.epi.class the number of states in the model
#' @param epi.class.label the names of the states
#' @param n.subpops xxx
#' @param value starting values, defaults to all 0
#'
#' @return an ID.state.matrix object
#' @export
#'

space.create.ID.state.matrix <- function(n.age.class,
                                         n.epi.class = 3,
                                         epi.class.label = c("S","I","R"),
                                         n.subpops = 1,
                                         value = rep(0,n.age.class*n.epi.class*n.subpops)) {

  rc <- new("ID.state.matrix",
            nrow = n.age.class*n.epi.class*n.subpops,
            ncol = 1,
            n.epi.class = n.epi.class,
            epi.class = rep(rep(1:n.epi.class, n.age.class),n.subpops),
            epi.class.label = epi.class.label,
            n.age.class = n.age.class)

  rc[,1] <- value
  return(rc)
}
