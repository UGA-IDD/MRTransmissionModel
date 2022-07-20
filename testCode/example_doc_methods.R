##############################################################
#' The title, in this case: Helloworld-ify the argument.
#'
#' Some additional details about this S4 generic and its methods.
#' The extra blank line between this section and the title is
#' critical for roxygen2 to differentiate the title from the
#' description section.
#'
#' @param x Description of \code{x}. The main argument in this
#'  example. Most often has such and such properties.
#'
#' @param y Description of \code{y}. An argument that is rarely
#'  used by \code{"helloworld"} methods.
#'
#' @param ... Additional argument list that might not ever
#'  be used.
#'
#' @return A helloworld-ified argument. Oh, you'll see.
#'
#' @seealso \code{\link{print}} and \code{\link{cat}}
#'
#' @export
#' @docType methods
#' @rdname helloworld-methods
#'
#' @examples
#' helloworld("thisismystring")
#' helloworld(char2helloworld("thisismystring"))
#' helloworld(matrix(0,3,3))
#' helloworld(list(0,0,0))
#' helloworld(integer(0))
setGeneric("helloworld", function(x, y, ...){
  cat("Hello World!")
  cat("\n")
  standardGeneric("helloworld")
})

#' @rdname helloworld-methods
#' @aliases helloworld,ANY,ANY-method
setMethod("helloworld", "ANY", function(x, y, ...){
  cat(class(x))
})

#' @rdname helloworld-methods
#' @aliases helloworld,character,ANY-method
setMethod("helloworld", "character", function(x){
  show(x)
})

#' @rdname helloworld-methods
#' @aliases helloworld,character,character-method
setMethod("helloworld", c("character", "character"), function(x, y){
  show(x)
})
