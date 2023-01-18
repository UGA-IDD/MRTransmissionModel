#' Convert the results to proportions of each age class in each compartment
#'
#' @param ob object of type sim.results.MSIRV
#' @param ... ignored
#'
#' @include setClasses.R
#' @return result object of proportions instead of absolute numbers
#'
#' @export
#' @docType methods
#' @rdname toProportions-methods

setGeneric("toProportions", function(ob, ...) standardGeneric("toProportions"))

#' @rdname toProportions-methods
#' @aliases toProportions,sim.results.MSIRV-method
setMethod("toProportions",
          "sim.results.MSIRV",
          function (ob) {
            ac <- c()
            for (i in ob@age.class) {
              ac <- c(ac, rep(i,5))
            }
            #print(cbind(round(ob@t,2),t(round(ob[ac==50,],2)),colSums(ob[ac==50,])))
            for (i in ob@age.class) {
              div <- t(matrix(rep(colSums(ob[ac==i,]),5), ncol=5))
              ob[ac==i,] <- ob[ac==i,]/div
            }
            #print(cbind(round(ob@t,2),t(round(ob[ac==50,],2)),colSums(ob[ac==50,])))
            return(ob)
          })
