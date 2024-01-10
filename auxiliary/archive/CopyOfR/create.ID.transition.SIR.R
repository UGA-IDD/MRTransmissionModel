#' Simple method to create an ID.transition.matrix
#'
#' @param n.age.class the number of age classes
#' @param aging.rate the rate which people move to the next age class
#' @param survival.rate the percent who survive at each time step
#' @param waifw whom acquires infection from whom matrix
#' @param age.class the actual age classes defined by upper age in class
#' @param n.epi.class the number of states in the model
#' @param birth.rate xxx
#' @param introduction.rate xxx
#' @param exponent xxx
#' @param frequency.dep xxx
#' @param is.stochastic xxx
#' @param get.births xxx
#'
#' @importFrom methods new
#'
#' @return ID transition matrix
#' @export
#'

create.ID.transition.SIR <- function(n.age.class,
                                     aging.rate,
                                     survival.rate,
                                     waifw,
                                     age.class = numeric(),
                                     n.epi.class = 3,
                                     birth.rate = 0,
                                     introduction.rate = 0,
                                     exponent = 1,
                                     frequency.dep=F,
                                     is.stochastic=F,
                                     get.births = function(pop,tran){return(c(tran@birth.rate,rep(0,length(pop)-1)))}) {
  rc <- new("ID.transition.SIR",
            n.epi.class = n.epi.class,
            epi.class = rep(1:n.epi.class, n.age.class),
            epi.class.label = c("S","I","R"),
            n.age.class = n.age.class,
            age.class = age.class,
            aging.rate = aging.rate,
            survival.rate = survival.rate,
            birth.rate = birth.rate,
            waifw = waifw,
            introduction.rate = introduction.rate,
            exponent = exponent,
            frequency.dep = frequency.dep,
            is.stochastic = is.stochastic,
            get.births = get.births)

  rc@s.inds <- which(rc@epi.class==1)
  rc@i.inds <- which(rc@epi.class==2)
  rc@r.inds <- which(rc@epi.class==3)

  rc@age.surv.matrix <- matrix(0,nrow = n.age.class*n.epi.class,
                               ncol = n.age.class*n.epi.class)


  #for zeroing out all the appropriate transitions
  template.mtrx <- matrix(c(1,1,0,0,1,1,0,0,1), nrow=3, ncol=3)

  #for each age class
  for (i in 1:n.age.class) {
    #first fill in the diagonal matrix
    inds <- (i-1)*n.epi.class+(1:n.epi.class)
    rc@age.surv.matrix[inds,inds] <- (1-aging.rate[i])*survival.rate[i] *
      template.mtrx

    #now fill in the off diagnal matrix
    if (i!=n.age.class) {
      #first fill in the diagonal matrix
      inds2 <- (i)*n.epi.class+(1:n.epi.class)
      rc@age.surv.matrix[inds2,inds] <- aging.rate[i]*survival.rate[i] *
        template.mtrx
    }
  }

  return(rc)
}

