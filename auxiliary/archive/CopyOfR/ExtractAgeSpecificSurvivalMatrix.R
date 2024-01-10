#' Function to produce age survival matrix
#'
#' @param tmp.trans transitions object of class ID.transition.MSIRV
#' @param maternal.obj maternal object of class maternal.exp.decay
#' @param surv.at.timestep.t the estimated age-specific survival rates at time t
#'
#' @return the age survival matrix per time t
#' @export
#'

ExtractAgeSpecificSurvivalMatrix <- function(tmp.trans, maternal.obj, surv.at.timestep.t) {

  age.surv.matrix <- matrix(0,nrow = tmp.trans@n.age.class*tmp.trans@n.epi.class,
                            ncol = tmp.trans@n.age.class*tmp.trans@n.epi.class)


  #for zeroing out all the appropriate transitions from MSIRV to MSIRV
  template.mtrx <- matrix(c(1,1,0,0,1,0,1,1,0,1,0,0,1,1,0,0,0,0,1,0,0,0,0,0,1), nrow=5, ncol=5)

  #useful low age
  low.age <-c(0, tmp.trans@age.class[2:length(tmp.trans@age.class)-1])

  #for each age class
  for (i in 1:tmp.trans@n.age.class) {
    #print(i)
    #first fill in the diagonal matrix
    tmp.template <- template.mtrx
    tmp.template[2,1] <- 0
    inds <- (i-1)*tmp.trans@n.epi.class+(1:tmp.trans@n.epi.class)
    age.surv.matrix[inds,inds] <- (1-tmp.trans@aging.rate[i])*surv.at.timestep.t[i] * tmp.template

    #now fill in the off diagnal matrix
    if (i!=tmp.trans@n.age.class) {
      #create a modified template matrix to move
      #people from the M class to the S class
      tmp.template <- template.mtrx

      #calculate the probability of losing maternal protection during
      p.mat.loss <- (pmaternal(low.age[i], maternal.obj)-
                       pmaternal(tmp.trans@age.class[i], maternal.obj))/pmaternal(low.age[i], maternal.obj)

      p.mat.loss[pmaternal(low.age[i], maternal.obj)==0] <- 1

      tmp.template[1,1] <- 1-p.mat.loss
      tmp.template[2,1] <- 1-tmp.template[1,1]

      inds2 <- (i)*tmp.trans@n.epi.class+(1:tmp.trans@n.epi.class)
      age.surv.matrix[inds2,inds] <- tmp.trans@aging.rate[i]*surv.at.timestep.t[i] * tmp.template
    }
  }
  return(age.surv.matrix)
}
