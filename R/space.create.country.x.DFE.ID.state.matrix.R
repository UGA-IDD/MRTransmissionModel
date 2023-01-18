#' Function to create a state matrix at the disease free equilibrium (DFE) from a transition object
#'
#' @param uncode numeric; UN country code
#' @param tot.subpop vector; sub-populations by which you can rescale if you want
#' @param tran the transition object
#' @param epi.class.label epi class label
#' @param year numeric; options are 1950-2100
#'
#' @importFrom stats smooth.spline predict
#'
#' @return a starting DFE vector of S and R individuals
#' @export
#'

space.create.country.x.DFE.ID.state.matrix <- function(uncode, tot.subpop=NULL, tran,
                                                       epi.class.label = c("M","S","I","R","V"),
                                                       year=1990){

  #make template
  rc <- space.create.ID.state.matrix(tran@n.age.class,
                                     tran@n.epi.class,
                                     epi.class.label = epi.class.label,
                                     tran@n.subpops)

  #fill in the template
  demog <- space.getDemography(uncode=uncode)
  demog.n.ages <- nrow(demog$pop.age.byageclasses.1950.2100)/tran@n.subpops

  #loop through the number of sub-populations
  for (s in 1:tran@n.subpops){

    pop.struct <- demog$pop.age.byageclasses.1950.2100[(demog.n.ages*s-demog.n.ages+1):(demog.n.ages*s),(year-1950+1)]
    age <- as.numeric(names(pop.struct))

    if (!is.null(tot.subpop)) {
      tot.pop <- tot.subpop[s]
    } else {
      tot.pop <- demog$pop.total.1950.2100[s,(year-1950+1)]
    }

    #Turn this into prop of desired age classes
    #note that NUMBER is imposed externally (since we might not want full countries... because interested in stochasticity)
    sp <- smooth.spline(age, pop.struct)  # gives "slope" between age (range 0-100) and population size per age (therefore 101 items)
    pred <- predict(sp, tran@age.class/12)$y # predicts population size per age class (total 280 age classes)
    #plot(age,pop.struct, ylim=c(0,3000000))
    #points(tran@age.class/12,pred,type="l",col=2)

    #Adjust for varying bin width
    pred <- pred*diff(c(0,tran@age.class))
    #added a 0 at the beginning of age classes
    #took the difference between each class n from n-1.   therefore difference is 1 (or 1 month) for item 1-241 and then 12 for items 242-281
    #then weighing each predicted population size per age class by the size of the age class

    #Fill in
    rc[tran@s.inds[(tran@n.age.class*s-tran@n.age.class+1):(tran@n.age.class*s)],1] <- tot.pop*pred/sum(pred)
    #if pop=NULL then assuming that the population w/in each state with missing ages are distributed within the population based on the pop structure that we do have ages for
    #filling in the susceptibles per age class

  }

  return(rc)

}
