#' Function to create a state matrix at the disease free equilibrium (DFE) from a transition object
#'
#' @param uncode UN country code
#' @param tot.pop population by which you can rescale if you want
#' @param tran the transition object
#' @param epi.class.label xxx ask Amy
#' @param year numeric; options are 1950 to 2100
#' @param use_montagu_demog logical; if TRUE use montagu demography data
#' @param demog_data optional demography data
#'
#' @importFrom stats smooth.spline predict
#'
#' @return a starting DFE vector of S and R individuals
#' @export
#'

create.country.x.DFE.ID.state.matrix <- function(uncode, tot.pop=NULL, tran,
                                                 epi.class.label = c("M","S","I","R","V"),
                                                 year=1990, use_montagu_demog=FALSE,
                                                 demog_data = NULL){

  # make template
  rc <- create.ID.state.matrix(tran@n.age.class,
                               tran@n.epi.class,  epi.class.label = epi.class.label)

  if (use_montagu_demog){
    demog <- demog_data

    pop.struct <- demog$pop.age.byageclasses.1950.2100[,(year-1950+1)]
    age <- as.numeric(rownames(demog$pop.age.byageclasses.1950.2100))
    if (is.null(tot.pop)) {
      tot.pop <- demog$pop.total.1950.2100[(year-1950+1)]
    }

  }else{
    demog <- getDemography.wpp2022.mini(uncode=uncode)

    #wpp2019 data is population/1000 so we need to multiply to get true population
    pop.struct <- demog$pop.age.byageclasses.1950.2100[,(year-1950+1)]*1000
    age <- as.numeric(rownames(demog$pop.age.byageclasses.1950.2100))
    if (is.null(tot.pop)) {
      tot.pop <- demog$pop.total.1950.2100[(year-1950+1)]*1000
    }

  }


  #Turn this into prop of desired age classes
  #note that NUMBER is imposed externally (since we might not want full countries... because interested in stochasticity)
  sp <- smooth.spline(age, pop.struct)  # gives "slope" between age (range 0-100) and popualtion size per age (therefore 101 items)
  pred <- predict(sp,tran@age.class/12)$y # predicts population size per age class (total 280 age classes)
  #plot(age,pop.struct, ylim=c(0,3000000))
  #points(tran@age.class/12,pred,type="l",col=2)

  #changed in Dec 2013 - adjust for varying bin width!
  pred <- pred*diff(c(0,tran@age.class))
  #added a 0 at the beginning of age classes - so now 281
  #took the difference between each class n from n-1.   therefore difference is 1 (or 1 month) for item 1-241 and then 12 for items 242-281
  #then weighing each predicted population size per age class by the size of the age class

  #Fill in
  rc[tran@s.inds[1:tran@n.age.class],1] <- tot.pop*pred/sum(pred)
  #if pop=NULL then assuming that the population w/in each state with missing ages are distributed within the popualtion based on the pop structure that we do have ages for
  #filling in the susceptibles per age class

  return(rc)

}
