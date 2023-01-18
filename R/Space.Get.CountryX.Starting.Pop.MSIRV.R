#' Function to get starting state vector and transition object - built originally from Get.CountryX.Starting.Pop.MSIRV
#'
#' @param uncode UN country code
#' @param generation.time the desired generation time in months
#' @param age.classes vector; the upper limit of the age classes we are interested in months
#' @param maternal.decay.rt rate of maternal decay of immunity to rubella
#' @param exponent exponent for the infected
#' @param frequency.dep boolean
#' @param is.stochastic boolean
#' @param tot.subpop numeric vector; population you want to scale you whole experiment by at the beginning by spatial unit, NULL or numeric
#' @param yr.births.per.1000.bysubpop numeric vector; crude birth rate per year per 1000 by spatial unit
#' @param intro.rate numeric
#' @param flat.WAIFW boolean; if FALSE then default waifw is polymod GB
#' @param space.asdr.object space.nMx object with rate and mid-age and year and space (age specific death rates)
#' @param targeted.intro boolean; FALSE means space introduction out over all age classes, TRUE means to concentrate introductions
#' @param year numeric; year to pull country DFE demography (age structure and population size)
#' @param get.births vector of births with length of 5*n.age.classes
#' @param routine.vac routine vaccination coverage for each subpop
#' @param routine.vac.age.index age index based on age classes for each subpop
#' @param n.subpops the number of subpopulations
#' @param coupling matrix; coupling of subpopulations
#'
#' @return  starting state vector and transition object
#' @export
#'

Space.Get.CountryX.Starting.Pop.MSIRV <- function(uncode,
                                                  generation.time = 0.5,  #generation time in months
                                                  age.classes = c(1:240, seq(241,720,12)),
                                                  maternal.decay.rt=0.95, #based on Metcalf, 2012 paper
                                                  exponent = 0.97,
                                                  frequency.dep=TRUE,
                                                  is.stochastic=FALSE,
                                                  tot.subpop=NULL,
                                                  yr.births.per.1000.bysubpop,
                                                  intro.rate=intro.rate,
                                                  flat.WAIFW=FALSE,
                                                  space.asdr.object=NULL,
                                                  targeted.intro=FALSE,
                                                  year=1980,
                                                  get.births,
                                                  routine.vac=0,
                                                  routine.vac.age.index=12,
                                                  n.subpops = 1,
                                                  coupling){

  ## Calculate the aging rate using the age classes and the generation time
  age.lows <- c(0,age.classes[2:length(age.classes)-1]) # makes it 0 to 697 rather than 1 to 709, so lowest age in wach age range
  ac.sz <- age.classes - age.lows # the size in months of each age range (1 month till age 20, then 12 months to age 59)
  aging.rate <- generation.time/ac.sz # converting generation time units into year/month time units based on size of age class
  aging.rate[length(aging.rate)] <- 0 # forcing the last aging range to 0

  ## Returns the age profile (cols) of survivorship (in units of the generation time) by spatial unit (rows)
  survs <- space.wrapper.surv.prob.over.age(age.classes, generation.time, space.nMx=space.asdr.object, year)

  ## Create WAIFW matrix age.classes X age.classes dimensions, default is POLYMOD based on Great Britain
  waifw <- get.polymod.WAIFW(age.classes/12)
  if (flat.WAIFW) waifw <- get.flat.WAIFW(age.classes/12)

  ## Set up maternal immunity
  maternal.obj = new("maternal.exp.decay", decay.rt=maternal.decay.rt)

  ## Create the transition object
  tran <- space.create.ID.transition.MSIRV(n.age.class = length(age.classes),
                                           aging.rate = aging.rate,
                                           survival.rate = survs,
                                           waifw = waifw,
                                           routine.vac = routine.vac,
                                           routine.vac.age.index = routine.vac.age.index,
                                           maternal.obj = maternal.obj,
                                           time.step = generation.time,
                                           age.class = age.classes,
                                           #birth.rate,
                                           exponent = exponent,
                                           frequency.dep=frequency.dep,
                                           is.stochastic=is.stochastic,
                                           get.births=get.births,
                                           #sia.vac = 0,
                                           #sia.vsucc =  new("vsucc.constant", success.rate=1),
                                           #introduction.rate = 0,
                                           n.subpops = n.subpops, #number of sub-populations
                                           coupling = matrix(1, nrow=n.subpops, ncol=n.subpops))


  ## Putting in starting state where everyone susceptible
  state <- space.create.country.x.DFE.ID.state.matrix(uncode=uncode, tot.subpop=tot.subpop, tran=tran,
                                                      epi.class.label = c("M","S","I","R","V"), year=year)

  ## Update starting births now because need sum(state) = pop for each subpop
  ## the expected number of total births for each time step
  for (s in 1:n.subpops){
    subpop.tmp <- sum(state[tran@s.inds[(tran@n.age.class*s-tran@n.age.class+1):(tran@n.age.class*s)],1])
    tran@birth.rate[s] <- (yr.births.per.1000.bysubpop[s]/1000*subpop.tmp)*generation.time/12
  }

  ## Input transition introduction of infection rate
  if (targeted.intro) {
    tran@introduction.rate <- rep(0, tran@n.age.class*tran@n.subpops)
    tran@introduction.rate[(60:71)*(rep(1:tran@n.subpops, each=length(60:71)))] <- intro.rate
  } else {
    tran@introduction.rate <- rep(intro.rate, tran@n.age.class*tran@n.subpops)
  }

  ## Assumes everyone has maternal protection
  state[tran@m.inds,1] <- state[tran@s.inds,1] * pmaternal(age.lows, maternal.obj)
  state[tran@s.inds,1] <- state[tran@s.inds,1] * (1-pmaternal(age.lows, maternal.obj))

  return(list(state = state, tran = tran))
}
