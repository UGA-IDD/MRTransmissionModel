#' Function to get starting state and transition object
#'
#' @param uncode UN country code
#' @param generation.time the desired generation time in months
#' @param age.classes vector; the upper limit of the age classes we are interested in months
#' @param maternal.decay.rt rate of maternal decay of immunity to rubella
#' @param exponent numeric; exponent for the infected
#' @param frequency.dep logical; TRUE
#' @param is.stochastic logical; FALSE
#' @param tot.pop numeric; population you want to scale you whole experiment by at the beginning, NULL or numeric
#' @param yr.births.per.1000 numeric; crude birth rate per year per 1000
#' @param intro.rate xxx
#' @param flat.WAIFW logical; if F then default waifw is polymod GB
#' @param asdr.object nMx object with rate and mid-age (age specific death rates)
#' @param targeted.intro logical; FALSE means space introduction out over all age classes, TRUE means to concentrate introductions
#' @param year numeric; year to pull country DFE demography (age structure and population size)
#' @param get.births xxx
#' @param use_montagu_demog xxx
#' @param routine.vac xxx
#' @param routine.vac.age.index xxx
#'
#' @return starting state and transition object
#' @export
#'

Get.CountryX.Starting.Pop.MSIRV <- function(uncode ,
                                            generation.time = 0.5,  #generation time in months
                                            age.classes = c(1:240, seq(241,720,12)),
                                            maternal.decay.rt=0.95, #based on Metcalf, 2012 paper
                                            exponent = 0.97,
                                            frequency.dep=TRUE,
                                            is.stochastic=FALSE,
                                            tot.pop=NULL,
                                            yr.births.per.1000,
                                            intro.rate=intro.rate,
                                            flat.WAIFW=F,
                                            asdr.object=NULL,
                                            targeted.intro=FALSE,
                                            year=1990,
                                            get.births,
                                            use_montagu_demog=T,
                                            routine.vac=0,
                                            routine.vac.age.index=12){

  ## Calculate the aging rate using the age classes and the generation time
  age.lows <- c(0,age.classes[2:length(age.classes)-1]) # makes it 0 to 697 rather than 1 to 709, so lowest age in wach age range
  ac.sz <- age.classes - age.lows # the size in months of each age range (1 month till age 20, then 12 months to age 59)
  aging.rate <- generation.time/ac.sz # converting generation time units into year/month time units based on size of age class
  aging.rate[length(aging.rate)] <- 0 # forcing the last aging range to 0

  ## Returns the age profile of survivorship in units of the generation time
  survs <- create.surv.prob.over.age(age.classes=age.classes, generation.time=generation.time, nMx=asdr.object, year=year)

  ## Create WAIFW matrix age.classes X age.classes dimensions, default is polymod basedo on great britain
  if (flat.WAIFW) {
    waifw <- get.flat.WAIFW(age.classes/12)
  }else{
    waifw <- get.polymod.WAIFW(age.classes/12)
  }


  ## Set up maternal immunity
  maternal.obj = new("maternal.exp.decay", decay.rt=maternal.decay.rt)

  ## Create the transition object
  tran <- create.ID.transition.MSIRV(n.age.class = length(age.classes),
                                     aging.rate = aging.rate,
                                     survival.rate = survs,
                                     waifw = waifw,
                                     routine.vac=routine.vac,
                                     routine.vac.age.index=routine.vac.age.index,
                                     maternal.obj = maternal.obj,
                                     time.step = generation.time,
                                     age.class = age.classes,
                                     birth.rate =
                                       (yr.births.per.1000*tot.pop/1000)*
                                       (generation.time/12), #the expected number of total births for each time step
                                     exponent = exponent,
                                     frequency.dep=frequency.dep,
                                     is.stochastic=is.stochastic,
                                     get.births=get.births)


  ## Putting in starting state where everyone susceptible
  state <- create.country.x.DFE.ID.state.matrix(uncode=uncode, tot.pop=tot.pop, tran=tran,
                                                epi.class.label = c("M","S","I","R","V"), year=year,
                                                use_montagu_demog=use_montagu_demog)

  ## Update births now because need sum(state) = pop
  ## the expected number of total births for each time step
  tran@birth.rate <- (yr.births.per.1000*sum(state)/1000)*generation.time/12

  ## Input transition introduction of infection rate
  if (targeted.intro) {
    tran@introduction.rate <- rep(0, tran@n.age.class)
    tran@introduction.rate[60:71] <- intro.rate
  } else {
    tran@introduction.rate <- rep(intro.rate, tran@n.age.class)
  }

  ## Assumes everyone has maternal protection
  state[tran@m.inds,1] <- state[tran@s.inds,1] * pmaternal(age.lows, maternal.obj)
  state[tran@s.inds,1] <- state[tran@s.inds,1] * (1-pmaternal(age.lows, maternal.obj))

  return(list(state = state, tran = tran))
}
