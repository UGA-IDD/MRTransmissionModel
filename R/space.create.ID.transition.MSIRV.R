#' Method to create an ID.transition.MSIRV.space object
#'
#' @param n.age.class number of age classes (not location specific)
#' @param aging.rate aging rates (includes all locations concatenated)
#' @param survival.rate vector with locations survival rates concatenated
#' @param waifw assume same across locations
#' @param routine.vac xxx
#' @param routine.vac.age.index xxx
#' @param maternal.obj maternal antibody object
#' @param time.step time step in months == generation time
#' @param age.class xxx
#' @param birth.rate.subpop xxx
#' @param sia.vac xxx
#' @param sia.vsucc  xxx
#' @param introduction.rate xxx
#' @param exponent xxx
#' @param frequency.dep xxx
#' @param is.stochastic xxx
#' @param get.births xxx
#' @param n.subpops number of sub-populations
#' @param subpop.class.label xxx
#' @param coupling matrix describing movement between sub-populations
#'
#' @return xxx
#' @export

space.create.ID.transition.MSIRV <- function(n.age.class, #number of age classes (not location specific)
                                             aging.rate,  #aging rates (includes all locations concatenated)
                                             survival.rate, #vector with locations survival rates concatenated
                                             waifw, # assume same across locations
                                             routine.vac = 0,
                                             routine.vac.age.index = 12,
                                             maternal.obj, #maternal antibody object
                                             time.step = 1/2, #time step in months == generation time
                                             age.class = numeric(),
                                             birth.rate.subpop = rep(0, n.subpops),
                                             sia.vac = 0,
                                             sia.vsucc =  new("vsucc.constant", success.rate=1),
                                             introduction.rate = 0,
                                             exponent = exponent,
                                             frequency.dep = F,
                                             is.stochastic = F,
                                             get.births = function(pop,tran){return(c(tran@birth.rate,rep(0,length(pop)-1)))},
                                             n.subpops = 1, #number of sub-populations
                                             subpop.class.label,
                                             coupling = 1 #matrix describing movement between sub-populations
) {


  #First do all of the standard logic
  rc <- new("ID.transition.MSIRV.space",
            n.age.class = n.age.class,
            n.epi.class = 5,
            epi.class.label = c("M","S","I","R","V"),
            epi.class = rep(rep(1:5, n.age.class),n.subpops),
            age.class = age.class,
            aging.rate = aging.rate,
            survival.rate = survival.rate,
            birth.rate = birth.rate.subpop,
            waifw = waifw,
            sia.vsucc = sia.vsucc,
            introduction.rate = introduction.rate,
            exponent = exponent,
            frequency.dep=frequency.dep,
            is.stochastic=is.stochastic,
            get.births = get.births,
            n.subpops=n.subpops,
            subpop.class.label=rep(1:n.subpops,each=n.age.class*5),
            coupling=coupling)

  rc@m.inds = which(rc@epi.class==1)
  rc@s.inds = which(rc@epi.class==2)
  rc@i.inds = which(rc@epi.class==3)
  rc@r.inds = which(rc@epi.class==4)
  rc@v.inds = which(rc@epi.class==5)

  rc@age.surv.matrix <- matrix(0,nrow = n.age.class*rc@n.epi.class*n.subpops,
                               ncol = n.age.class*rc@n.epi.class)

  #for zeroing out all the appropriate transitions
  template.mtrx <- matrix(c(1,1,0,0,1,0,1,1,0,1,0,0,1,1,0,0,0,0,1,0,0,0,0,0,1), nrow=5, ncol=5)

  #useful low age
  low.age <-rep(c(0, age.class[2:length(age.class)-1]), n.subpops)

  #for each subpop and age class
  for (j in 1) { # (j in 1:n.subpops) #xxxamy tmp hack
    for (i in 1:(n.age.class)) {
      #first fill in the diagonal matrix
      tmp.template <- template.mtrx
      tmp.template[2,1] <- 0
      #move into right spatial location
      adj <- (j-1)*n.age.class*rc@n.epi.class
      #index ages
      inds <- adj+(i-1)*rc@n.epi.class+(1:rc@n.epi.class)
      rc@age.surv.matrix[inds,inds-adj] <- (1-aging.rate[i])*survival.rate[j,i]*tmp.template

      #now fill in the off diagonal matrix
      if ((i%%n.age.class)!=0) {
        #create a modified template matrix to move
        #people from the M class to the S class
        tmp.template <- template.mtrx

        #calculate the probability of losing maternal protection duringan age class
        p.mat.loss <- (pmaternal(low.age[i], maternal.obj)-
                         pmaternal(age.class[i], maternal.obj))/pmaternal(low.age[i],maternal.obj)
        p.mat.loss[pmaternal(low.age[i],maternal.obj)==0] <- 1

        tmp.template[1,1] <- 1-p.mat.loss
        tmp.template[2,1] <- 1-tmp.template[1,1]

        inds2 <- adj + (i)*rc@n.epi.class+(1:rc@n.epi.class)
        #print(inds2)
        rc@age.surv.matrix[inds2,inds-adj] <- aging.rate[i]*survival.rate[j,i]*tmp.template
      }
    }
  }
  #xxxamy tmp hack
  for (j in 1:n.subpops){
    rc@age.surv.matrix[(1500*j-1499):(1500*j),1:1500] <- rc@age.surv.matrix[1:1500,1:1500]
  }


  #Next do all of the vaccination logic, concatenating two sites together, with desired coverage
  #get the %vaccinated in each age group per time step time step must be in years!
  vac.per <- c()
  for (l in 1:n.subpops) {
    vac.per.loc <- new("vacc.per.time.step.by.age",
                       pvacc.in.age.class = numeric(length=n.age.class))
    vac.per.loc@pvacc.in.age.class[routine.vac.age.index[l]] <- routine.vac[l]
    vac.per <- c(vac.per, vac.per.loc@pvacc.in.age.class)
  }

  #rc@vac.per <- vac.per #non-patch uses this
  rc@vac.per@pvacc.in.age.class <- vac.per
  rc@sia.vac <- sia.vac

  return(rc)
}
