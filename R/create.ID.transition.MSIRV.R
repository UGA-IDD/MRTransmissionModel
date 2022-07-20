#' Method to create an ID.transition.MSIRV object
#'
#' @param n.age.class xxx
#' @param aging.rate xxx
#' @param survival.rate xxx
#' @param waifw xxx
#' @param routine.vac xxx
#' @param routine.vac.age.index xxx
#' @param maternal.obj xxx
#' @param time.step xxx
#' @param age.class xxx
#' @param birth.rate xxx
#' @param introduction.rate xxx
#' @param exponent xxx
#' @param frequency.dep xxx
#' @param is.stochastic xxx
#' @param get.births xxx
#'
#' @importFrom methods new
#'
#' @return xxx
#' @export
#'

create.ID.transition.MSIRV <- function(n.age.class,
                                       aging.rate,
                                       survival.rate,
                                       waifw,
                                       routine.vac = 0, #xxamy
                                       routine.vac.age.index = 12, #xxamy age index based on age classes
                                       maternal.obj, #maternal antibody object
                                       time.step = 1/2, #time step in months == generation time
                                       age.class = numeric(),
                                       birth.rate = 0,
                                       introduction.rate = 0,
                                       exponent = exponent,
                                       frequency.dep = F,
                                       is.stochastic = F,
                                       get.births = function(pop,tran){return(c(tran@birth.rate,rep(0,length(pop)-1)))}) {


  #First do all of the standard logic
  rc <- new("ID.transition.MSIRV",
            n.age.class = n.age.class,
            n.epi.class = 5,
            epi.class.label = c("M","S","I","R","V"),
            epi.class = rep(1:5, n.age.class),
            age.class = age.class,
            aging.rate = aging.rate,
            survival.rate = survival.rate,
            birth.rate = birth.rate,
            waifw = waifw,
            introduction.rate = introduction.rate,
            exponent = exponent,
            frequency.dep=frequency.dep,
            is.stochastic=is.stochastic,
            get.births = get.births)

  rc@m.inds = which(rc@epi.class==1)
  rc@s.inds = which(rc@epi.class==2)
  rc@i.inds = which(rc@epi.class==3)
  rc@r.inds = which(rc@epi.class==4)
  rc@v.inds = which(rc@epi.class==5)

  rc@age.surv.matrix <- matrix(0,nrow = n.age.class*rc@n.epi.class,
                               ncol = n.age.class*rc@n.epi.class)


  #for zeroing out all the appropriate transitions
  template.mtrx <- matrix(c(1,1,0,0,1,0,1,1,0,1,0,0,1,1,0,0,0,0,1,0,0,0,0,0,1), nrow=5, ncol=5)

  #useful low age
  low.age <-c(0, age.class[2:length(age.class)-1])

  #for each age class
  for (i in 1:n.age.class) {
    #first fill in the diagonal matrix
    tmp.template <- template.mtrx
    tmp.template[2,1] <- 0
    inds <- (i-1)*rc@n.epi.class+(1:rc@n.epi.class)
    rc@age.surv.matrix[inds,inds] <- (1-aging.rate[i])*survival.rate[i] *
      tmp.template

    #now fill in the off diagnal matrix
    if (i!=n.age.class) {
      #create a modified template matrix to move
      #people from the M class to the S class
      tmp.template <- template.mtrx

      #calculate the probability of losing maternal protection during
      #an age class
      p.mat.loss <- (pmaternal(low.age[i], maternal.obj)-
                       pmaternal(age.class[i], maternal.obj))/pmaternal(low.age[i],maternal.obj)

      #debug - Justin?
      p.mat.loss[pmaternal(low.age[i],maternal.obj)==0] <- 1
      #print(c(age.class[i],p.mat.loss))

      tmp.template[1,1] <- 1-p.mat.loss
      tmp.template[2,1] <- 1-tmp.template[1,1]

      #cat("age.class.i",age.class[i],"\n")
      #print(tmp.template)


      inds2 <- (i)*rc@n.epi.class+(1:rc@n.epi.class)
      rc@age.surv.matrix[inds2,inds] <- aging.rate[i]*survival.rate[i] *
        tmp.template
    }
  }

  #Next do all of the vaccination logic

  #get the %vaccinated in each age group per time step time step must be in years!
  #vac.per <- create.vacc.per.time.step.by.age(vpdf, age.class, time.step) #xxamy
  vac.per <- new("vacc.per.time.step.by.age", #xxamy
                 pvacc.in.age.class = numeric(length=n.age.class))
  vac.per@pvacc.in.age.class[routine.vac.age.index] <- routine.vac #xxamy
  rc@vac.per <- vac.per

  #print(range(vac.per@pvacc.in.age.class))

  return(rc)
}
