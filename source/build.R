## Functions for building out inputs, transitions, and outputs
## Authors: Amy Winter & Justin Lessler & Jess Metcalf


#### Maternal Immunity ####

#assumes that protection is the same as exponetial decay
setClass("maternal.exp.decay",
         representation(decay.rt = "numeric")#the exponential decay rate
)


#assume that there is some threshold age at which
#you are 100% protected if this age or younger
setClass("maternal.thresh",
         representation(thresh="numeric"))


#Method to get the percent immune by maternal antibodies.
#
#Parameters -
#    age - the age we want the decay rate for
#    mat.obj - the object describing maternally acquired immunity
setGeneric("pmaternal",
           function(age, mat.obj) standardGeneric("pmaternal"))

setMethod("pmaternal",
          c("numeric", "maternal.exp.decay"),
          function (age, mat.obj) {
            return(exp(-age*mat.obj@decay.rt))
          })


setMethod("pmaternal",
          c("numeric", "maternal.thresh"),
          function (age, mat.obj) {
            rc <- sapply(age, function(age) {
              if(age<=mat.obj@thresh) {return(1)}
              return(0)
            })
            return(rc)
          })


#### WAIFW ####

#Make a WAIFW matrix to mimic Emilia's - currently data-free
#
#Parameters -
#   age class boundries - the upper age limit for each age class in years
#   beta_young - beta on young individuls (<=13 years old)
#   beta_old - beta on old individuals (>13 years old)
#Returns -
#   a WAIFW matrix that mimics Emilia's 
get.vynnycky.WAIFW <- function (age.class.boundries = (1:120/12), beta_young, beta_old) {
  #get the force of infection for each age class.
  #make sure we use the mid point of each category
  n.age.cats <- length(age.class.boundries)
  ages.to.use <- (age.class.boundries +
                    c(0,age.class.boundries[2:n.age.cats-1]))/2
  
  #make the unscaled matrix version of the transmission over age
  foi.matrix <- matrix(1,nrow= n.age.cats,
                       ncol=n.age.cats)
  
  #Emilia's WAIFW has old and young betas, and contact between then is beta_old*0.7
  foi.matrix[ages.to.use<=13,ages.to.use<=13] <- beta_young
  foi.matrix[ages.to.use>13,ages.to.use>13] <- beta_old
  foi.matrix[ages.to.use<=13,ages.to.use>13] <- beta_old*0.7
  foi.matrix[ages.to.use>13,ages.to.use<=13] <- beta_old*0.7
  
  colnames(foi.matrix) <- age.class.boundries
  rownames(foi.matrix) <- age.class.boundries
  
  return(foi.matrix)
}

#Make a WAIFW matrix based on Polymod (ignoring rescaling by population
#size that might be sensible - just use raw contacts)
#
#Parameters -
#   age class boundries - the upper age limit for each age class in YEARS
#   the desired country - levels possible "all", "IT" "DE" "LU" "NL" "PL" "GB" "FI" "BE",
#         default is "GB"
#   bandwidth - desired smooth bandwidth - default=c(3,3)
#   do.touch - do you want to only include contacts involving touching? default is FALSE
#Returns -
#   a WAIFW matrix based on the Polymod results from chosen location with row and col
#   names indicating age classes
get.polymod.WAIFW <- function (age.class.boundries = (1:90),
                               country="GB",
                               bandwidth=c(3,3),
                               do.touch=FALSE) {
  require(KernSmooth)
  
  #bring in the polymod data on contacts for UK
  polymod <- read.csv("data/polymodRaw.csv")
  if (country!="all") polymod <- polymod[polymod$country==country,]
  
  #do touch only
  if (do.touch) polymod <- polymod[polymod[,"cnt_touch"]==1,]
  
  #obtain contacts, remove NAs, and make symmetrical by doubling up
  x <- cbind(polymod$participant_age, polymod$cnt_age_mean)
  x <- x[!is.na(x[,1]),]
  x <- x[!is.na(x[,2]),]
  x <- rbind(x,cbind(x[,2],x[,1]))
  
  #smooth monthly from 0 to 91 yrs of age and get the corresponding ages
  est<- bkde2D(x, bandwidth=bandwidth,gridsize=c(12*101, 12*101),range.x=list(c((1/24),100),c((1/24),100)))
  ages.polymod.smooth <- est$x1
  
  #image(ages.polymod.smooth,ages.polymod.smooth,est$fhat,xlim=c(0,10),ylim=c(0,10))
  
  #find which fitted ages (ages.polymod.smooth) each of the desired lower boundaries is between
  n.age.cats <- length(age.class.boundries)
  lowpoints <-  c(0,age.class.boundries[2:n.age.cats-1])
  index <- (findInterval(lowpoints,ages.polymod.smooth));
  #set very large ages to all be the same as the largest age
  index[index>=length(ages.polymod.smooth)]=length(ages.polymod.smooth)-1
  
  #extract appropriate matrix, taking smoothed estimates for the ranges
  foi.matrix <- est$fhat[index+1,index+1]
  
  #adjust for fact that the width of your age class should not affect the number of contacts you make
  #foi.matrix <- foi.matrix/diff(c(0,age.class.boundries))
  
  colnames(foi.matrix) <- age.class.boundries
  rownames(foi.matrix) <- age.class.boundries
  return(foi.matrix)
}

#Make a flat WAIFW matrix
#
#Parameters -
#   age class boundries - the upper age limit for each age class in years,
#Returns -
#   a flat WAIFW matrix
get.flat.WAIFW <- function (age.class.boundries = (1:120/12)) {
  n.age.cats <- length(age.class.boundries)
  ages.to.use <- (age.class.boundries +
                    c(0,age.class.boundries[2:n.age.cats-1]))/2
  
  foi.matrix <- matrix(1,nrow= n.age.cats,
                       ncol=n.age.cats)
  
  colnames(foi.matrix) <- age.class.boundries
  rownames(foi.matrix) <- age.class.boundries
  
  return(foi.matrix)
}


#Make a diagonal WAIFW matrix
#
#Parameters -
#   age class boundries - the upper age limit for each age class in years,
#Returns -
#   a diagonal WAIFW matrix
get.diagonal.WAIFW <- function (age.class.boundries = (1:120/12)) {
  n.age.cats <- length(age.class.boundries)
  ages.to.use <- (age.class.boundries +
                    c(0,age.class.boundries[2:n.age.cats-1]))/2
  
  foi.matrix <- matrix(0,nrow= n.age.cats,
                       ncol=n.age.cats)
  diag(foi.matrix) <- 1
  colnames(foi.matrix) <- age.class.boundries
  rownames(foi.matrix) <- age.class.boundries
  
  # fully diagonal perhaps a bit extreme?
  #  here can infect 6months above and below (at least with current age config)
  #foi.matrix[cbind(2:n.age.cats,1:(n.age.cats-1))] <- 1
  #foi.matrix[cbind(1:(n.age.cats-1),2:n.age.cats)] <- 1
  
  ## Updated:
  # age structure for < 15y is monthly, so let them infect if they are the same year
  # ie: 1-12m with 1-12m, 13-24m with 13-24m, etc.
  
  for(i in 1:15){
    for(j in 1:12){
      
      foi.matrix[cbind(rep((12*(i-1)+1):(12*(i-1)+12),12),c(rep((12*(i-1)+j),12)))] <- 1
      
    }}
  
  
  
  return(foi.matrix)
}


# Make a parametric smooth WAIFW
# from Farrington, J. Amer. Stat. Assoc.; 2005, 100 p370;
#
#  parameters -
#     age class boundries - the upper age limit for each age class in years,
#     parameters defining the shape:
#        mu: age of highest contact increases with mu
#        gam: width around equal age diagonal increases with gam
#        sig: decreases strength in other diagonal (shrinks high trans int)
#        delta:  background homogeneous contact rate
#
#Returns -
#   a smooth WAIFW matrix
get.smooth.WAIFW<-function(age.class.boundries = (1:120/12),
                           mu=12.71,sig=0.69, gam=0.17, delta=0){
  
  n.age.cats <- length(age.class.boundries)
  ages.to.use <- (age.class.boundries +
                    c(0,age.class.boundries[2:n.age.cats-1]))/2
  
  #gamma (p371)
  gam.func <- function(x,y,mu,sig){
    u <- (x+y)/(sqrt(2))
    vee <- 1/(sig^2)
    cval <- (sqrt(2)*mu*(1-(1/vee)))^(vee-1)
    cval <- cval*exp(1-vee)
    if (vee<1) cval <- 1
    gamma <- (1/cval)*(u^(vee-1))*exp((-vee*u)/(sqrt(2)*mu))
    return(gamma)
  }
  
  #b (p371); set alpha=beta here since always want symmetrical matrix
  b.func <- function(x,y,gam){
    u <- (x+y)/(sqrt(2))
    v <- (x-y)/(sqrt(2))
    alpha<-beta<-(1-gam)/(2*gam)
    b <- (((u+v)^(alpha-1))*((u-v)^(beta-1)))/(u^(alpha+beta-2))
    return(b)
  }
  
  #smooth value
  beta.func <- function(x,y,mu,sig, gam, delta){
    gam.func(x,y,mu,sig)*b.func(x,y,gam)+delta
  }
  
  betas <- outer(X=ages.to.use,Y=ages.to.use,
                 FUN=beta.func,mu=mu,sig=sig,gam=gam,delta=delta)
  
  return(betas)
}

#Make a WAIFW matrix based on Prem et al. 2021 
#using pakistan for afghanistan
#
#Parameters -
#   age class boundries - the upper age limit for each age class in YEARS
#   uncode - country code
#   bandwidth - desired smooth bandwidth - default=c(3,3)
#   adjustment_1980 - for 202110gavi_v3 given that Shaun's R0 estimates are adjusted to 1980 pop structure
#Returns -
#   a WAIFW matrix based on the Polymod results from chosen location with row and col
#   names indicating age classes
get.prem.WAIFW <- function (age.class.boundries = (1:90),
                            uncode, other.contact.matrix=F,
                            bandwidth=c(3,3), adjustment_start_year=F, year=1980) {
  
  if (uncode==583) iso3code <- "ETH" #FSM - Somalia becomes Ethiopia
  if (uncode==332) iso3code <- "DOM" #HTI - Haiti becomes Dominican Republic
  if (uncode==296) iso3code <- "FJI" #KIR - Kiribati becomes Fuji
  if (uncode==584) iso3code <- "SLB" #MHL - Marshall Islands becomes Solomon Islands
  if (uncode==706) iso3code <- "ETH" #SOM - Somalia becomes Ethiopia
  if (uncode==798) iso3code <- "TON" #TUV - Tuvalu becomes Tonga
  if (uncode==999) iso3code <- "ALB" #XK - Kosovo becomes Albania
  if (!uncode %in% c(583, 332, 296, 584, 706, 798, 999)) {
    iso3code <- countrycode::countrycode(uncode, origin="un", destination="iso3c")
    #name <- countrycode::countrycode(iso3code, origin="iso3c", destination="country.name")
  }
  
  load("./data/prem_contact_matrices_2021/contact_all.rdata")
  index <- which(names(contact_all)==iso3code)
  if (length(index)==0) {
    stop(paste("missing contact matrix for iso3code:", iso3code, "see get.prem.WAIFW()"))
  } else {
    df.prem <- contact_all[[index]]
  }
  
  #turn data into long form (2 columns of ages)
  prem <- round(df.prem*1000)
  prem.v <- as.numeric(unlist(c(prem)))
  if (prem.v[length(prem.v)]==0) prem.v[length(prem.v)] <- 1 #hack, code breaks if end on a zero.
  prem.sum <- as.numeric(cumsum(prem.v))
  prem.ages <- seq(3, 78, 5)
  x <- matrix(NA, sum(prem.v), 2)
  for (c in 1:length(prem.v)){
    index <- min(which(is.na(x[,1])))
    index2 <- prem.sum[c]
    x[(index:index2), 1] <- prem.ages[ceiling(c/16)]
    if (!any(c==seq(16,256,16))) x[(index:index2), 2] <- prem.ages[c-(floor(c/16)*16)] #prem.ages[c-(floor(c/16)*16)+1]
    if (c%in%seq(16,256,16)) x[(index:index2), 2] <- prem.ages[(floor(c/16)*16)-(16*(c/16-1))]
  }
  
  if (other.contact.matrix){
    stop("other.contact.matrix=T disabled beginning 202110gavi_v3 - need to fix")
    #if (name=="China") name <- "China_Read2014"
    #if(substr(name, 1,3)<"Mos") df.prem <- read_excel("./data/prem_contact_matrices/MUestimates_all_locations_1.xls", sheet=name)
    #if(substr(name, 1,3)>"Mos") df.prem <- read_excel("./data/prem_contact_matrices/MUestimates_all_locations_2.xls", sheet=name)
    #prem <- round(df.prem*100)
    #prem.v <- as.numeric(unlist(c(prem)))
    #if (prem.v[length(prem.v)]==0) prem.v[length(prem.v)] <- 1 #hack, code breaks if end on a zero.
    #prem.sum <- as.numeric(cumsum(prem.v))
    #ages.min <- c(1,6,20,65)
    #ages.max <- c(5,19,64,80)
    #x <- matrix(NA, sum(prem.v), 2)
    #for (c in 1:length(prem.v)){
    #    index <- min(which(is.na(x[,1])))
    #    index2 <- prem.sum[c]
    #    x[(index:index2), 1] <- ages.min[c]
    #    x[(index:index2), 2] <- ages.max[c]
    #}
  }
  
  #smooth monthly from 0 to 91 yrs of age and get the corresponding ages
  est <- bkde2D(x, bandwidth=bandwidth, gridsize=c(12*101, 12*101), range.x=list(c((1/24),100),c((1/24),100)))
  ages.polymod.smooth <- est$x1
  
  #image(ages.polymod.smooth,ages.polymod.smooth,est$fhat,xlim=c(0,100),ylim=c(0,100))
  
  #find which fitted ages (ages.polymod.smooth) each of the desired lower boundaries is between
  n.age.cats <- length(age.class.boundries)
  lowpoints <-  c(0,age.class.boundries[2:n.age.cats-1])
  index <- (findInterval(lowpoints,ages.polymod.smooth));
  #set very large ages to all be the same as the largest age
  index[index>=length(ages.polymod.smooth)]=length(ages.polymod.smooth)-1
  
  #extract appropriate matrix, taking smoothed estimates for the ranges
  foi.matrix <- est$fhat[index+1,index+1]
  #lattice::levelplot(foi.matrix*100)
  
  #adjust for fact that the width of your age class should not affect the number of contacts you make
  #foi.matrix <- foi.matrix/diff(c(0,age.class.boundries))
  
  colnames(foi.matrix) <- age.class.boundries
  rownames(foi.matrix) <- age.class.boundries
  
  #Prem_2020_Age_Contacts* startyear_Age_Structure / Prem_2020_Age_Structure = startyear_Age_Contacts
  if (adjustment_start_year) {
    demog.tmp <- getDemography.wpp2017(uncode, age.classes=age.class.boundries, if.missing.use.region=T)
    pop.start.year <- demog.tmp$pop.age.byageclasses.1950.2100[,which(colnames(demog.tmp$pop.age.byageclasses.1950.2100)==as.character(year))]
    pop.2020 <- demog.tmp$pop.age.byageclasses.1950.2100[,which(colnames(demog.tmp$pop.age.byageclasses.1950.2100)=="2020")]
    foi.matrix <- foi.matrix*(pop.start.year/pop.2020)
    #lattice::levelplot(foi.matrix*100)
  }
  
  return(foi.matrix)
}

#Function to scale the WAIFW matrix to a particular
#R0 given a particular population
#
#Parameters -
#    R0 - the reproductive rate to scale to
#    state - a population at disease free equilibrium
#    waifw - the matrix to scale
#
#Returns -
#    a scaled vertions of waifw
scale.WAIFW <- function(R0, state, waifw, frequency.dep=F,suscept.state = 1) {
  
  if (frequency.dep) denom <- sum(state[,1]) else denom <- 1
  
  #print(length(state[,1]))
  
  #move everyone into susceptible category for DFE
  if (length(state@epi.class)==length(state[,1])) {
    DFE.state <- matrix(0,length(state[,1]),1)
    for (k in  1:state@n.epi.class)
      DFE.state[state@epi.class==suscept.state,1] <-
        DFE.state[state@epi.class==suscept.state,1]+state[state@epi.class==k,1]
  } else {
    epi.class.here <- rep(state@epi.class,2)
    DFE.state <- matrix(0,length(epi.class.here),1)
    for (k in  1:state@n.epi.class) {
      DFE.state[epi.class.here == suscept.state,] <- DFE.state[epi.class.here == suscept.state,]+
        state[epi.class.here == k,]
      # print(k)
    }
    #add the genders together
    DFE.state <- DFE.state[1:length(state@epi.class),]+DFE.state[(length(state@epi.class)+1):(2*length(state@epi.class)),]
    DFE.state <- matrix(DFE.state,length(state@epi.class),1)
  }
  
  
  #worked kinda
  #next.gen <- (state[state@epi.class==suscept.state,1]*waifw)/denom
  
  #plot(DFE.state[state@epi.class==suscept.state,1],type="l")
  
  
  #more correct
  next.gen <- DFE.state[state@epi.class==suscept.state,1]*(1-exp(-waifw/denom))
  
  #get the first eigen value
  cur.R0 <- Re(eigen(next.gen)$value[1])
  
  #More correct transform
  R.ratio <- R0/cur.R0; #print(R0); #print(cur.R0); #print(R.ratio)
  waifw <- -log(1-R.ratio*(1-exp(-waifw/denom)))*denom
  
  return(waifw)
}


#### Vaccination ####

#Class that holds the information on how many people should be
#vaccinated during each time step in each age class
#in the TSIR framework.
setClass("vacc.per.time.step.by.age",
         representation(pvacc.in.age.class = "numeric") #vector of how many are
         #vaccinated in the age
         #class per time step
)


## object for the CDF of vaccination by age (non scaled)
setClass("vaccine.cdf.byage",
         representation(cdf = "numeric", #vector of length ages
                        ages = "numeric" #vector of length cdf
         )
)

#Make the generic vaccine success object
setClass("vsucc")

#make a constants vaccine success object
setClass("vsucc.constant",
         representation(success.rate="numeric"),
         contains="vsucc")

#Make an object to hold a logistic function
#for vaccine success
setClass("vsucc.Logistic",
         representation(intercept = "numeric",
                        mo.eff = "numeric",
                        full.efficacy = "numeric" #the max reachable efficacy
         ),
         contains = "vsucc"
)

#register the pvaccsuccess generic
#
#Parameters -
#   age - the age in months we are checking
#   vs.obj - the vaccince success object
#
#Returns -
#   the probability of successful vaccination at that
#   age
setGeneric("pvacsuccess",
           function(age, vsucc.obj) standardGeneric("pvacsuccess"))

#default method for vaccine success, just returns 1 (vaccine is always successful)
setMethod("pvacsuccess",
          c("numeric","vsucc.constant"),
          function (age, vsucc.obj) {
            
            prob.vsucc <- new("prob.vsucc.byage",
                              prob.vsucc=rep(vsucc.obj@success.rate, length(age)),
                              ages=age)
            
            return(prob.vsucc)
          })

#set the method for the vsucc logistic
setMethod("pvacsuccess",
          c("numeric","vsucc.Logistic"),
          function (age, vsucc.obj) {
            rc <- plogis(vsucc.obj@mo.eff*age,
                         -vsucc.obj@intercept)
            rc[rc>vsucc.obj@full.efficacy] <- vsucc.obj@full.efficacy
            
            prob.vsucc <- new("prob.vsucc.byage",
                              prob.vsucc=rc,
                              ages=age)
            return(prob.vsucc)
          }
)


## object for the prob of vaccine success by age
setClass("prob.vsucc.byage",
         representation(prob.vsucc = "numeric", #vector of length ages
                        ages = "numeric" #vector of length prob
         )
)



#Function to find the non-scaled vaccination CDF (assuming uniform pdf)
#
#Parameters -
#     min.age - numeric - age in months  
#     max.age - numeric - age in months  
#
#Outputs - CDF vaccinated across age in months 
get.vcdf.uniform <- function(min.age=25, max.age=36){
  
  ages <- min.age:max.age
  cdf <- cumsum(rep(1/length(ages), length(ages)))
  
  vcdf <- new("vaccine.cdf.byage",
              cdf=cdf, ages=ages)
  
  return(vcdf)
  
}


#Function to find the non-scaled vaccination CDF (assuming normal pdf)
#
#Parameters -
#     min.age - numeric - age in months  
#     max.age - numeric - age in months  
#
#Outputs - CDF vaccinated across age in months 
get.vcdf.normal <- function(min.age=25, max.age=36){
  
  ages <- min.age:max.age
  age.range <- length(ages)
  cdf <- c(pnorm(seq(-2,2,4/age.range),0,1)[-1])
  
  vcdf <- new("vaccine.cdf.byage",
              cdf=cdf, ages=ages)
  
  return(vcdf)
  
}

#Data sets on vaccine efficacy
vdata.gans1998 <- data.frame(age=c(6,9,12), successes=c(10,17,21),
                             N=c(23,20,22))


#rubella with MMR
vdata.MMR.boulianne1995 <- data.frame(age=c(12,13,14,15), successes=c(62,74,50,47),
                                      N=c(64,79,50,48))


#get the Boulianne vauccine success object (rubella)
get.boulianne.vsucc <- function() {
  return(fit.vsucc.Logistic(vdata.MMR.boulianne1995,.97))
}

#get the gans vauccine success object
get.gans.vsucc <- function() {
  return(fit.vsucc.Logistic(vdata.gans1998,.97))
}

#fits a vsucc.Logistic object to
#some data using logistic regression
#
#Parameters -
#   data - a data frame of age (in months), number of
#          success and total tests
#   full.efficacy - the maximum efficacy we wish to allow
#Returns -
#   a vsucc.Logistic object
fit.vsucc.Logistic <- function(data, full.efficacy) {
  succ.fail <- c()
  age <- c()
  for (i in 1:nrow(data)) {
    succ.fail <- c(succ.fail, rep (1,data$successes[i]))
    succ.fail <- c(succ.fail, rep (0,data$N[i]-data$successes[i]))
    age <- c(age, rep(data$age[i], data$N[i]))
  }
  
  mdl<-glm(succ.fail~age, family=binomial(link="logit"))
  
  res <- new("vsucc.Logistic",
             intercept=as.numeric(mdl$coef[1]),
             mo.eff=as.numeric(mdl$coef[2]),
             full.efficacy=as.numeric(full.efficacy))
  
  return(res)
}

setClass("cov.estimates",
         representation(sia.cov = "numeric", 
                        sia.min.age = "numeric", 
                        sia.max.age = "numeric", 
                        MR1.cov = "numeric", 
                        MR2.cov = "numeric",
                        inaccessible.prop = "numeric"
         ))


get.routine.time.age.specific <- function(time.step=0.5, age.classes = c(1:240, seq(252,1212,12)),
                                          time.specific.MR1cov = rep(0.6, 50), age.min.MR1=rep(12, 50), age.max.MR1=rep(23,50), 
                                          time.specific.MR2cov=rep(0.5, 50), age.min.MR2=rep(24, 50), age.max.MR2=rep(35, 50), 
                                          obj.vcdf.MR1=get.vcdf.uniform(12, 23), obj.vcdf.MR2=get.vcdf.uniform(24, 35), 
                                          obj.prob.vsucc = pvacsuccess(1:(14*12), new("vsucc.constant", success.rate=0.8)), 
                                          MR1MR2correlation=F){
  
  time.specific.routine <- prop.fail.MR1.byage <- prop.fail.MR2.byage <- matrix(0, length(time.specific.MR1cov), length(age.classes))
  prop.fail.MR1 <- one.minus.ve1 <-  prop.fail.MR2 <- one.minus.ve2 <- rep(0, length(time.specific.MR1cov))
  
  #get probability of successful MR1
  for (j in 1:length(time.specific.MR1cov)) {
    if (time.specific.MR1cov[j]!=0){
      
      cdf.vaccination.MR1 <- obj.vcdf.MR1@cdf[obj.vcdf.MR1@ages %in% c(age.min.MR1[j]:age.max.MR1[j])]
      prob.success.MR1 <- obj.prob.vsucc@prob.vsucc[obj.prob.vsucc@ages %in% c(age.min.MR1[j]:age.max.MR1[j])]
      if (length(cdf.vaccination.MR1)!=length(prob.success.MR1)) stop("problem matching age in vcdf and vsucc, MR1")
      cdf.scaled.vaccination <- time.specific.MR1cov[j]*(cdf.vaccination.MR1/cdf.vaccination.MR1[length(cdf.vaccination.MR1)])
      cdf.scaled.vaccination <- c(0,cdf.scaled.vaccination)
      pdf <- diff(cdf.scaled.vaccination)
      h <- pdf/(1-cdf.scaled.vaccination[2:length(cdf.scaled.vaccination)-1]) #hazard
      low.age <- (age.min.MR1[j]-1):(age.max.MR1[j]-1)
      high.age <- age.min.MR1[j]:age.max.MR1[j]
      age.sz <- high.age-low.age
      ts.per.class <- age.sz/time.step
      final.pdf <- 1-(1-h)^(1/ts.per.class)
      pdf.scaled.vaccination <- pmax(final.pdf,0)
      time.specific.routine[j,age.classes %in% c(age.min.MR1[j]:age.max.MR1[j])] <-  pdf.scaled.vaccination*prob.success.MR1 #pdf*prob.success.MR1 = rep(0.04, 12); sum(rep(0.04, 12)) = 0.48 = ve1*mcv1, check
      #time.specific.routine[j,age.classes %in% c(age.min.MR1[j]:age.max.MR1[j])] <-  pdf*prob.success.MR1
      one.minus.ve1[j] <- sum(pdf/sum(pdf)*(1-prob.success.MR1)) #using pdf/sum(pdf) to get weighted average 1-VE1 = 0.2, check 
      prop.fail.MR1[j] <- sum(pdf*(1-prob.success.MR1)) # = 0.12 = (1-ve1)*mcv1, proportion of total population with primary failure from MR1, check
      prop.fail.MR1.byage[j,age.classes %in% c(age.min.MR1[j]:age.max.MR1[j])] <- pdf*(1-prob.success.MR1)
      
    }
  }
  
  one.minus.ve1.lagged <- c(one.minus.ve1[1], one.minus.ve1[-length(one.minus.ve1)]) #lag one year so that it applies to the correct birth cohort
  prop.MR1failANDnoMR1.laggged <- (1-c(time.specific.MR1cov[1], time.specific.MR1cov[-length(time.specific.MR1cov)])) + c(prop.fail.MR1[1], prop.fail.MR1[-length(prop.fail.MR1)]) #lag one year so that it applies to the correct birth cohort = 0.52, check
  
  
  #get probability of successful MR2
  for (j in 1:length(time.specific.MR2cov)) {
    if (time.specific.MR2cov[j]!=0){
      cdf.vaccination.MR2 <- obj.vcdf.MR2@cdf[obj.vcdf.MR2@ages %in% c(age.min.MR2[j]:age.max.MR2[j])]
      prob.success.MR2 <- obj.prob.vsucc@prob.vsucc[obj.prob.vsucc@ages %in% c(age.min.MR2[j]:age.max.MR2[j])]
      if (length(cdf.vaccination.MR2)!=length(prob.success.MR2)) stop("prablem matching age in vcdf and vsucc, MR2")
      age.range <- length(age.min.MR2[j]:age.max.MR2[j])
      #if (MR1MR2correlation) cdf.scaled.vaccination <- one.minus.ve1.lagged[j]*time.specific.MR2cov[j]*(cdf.vaccination.MR2/cdf.vaccination.MR2[length(cdf.vaccination.MR2)])/prop.MR1failANDnoMR1.laggged[j] #last index is 0.19230769 = 0.1/0.52 = (mcv2*(1-ve1)) / (1-ve1*mcv1), proportion of susceptibles vaccinated with MCV2, assuming correlation, gtg
      if (MR1MR2correlation) cdf.scaled.vaccination <- one.minus.ve1.lagged[j]*time.specific.MR2cov[j]*(cdf.vaccination.MR2/cdf.vaccination.MR2[length(cdf.vaccination.MR2)]) #last index is 0.1 = (mcv2*(1-ve1)), proportion of total population vaccinated with MCV2, assuming correlation, gtg
      if (!MR1MR2correlation) cdf.scaled.vaccination <- time.specific.MR2cov[j]*(cdf.vaccination.MR2/cdf.vaccination.MR2[length(cdf.vaccination.MR2)]) #last index is 0.5 = mcv2, proportion of total population (or susceptibles) vaccinated with MCV2, no correlation, gtg
      cdf.scaled.vaccination <- c(0,cdf.scaled.vaccination) 
      pdf <- diff(cdf.scaled.vaccination)
      h <- pdf/(1-cdf.scaled.vaccination[2:length(cdf.scaled.vaccination)-1]) #hazard
      low.age <- (age.min.MR2[j]-1):(age.max.MR2[j]-1)
      high.age <- age.min.MR2[j]:age.max.MR2[j]
      age.sz <- high.age-low.age
      ts.per.class <- age.sz/time.step
      final.pdf <- 1-(1-h)^(1/ts.per.class)
      pdf.scaled.vaccination <- pmax(final.pdf,0)
      time.specific.routine[j,age.classes %in% c(age.min.MR2[j]:age.max.MR2[j])] <- pdf.scaled.vaccination*prob.success.MR2 
      #time.specific.routine[j,age.classes %in% c(age.min.MR2[j]:age.max.MR2[j])] <- pdf*prob.success.MR2 
      #IF assuming correlation, pdf*prob.success.MR2 = rep(0.1282051, 12); sum(rep(0.1282051, 12)) = 0.1538462 = 0.08/0.52 = (mcv2*ve2*(1-ve1)) / (1-ve1*mcv1), proportion of susceptibles successfully vaccinated with MCV2, gtg
      #IF assuming no correlation pdf*prob.success.MR2 = rep(0.333, 12); sum(rep(0.333, 12)) = 0.4 = 0.5*0.8 = mcv2*ve2, check
      one.minus.ve2[j] <- sum(pdf/sum(pdf)*(1-prob.success.MR2)) # = 0.2 = 1-ve2, check
      #prop.fail.MR2[j] <- sum(pdf*(1-prob.success.MR2))*prop.MR1failANDnoMR1.laggged[j] # 0.02, proportion of total population with primary failure from MR2, gtg
      prop.fail.MR2[j] <- sum(pdf*(1-prob.success.MR2)) # 0.02, proportion of total population with primary failure from MR2, gtg
      prop.fail.MR2.byage[j,age.classes %in% c(age.min.MR2[j]:age.max.MR2[j])] <- pdf*(1-prob.success.MR2)
    }
  }
  
  
  return(list(
    age.time.specific.routine=time.specific.routine, 
    one.minus.ve1=one.minus.ve1,
    one.minus.ve2=one.minus.ve2,
    prop.fail.MR1=prop.fail.MR1, 
    prop.fail.MR1.byage = prop.fail.MR1.byage,
    prop.fail.MR2=prop.fail.MR2,
    prop.fail.MR2.byage = prop.fail.MR2.byage))
  
}

get.sia.time.age.specific <- function(age.classes=c(1:240, seq(252,1212,12)), 
                                      time.specific.SIAcov=c(rep(0,29),0.5,rep(0,19),0.7),
                                      age.min.sia=rep(12,50), 
                                      age.max.sia=rep(60,50), 
                                      obj.prob.vsucc=pvacsuccess(1:240, new("vsucc.constant", success.rate=0.85))){
  
  ## Get SIA age specific estimates for each year
  age.time.specific.SIA <-  matrix(0, length(time.specific.SIAcov), length(age.classes))
  prop.fail.SIA <- rep(0, length(time.specific.SIAcov))
  
  for (j in 1:length(time.specific.SIAcov)){
    if (time.specific.SIAcov[j]!=0){
      prob.success.SIA <- obj.prob.vsucc@prob.vsucc[obj.prob.vsucc@ages %in% c(age.min.sia[j]:age.max.sia[j])]
      vacc.SIA <- rep(time.specific.SIAcov[j], length(age.min.sia[j]:age.max.sia[j]))
      age.time.specific.SIA[j,(age.min.sia[j]:age.max.sia[j])] <- vacc.SIA*prob.success.SIA
      prop.fail.SIA[j] <- (vacc.SIA*(1-prob.success.SIA))[1]
    }
  }
  
  return(list(age.time.specific.SIA=age.time.specific.SIA,
              prop.fail.SIA=prop.fail.SIA))
  
}


#Function to find the non-scaled MR1 CDF (uses most recent DHS data as of 2019)
#
#Parameters -
#     uncode- numeric - country UN code
#     min.age - numeric - age in months
#     max.age - numeric - age in months
#
#Outputs - CDF vaccinated across age in months 
get.MR1cdf.survival <- function(uncode, min.age=1, max.age=24){
  
  data <- read.csv("./data/DHS_ageMCV1_Data.csv")
  dat <- data[which(data$uncode==uncode),]
  ages <- min.age:max.age
  if (nrow(dat)>0){ #if DHS data for specific country - use country data
    obj <- survival::Surv(time=dat$age, time2=rep(max(dat$age), length(dat$age)), event=as.numeric(dat$censoring), type="interval")
    out <- summary(survfit(obj~1))
    sm <- smooth.spline(c(0,out$time), c(0,(1-out$surv)))
    cdf <- predict(sm, ages)$y
  } else { #otherwise use region data
    cc <- read.csv("./data/country_codes.csv")
    reg <- cc$Region_Code[which(cc$uncode==uncode)]
    if (reg=="AFRO"){
      dat <- read.csv("./data/DHS_ageMCV1_Data_AFRO.csv")
    } else if (reg=="EMRO"){
      dat <- read.csv("./data/DHS_ageMCV1_Data_EMRO.csv")
    } else if (reg=="SEARO"){
      dat <- read.csv("./data/DHS_ageMCV1_Data_SEARO.csv")
    } else{
      dat <- data[data$region==reg,]
    }
    obj <- survival::Surv(time=dat$age, time2=rep(max(dat$age), length(dat$age)), event=as.numeric(dat$censoring), type="interval")
    out <- summary(survfit(obj~1))
    sm <- smooth.spline(c(0,out$time), c(0,(1-out$surv)))
    cdf <- predict(sm, ages)$y
  }
  
  vcdf <- new("vaccine.cdf.byage",
              cdf=cdf, ages=ages)
  
  return(vcdf)
}



#### Seasonality ####

#register a function for getting the seasonal adjustment for
#a time or series of times based on a seasonality object
#
#Parameter:
#     t.in.yr - time of the year in fraction of 1 year
#     season.obj - an object encoding the seasonality
#
#
#Returns:
#   the seasonal multiplier, or series thereof
setGeneric("get.seasonal.mult",
           function(t.in.yr, season.obj) standardGeneric("get.seasonal.mult"))


#generic virtual class for seasonal forcing
setClass("seasonal.force")

setClass("seasonal.cosine",
         representation(amplitude = "numeric", #amplitude of seasonal forcing
                        offset = "numeric", #offset from t=0 for nadir
                        period = "numeric"), #number of peaks per year
         prototype = list(amplitude=.2, offset = 0, period = 1),
         contains="seasonal.force"
)

setClass("seasonal.periodic.bspline",
         representation(parameters = "numeric",
                        nbasis = "numeric",
                        degree = "numeric",
                        period = "numeric"
         ),
         prototype = list(parameters= c(0.5,1,1), nbasis=3, degree=3, period=1),
         contains="seasonal.force",
)

setClass("seasonal.age.specific",
         representation(amplitude = "numeric", #amplitude of seasonal forcing
                        offset = "numeric", #offset from t=0 for nadir
                        period = "numeric", #number of peaks per year
                        age.class = "numeric",#age classes in the WAIFW
                        chosen.age = "numeric"),#vector of length 2, indicating start and end of seasonally forced age class
         prototype = list(amplitude=.2, offset = 0, period = 1, age.class=c(1:60,60*12),chosen.age=c(50,60)),
         contains="seasonal.force",
)

#designed to take output from TSIR-type analysis
setClass("seasonal.piecewise.constant",
         representation(parameters = "numeric",
                        time.in.year = "numeric"),
         prototype = list(parameters = rep(1,24),time.in.year = (0:24)/24),
         contains="seasonal.force",
)

#designed to take output from TSIR-type analysis
setClass("seasonal.age.specific.piecewise.constant",
         representation(parameters = "numeric",
                        time.in.year = "numeric",
                        age.class = "numeric",#age classes in the WAIFW
                        chosen.age = "numeric"),#vector of length 2, indicating start and end of seasonally forced age class
         prototype = list(parameters = rep(1,24),time.in.year = (0:24)/24),
         contains="seasonal.force",
)


setMethod("get.seasonal.mult",
          c("numeric","seasonal.cosine"),
          function(t.in.yr, season.obj) {
            mult <- 1+ season.obj@amplitude *
              cos(2*season.obj@period*pi*(t.in.yr+season.obj@offset))
            return(mult)
          })

setMethod("get.seasonal.mult",
          c("numeric","seasonal.periodic.bspline"),
          function(t.in.yr, season.obj) {
            require(pomp)
            y <- periodic.bspline.basis(t.in.yr,nbasis=season.obj@nbasis,
                                        degree=season.obj@degree,period=season.obj@period)
            mult <- y%*%season.obj@parameters
            return(mult)
          })


setMethod("get.seasonal.mult",
          c("numeric","seasonal.piecewise.constant"),
          function(t.in.yr, season.obj) {
            mult <- season.obj@parameters[findInterval(t.in.yr%%1,season.obj@time.in.year)]
            return(mult)
          })


setMethod("get.seasonal.mult",
          c("numeric","seasonal.age.specific"),
          function(t.in.yr, season.obj) {
            mult.vec <- 1+ season.obj@amplitude *
              cos(2*season.obj@period*pi*(t.in.yr+season.obj@offset))
            
            nage <- length(season.obj@age.class)
            change.age <- which(season.obj@age.class>=season.obj@chosen.age[1] &
                                  season.obj@age.class<=season.obj@chosen.age[2])
            
            mult <- array(dim=c(nage,nage,length(t.in.yr)))
            for (j in 1:length(t.in.yr)) {
              mult[,,j] <- matrix(1,nage,nage)
              mult[,,j][change.age,change.age] <- mult.vec[j]
            }
            return(mult)
          })


setMethod("get.seasonal.mult",
          c("numeric","seasonal.age.specific.piecewise.constant"),
          function(t.in.yr, season.obj) {
            mult.vec <- season.obj@parameters[findInterval(t.in.yr%%1,season.obj@time.in.year)]
            
            nage <- length(season.obj@age.class)
            change.age <- which(season.obj@age.class>=season.obj@chosen.age[1] &
                                  season.obj@age.class<=season.obj@chosen.age[2])
            
            mult <- array(dim=c(nage,nage,length(t.in.yr)))
            for (j in 1:length(t.in.yr)) {
              mult[,,j] <- matrix(1,nage,nage)
              mult[,,j][change.age,change.age] <- mult.vec[j]
            }
            
            return(mult)
          })

#### Demography ####

#Class definition for transition matrix. This matrix knows how to shift the
#ID.state.matrix from one time step to another, and contains all the fun stuff
#we need to do the transition (i.e., the matrix)
setClass("nMx",
         representation(rates="data.frame", #age specific death rates per 1 from 1950 to 2100 by 5 year increments
                        mid.age="numeric" #mid-age associated with rates,
         ))

#Function to pull in Age Specific Death Rates and return survivorship for running out transients
#
#Paramaters -
#   age.classes - a set of age classes for which tran is being built
#   generation time - generation time
#   nMx - object of age specific death rates over time and mid-age per rate
#   year - year to pull age specific death rates that will be used to run out the transients (should coincide with DFE year)
#
#Returns-
#   returns the age profile of survivorship in units of the generation time
create.surv.prob.over.age <- function(age.classes, generation.time, check=F, nMx=NULL, year=1990){
  
  rate.years <- seq(1950,2100,5)
  index <- min(which(findInterval(rate.years, year)==1))
  rates.overtime <- nMx@rates
  mid.age <- nMx@mid.age
  generation.time.year <- generation.time/12   # put generation.time in terms of years rather than months (0.5 in months) now (0.041667 in years)
  
  yearly.mortality <- rates.overtime[,index]
  
  # Fit a smooth spline to this
  fit <- smooth.spline(mid.age, log(yearly.mortality), df=length(mid.age))
  
  # Convert the fit to our ages
  survs <- 1-exp(predict(fit, age.classes/12)$y) 
  # gives me only the $y or estiamted m.x.n for age.classes/12 (puts age.classes in years rather than months)
  # for each age predict the survival rate based on the smooth.spline
  # exp to get back to m.x.n from it's previous log state and take 1-m.x.n to get survival rate
  #plot(survs)
  
  # Convert yearly rates to generation time rates
  survs <- 1+(log(survs)/(1/generation.time.year)) 
  
  # Need to adjust last age class to deplete
  survs[length(survs)] <-  0.5^(12*1/generation.time)
  
  if (check) {
    plot(mid.age, yearly.mortality, type="b", xlim=c(0,100) ,ylim=range(yearly.mortality[mid.age<100]),pch=19)
    points(age.classes/12,(1-survs^(1/generation.time.year)), col=2,type="b")
    points(fit$x,exp(fit$y),col=4, type="l")
  }
  
  return(survs)
}

#Function to get starting state and transition object
#
#Parameters -
#     uncode - UN country code
#     generation.time - the desired generation time in months
#     age.classes - vector - the upper limit of the age classes we are interested in months
#     maternal.decay.rt -rate of maternal decay of immunity to rubella
#     exponent - numeric - exponent for the infected
#     frequency.dep - boolean - TRUE
#     is.stochastic - boolean - FALSE
#     tot.pop - numeric - population you want to scale you whole experiment by at the beginning, NULL or numeric
#     yr.births.per.1000 - numeric - crude birth rate per year per 1000
#     coverage - proportion - vaccination coverage
#     intro.rate - numeric - 
#     flat.WIAFW - boolean - if F then default waifw is polymod GB
#     asdr.object - nMx object with rate and mid-age (age specific death rates)
#     targeted.intro - boolean - FALSE means space introduction out over all age classes, TRUE means to concentrate introductions 
#     year - numeric - year to pull country DFE demography (age structure and population size)
#     boulianne.vsucc - boolean - TRUE means use boulianne for MR VE, FALSE means use gans for M VE
#
#Returns -
#     starting state and transition object
#was Get.Country.Tran.Objects
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
  waifw = get.polymod.WAIFW(age.classes/12)
  if (flat.WAIFW) waifw <- get.flat.WAIFW(age.classes/12)
  
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

#This function is an argument of the experiment builder
#  to ensure fertility from susceptible mothers is susceptible
#
#Parameters - 
#   state - state vector at this stage in simulation - length(n.epi.class*n.age.class)
#   tran - transition matrix at this stage in simulation
#   fert.curve - age specific fertlity rate used to establish fraction susceptible mothers
#   lower.age.boundary - corresponding lower age boundaries for the fert.curve
#
#Returns - 
#     vector of births of length(state) 
get.births.fun <- function(state,tran,
                           fert.curve=c(0.0,36.3,210.6,164.8,68.3,27.5,10.1,2.8,0.0),
                           lower.age.boundary = c(0,15,20,25,30,35,40,45,50,101)){ #xxamy - age.class change from lower.age.boundary = c(0,15,20,25,30,35,40,45,50,100)){
  
  # Births based on transition object which is changing at each time point
  births <- tran@birth.rate
  #print(births)
  
  # Establish fertility curve
  n.age.cats <- length(tran@age.class)
  ages.to.use <- (tran@age.class + c(0, tran@age.class[2:n.age.cats-1]))/2
  fertility.category <- findInterval(ages.to.use/12,lower.age.boundary) #tells which fert.cat to pull from for each age cat for all 280 age cats.  so 1 through 180 (12 months * 15 years) is 1 or fertilty rate of 0
  
  # Get the proportion of susceptible mothers
  #the age specific fertility rate is helping with is the determining the proportion of births that are susceptible
  #otherwise the crude birth rate is used to determine number of births
  totpop <- state[tran@m.inds]+state[tran@s.inds]+state[tran@i.inds]+state[tran@r.inds]+state[tran@v.inds] #total pop per age category (280 age cats)
  prop.susc.mothers <- sum(state[tran@s.inds]*fert.curve[fertility.category])/sum(totpop*fert.curve[fertility.category]) #proportion from 0 to 1
  #print(prop.susc.mothers)
  
  # Make the new baby vector - rescaling birth rate for change in births over time
  new.babies <- rep(0,length(state))
  new.babies[1] <- (1-prop.susc.mothers)*births #just dividing the births into S class and M class
  new.babies[2] <- prop.susc.mothers*births
  
  #print(new.babies[1:10])
  return(new.babies)
  
}


#Function to create a state matrix at the disease free equilibrium (DFE)
#from a transition object 
#
#Parameters -
#     uncode - UN country code
#     tot.pop - population by which you can rescale if you want
#     tran - the transition object
#     year - numeric - options are 1950 to 2100
#
#Return -
#     a starting DFE vector of S and R individuals
create.country.x.DFE.ID.state.matrix <- function(uncode, tot.pop=NULL, tran,
                                                 epi.class.label = c("M","S","I","R","V"),
                                                 year=1990, use_montagu_demog=T, ...){
  
  #make template
  rc <- create.ID.state.matrix(tran@n.age.class,
                               tran@n.epi.class,  epi.class.label = epi.class.label)
  
  if (use_montagu_demog) demog <- getDemography(uncode=uncode)
  if (!use_montagu_demog) demog <- getDemography.wpp2019(uncode=uncode)
  pop.struct <- demog$pop.age.byageclasses.1950.2100[,(year-1950+1)]*1000
  age <- as.numeric(rownames(demog$pop.age.byageclasses.1950.2100))
  
  if (is.null(tot.pop)) {
    tot.pop <- demog$pop.total.1950.2100[(year-1950+1)]*1000
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


# Function to produce age survival matrix - xxawinter
#
# Parameters
#   tmp.trans - transitions object of class ID.transition.MSIRV
#   maternal.obj - maternal object of class maternal.exp.decay
#   surv.at.timestep.t -  the estimated age-specific survival rates at time t
#
# outputs the age survival matrix per time t
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


#Function to pull Age Specific Death Rates and return survivorship over time
#
#Paramaters -
#   age.classes - a set of age classes for which tran is being built
#   generation time - generation time
#   nMx - object of age specific death rates over time and mid-age per rate
#   nMx.years - year to start pulling age specific death rates
#   no.gens.in.year - number of generations in a year
#
#Returns-
#   returns the age profile of survivorship in units of the generation time
create.surv.prob.over.age.time <- function(age.classes, generation.time, nMx=NULL, nMx.years, check=F){
  
  mid.age <- nMx@mid.age
  generation.time.year <- generation.time/12   # put generation.time in terms of years rather than months (0.5 in months) now (0.041667 in years)
  no.gens.in.year <- 12*1/generation.time
  
  rate.years <- seq(1950,2100,5)
  indexes <- findInterval(nMx.years, rate.years)
  rates.years <- data.frame(nMx@rates[,indexes])
  rates.generations <- rates.years[,rep(seq_len(ncol(rates.years)), each=no.gens.in.year)]
  rates.generations <- cbind(rates.generations[,1], rates.generations) #adding additional column to match the length of the experiment time steps
  
  survs <- matrix(NA, nrow=length(age.classes), ncol=ncol(rates.generations))
  
  for (t in 1:ncol(rates.generations)){
    yearly.mortality <- rates.generations[,t]
    
    # Fit a smooth spline to this
    fit <- smooth.spline(mid.age, log(yearly.mortality), df=length(mid.age))
    
    # Convert the fit to our ages
    survs[,t] <- 1-exp(predict(fit, age.classes/12)$y) 
    # gives me only the $y or estiamted m.x.n for age.classes/12 (puts age.classes in years rather than months)
    # for each age predict the survival rate based on the smooth.spline
    # exp to get back to m.x.n from it's previous log state and take 1-m.x.n to get survival rate
    #plot(survs)
    
    # Convert yearly rates to generation time rates
    survs[,t] <- 1+(log(survs[,t])/(1/generation.time.year)) 
    
    # Need to adjust last age class to deplete
    survs[,t][length(survs[,t])] <-  0.5^no.gens.in.year
  }
  
  if (check) {
    plot(mid.age, yearly.mortality, type="b", xlim=c(0,50) ,ylim=range(yearly.mortality[mid.age<60]),pch=19)
    points(age.classes/12,(1-survs[,t]^(1/generation.time.year)), col=2,type="b")
    points(fit$x,exp(fit$y),col=4, type="l")
  }
  return(survs)
  
}

### Function to get population size, age structure, and fertility from 1950 to 2105 
### the last five years I just keep all 2100 rates constant
### based on unpd wpp2017 package estimates
###
##' @param - uncode - UN country code
##' @param - age.classes  - in years
##' @param - if.missing.use.region - boolean (MHL, XK, TUV)
###
### Returns demography for each country
### Use to be getDemography.original
getDemography.wpp2017 <- function(uncode=NA, age.classes=c(1:101), if.missing.use.region=F){
  
  time.by5 <- seq(1950,2100,5)
  time <- seq(1950,2100,1)
  
  #UNPD data
  data(pop)
  
  if (any(pop$country_code==uncode)){
    cc <- uncode 
  } else if (if.missing.use.region) {
    if (uncode == 999) {
      region <- "Europe"
      cc <- pop$country_code[tolower(pop$name)==tolower(region)]
      print(paste("No UNPD data for uncode ", uncode, "; using region ", region, " instead", sep=""))
    }
    if (uncode != 999) {
      region <- countrycode::countrycode(uncode, origin="un", destination="continent")
      cc <- pop$country_code[tolower(pop$name)==tolower(region)]
      print(paste("No UNPD data for uncode ", uncode, "; using region ", region, " instead", sep=""))
    }
  } else {
    print(paste("No UNPD data for uncode ", uncode, "; NAs produced", sep=""))
    pop.total.1950.2100=NA
    pop.age.byageclasses.1950.2100=NA
    tfr.1950.2100=NA
    births.1950.2100=NA
    cbr.1950.2100=NA
    e0.1950.2100=NA
    asfr.1950.2100=NA
    repro.age.sex.dist.1950.2100=NA
    asdr.1950.2100.by5=NA
    asdr.object=NA
  }
  
  #Population total in years 1950 to 2100
  data(pop)
  data(popproj)
  pop.total.1950.2100.by5 <-  as.numeric(cbind(pop[pop$country_code==cc,3:ncol(pop)],
                                               popproj[popproj$country_code==cc,3:ncol(popproj)]))
  f <- smooth.spline(time.by5, pop.total.1950.2100.by5)
  pop.total.1950.2100 <- predict(f,time)$y 
  
  #Population Age Structure by Age over time
  data(popFprojMed)
  data(popMprojMed)
  popF.age.1950.2100.by5 <- cbind(popF[popF$country_code==cc,4:ncol(popF)], popFprojMed[popFprojMed$country_code==cc,4:ncol(popFprojMed)])
  popM.age.1950.2100.by5 <- cbind(popM[popM$country_code==cc,4:ncol(popM)], popMprojMed[popMprojMed$country_code==cc,4:ncol(popMprojMed)])
  pop.age.1950.2100.by5 <- popF.age.1950.2100.by5+popM.age.1950.2100.by5
  ages.by5 <- popF[popF$country_code==cc,3]
  mid.age.by5 <- seq(2.5, 102.5, 5)
  popF.age.1950.2100 <- popM.age.1950.2100 <- matrix(NA, length(ages.by5), length(time))
  for (age in 1:length(ages.by5)){
    #new with wpp2017 because there are NA in older ages and oldest cohorts
    time.by5.index <- which(!is.na(popF.age.1950.2100.by5[age,]))
    time.by5.tmp <- time.by5[time.by5.index] 
    time.index <- ((time.by5.index[1]-1)*5+1):length(time)
    time.tmp <- time[time.index]
    f <- smooth.spline(time.by5.tmp, popF.age.1950.2100.by5[age,time.by5.index])
    popF.age.1950.2100[age,time.index] <- predict(f,time.tmp)$y 
    m <- smooth.spline(time.by5.tmp, popM.age.1950.2100.by5[age,time.by5.index])
    popM.age.1950.2100[age,time.index] <- predict(m,time.tmp)$y 
  }
  colnames(popF.age.1950.2100) <- colnames(popM.age.1950.2100) <- time
  rownames(popF.age.1950.2100) <- rownames(popM.age.1950.2100) <- ages.by5
  #forcing 100+ to be 0 if NA - new with wpp2017 - then using na.spline in the zoo package to fill them in
  popF.age.1950.2100[21,is.na(popF.age.1950.2100[21,])] <- 0
  popM.age.1950.2100[21,is.na(popM.age.1950.2100[21,])] <- 0
  popF.age.byageclasses.1950.2100 <- apply(popF.age.1950.2100, 2, function(x) predict(smooth.spline(mid.age.by5, na.spline(x)), age.classes)$y)
  popM.age.byageclasses.1950.2100 <- apply(popM.age.1950.2100, 2, function(x) predict(smooth.spline(mid.age.by5, na.spline(x)), age.classes)$y)
  #adjust for varying bin width!
  popF.age.byageclasses.1950.2100 <- popF.age.byageclasses.1950.2100*diff(c(0,age.classes)) 
  popM.age.byageclasses.1950.2100 <- popM.age.byageclasses.1950.2100*diff(c(0,age.classes))
  pop.tmp <- popF.age.byageclasses.1950.2100+popM.age.byageclasses.1950.2100
  #rescale the age distibution to the estimated population each year (pop.total.1950.2100)
  pop.age.byageclasses.1950.2100 <- matrix(NA, nrow(pop.tmp), ncol(pop.tmp))
  for (t in 1:ncol(pop.tmp)){
    pop.age.byageclasses.1950.2100[,t] <- pop.total.1950.2100[t]*pop.tmp[,t]/sum(pop.tmp[,t])
  }
  rownames(pop.age.byageclasses.1950.2100) <- age.classes
  colnames(pop.age.byageclasses.1950.2100) <- time
  
  #Sex Distribution by age and over time
  pop.tmp <- popF.age.byageclasses.1950.2100+popM.age.byageclasses.1950.2100
  sex.dist <- matrix(NA, nrow(pop.tmp), ncol(pop.tmp))
  for (c in 1:ncol(pop.tmp)){
    sex.dist[,c] <- popF.age.byageclasses.1950.2100[,c]/pop.tmp[,c]
  }
  #Sex Distribution grouped by 5 year reproductive ages over time
  repro.ages.index <- findInterval(age.classes, seq(15,50,5))+1
  repro.age.sex.dist.1950.2100 <- matrix(NA, nrow=length(unique(repro.ages.index)), ncol=ncol(sex.dist))
  for (i in 1:length(unique(repro.ages.index))){
    repro.age.sex.dist.1950.2100[i,] <- apply(sex.dist, 2, function(x) mean(x[repro.ages.index==i]))
  }
  repro.age.sex.dist.1950.2100 <- repro.age.sex.dist.1950.2100[-c(1,9),]
  rownames(repro.age.sex.dist.1950.2100) <- c("15-19","20-24","25-29","30-34","35-39","40-44","45-49")
  colnames(repro.age.sex.dist.1950.2100) <- time
  
  #TFR over time
  data(tfr)
  data(tfrprojMed)
  mid.time.by5 <- seq(1952, 2097, 5) #TFR is given over a range, therefore assume it is the mid-period
  tfr.1950.2100.by5 <-  as.numeric(cbind(tfr[tfr$country_code==cc,3:15],
                                         tfrprojMed[tfrprojMed$country_code==cc,3:ncol(tfrprojMed)]))
  f <- smooth.spline(mid.time.by5, tfr.1950.2100.by5)
  tfr.1950.2100 <- predict(f,time)$y 
  names(tfr.1950.2100) <- time
  
  #Number of women of reproductive age (in five year age groups 15 to 50) over time
  repro.ages <- seq(15,45,5)
  popF.15to50.1950.2100 <- matrix(NA,7,length(1950:2100))
  for (c in 1:ncol(pop.age.byageclasses.1950.2100)){
    for (a in 1:length(repro.ages)){
      index <- which(age.classes>=repro.ages[a] & age.classes<(repro.ages[a]+5))
      popF.15to50.1950.2100[a,c] <- sum(pop.age.byageclasses.1950.2100[index,c]*sex.dist[index,c])
    }
  }
  colnames(popF.15to50.1950.2100) <- time
  rownames(popF.15to50.1950.2100) <- repro.ages
  
  #ASFR overtime
  data(percentASFR)
  p.asfr <- (percentASFR[percentASFR$country_code==cc,])
  p.asfr.1950.2100 <- matrix(NA, nrow=nrow(p.asfr), ncol=(length(1950:2100)))
  for (t in 1:(ncol(p.asfr)-3)){
    p.asfr.1950.2100[,(t*5-4):(t*5)] <- p.asfr[,(t+3)]/100
  }
  p.asfr.1950.2100[,ncol(p.asfr.1950.2100)] <-  p.asfr[,ncol(p.asfr)]/100
  asfr.1950.2100 <- matrix(NA, nrow=nrow(p.asfr), ncol=(length(1950:2100)))
  for (t in 1:ncol(asfr.1950.2100)){
    asfr.1950.2100[,t] <- tfr.1950.2100[t]/5*p.asfr.1950.2100[,t]
  }
  colnames(asfr.1950.2100) <- time
  
  #Births over time
  births.1950.2100 <- rep(NA, length(1950:2100))
  for (t in 1:length(births.1950.2100)){
    births.1950.2100[t] <- sum(asfr.1950.2100[,t]*popF.15to50.1950.2100[,t])
  }
  names(births.1950.2100) <- time
  
  #Crude Birth Rates over time
  cbr.1950.2100 <- rep(NA, length(1950:2100))
  cbr.1950.2100 <- births.1950.2100/pop.total.1950.2100*1000
  names(cbr.1950.2100) <- time
  
  #Age Specific Death Rate over time
  data(mxF)
  asdr.1950.2100.by5 <- (mxF[mxF$country_code==cc,4:ncol(mxF)])
  asdr.1950.2100.by5 <- cbind(asdr.1950.2100.by5, "2100-2105" = asdr.1950.2100.by5[,ncol(asdr.1950.2100.by5)]) #revised to 2105 xxamy
  asdr.maxage <- 5*(length((mxF[mxF$country_code==cc,3]))-2)
  asdr.object <- new("nMx",
                     rates = asdr.1950.2100.by5,
                     mid.age = c(0.5, seq(2.5, (asdr.maxage+2.5), 5)))
  
  #Life Expectancy over time
  data(e0F)
  data(e0Fproj)
  mid.time.by5 <- seq(1952, 2097, 5) #e0 is given over a range, therefore assume it is the mid-period
  e0.1950.2100.by5 <-  as.numeric(cbind(e0F[e0F$country_code==cc,3:15],
                                        e0Fproj[e0Fproj$country_code==cc,3:ncol(e0Fproj)]))
  f <- smooth.spline(mid.time.by5, e0.1950.2100.by5)
  e0.1950.2100 <- predict(f,time)$y 
  names(e0.1950.2100) <- time
  
  
  return(list(pop.total.1950.2100=pop.total.1950.2100, 
              pop.age.byageclasses.1950.2100=pop.age.byageclasses.1950.2100, 
              tfr.1950.2100=tfr.1950.2100,
              e0.1950.2100=e0.1950.2100,
              asfr.1950.2100=asfr.1950.2100,
              repro.age.sex.dist.1950.2100=repro.age.sex.dist.1950.2100,
              births.1950.2100=births.1950.2100, 
              cbr.1950.2100=cbr.1950.2100,
              asdr.1950.2100.by5=asdr.1950.2100.by5,
              asdr.object=asdr.object)) 
}


#### Manipulate Experiment Output ####

#Function to get number of individuals per age group
#
#Parameters - 
#     trans - transition object
#     state - vector of state of the population at one time point
#
#Retunrs - 
#     vector size n.age.classes with the number of individuals in each age
GetNumber.per.AgeGroup <- function(trans, state, tval=1) {
  epi.class <- trans@epi.class
  pop.per.agegroup <- rep(0, length(trans@age.class))
  for (k in 1:trans@n.epi.class) {
    pop.per.agegroup <- pop.per.agegroup+state[epi.class==k]
  }
  return(pop.per.agegroup)
}


#Function to get age-specific population size at mid-year
#
#Parameters - 
#     res - result data matrix
#     trans - transition object
#
#Retunrs - 
#     matrix of pop size by age (in years) per year
GetPop.per.Year <- function(res, trans, no.gens.in.year){
  
  mid.year.index <- seq(13,ncol(res)-1, no.gens.in.year) #mid-year every year
  pop.year <- matrix(NA, 101, length(mid.year.index))
  for (i in 1:length(mid.year.index)){
    t <- mid.year.index[i]
    state <- res[,t]
    pop.all <- GetNumber.per.AgeGroup(state=state, trans=trans)
    pop.year[,i] <- GetNumber.per.AgeYear(vec=pop.all, age.classes=trans@age.class) 
  }
  
  return(pop.year)
}

#Function to get number of individuals per age in years
#
#Parameters -
#       vec - vector of length(age.classes)
#       age.classes - age.classes
#
#Return - 
#     number of individuals per age in years
GetNumber.per.AgeYear <- function(vec, age.classes) {
  
  upper.age.year <- (age.classes[length(age.classes)]-12)/12  #xxamy - age.class change from upper.age.year <- (age.classes[length(age.classes)]-1)/12 #using only to age 58, because by 59 all mostly dead
  diff <- diff(age.classes)
  top.one.month.age <- age.classes[which(diff==unique(diff)[2])[1]]
  top.one.month.age.inyears <- (top.one.month.age)/12 #xxamy - age.class change from top.one.month.age.inyears <- (top.one.month.age-1)/12
  
  pop.per.youngage.year <- rep(0,top.one.month.age.inyears)
  for (y in 1:top.one.month.age.inyears) {
    for (u in 0:11){
      pop.per.youngage.year[y] <- pop.per.youngage.year[y]+vec[(12*y)-u] #xxamy - age.class change from pop.per.youngage.year[y] <- pop.per.youngage.year[y]+vec[(12*y+1)-u] #xxamy - leaves out age group 1
    }
  }
  
  pop.per.age.year <- c(pop.per.youngage.year, vec[(top.one.month.age+1):length(vec)])
  return(pop.per.age.year)
}

#Function to get age-specific rubella cases size per year
#
#Parameters - 
#     res - result data matrix
#     trans - transition object
#     epi.state - epidemiological state indexes
#
#Retunrs - 
#     matrix of pop size by age (in years) per year
GetRubellaCases.per.Year <- function(res, trans, epi.state, no.gens.in.year){
  
  #get number of cases per age in year
  inf <- (res[epi.state==3,])  #cases only
  inf.age <- matrix(NA, 101, ncol(inf))
  for (t in 1:ncol(inf)){
    inf.age[,t] <- GetNumber.per.AgeYear(vec=inf[,t], age.classes=trans@age.class) 
  }
  
  #sum over each year
  inf.pop <- GetSumInYear(inf.age, no.gens.in.year=no.gens.in.year)
  
  return(inf.pop)
}

#Sums over no.gens.in.year to get number per year
#mat - matrix with ncol equal to no.gens.year*year
GetSumInYear <- function(mat, no.gens.in.year){
  max.year <- ncol(mat)/no.gens.in.year
  out <- matrix(0, nrow(mat), max.year)
  for (y in 1:max.year) {
    for (u in 0:(no.gens.in.year-1)){
      #out[,y] <- out[,y]+mat[,(no.gens.in.year*y+1)-u] #leaves out the first time point, ok through b/c only starting population
      out[,y] <- out[,y]+mat[,(no.gens.in.year*y)-u] #CORRECTION NEED TO MAKE, OTHERWISE VACCINE TIMING BEING MOVED UP A YEAR
    }
  }
  return(out)
}

#Function to get age-specific rubella seroprevalence per year
#
#Parameters - 
#     res - result data matrix
#     trans - transition object
#     epi.state - epidemiological state indexes
#
#Retunrs - 
#     matrix of seroprevalence by age (in years) per year
GetRubellaSeroprevalence.per.Year <- function(res, trans, epi.state, no.gens.in.year){
  
  #get number of cases per age in year
  sus <- (res[epi.state==2,])  #susceptibles only
  sus.age <- matrix(NA, floor(max(trans@age.class)/12), ncol(sus))
  for (t in 1:ncol(sus)){
    sus.age[,t] <- GetNumber.per.AgeYear(vec=sus[,t], age.classes=trans@age.class) 
    sus.age[1,t] <- sum(sus[(9:12),t])
  }
  
  #sum over each year
  sus.pop <- sus.age[,seq(12,ncol(res),no.gens.in.year)] #number of susceptibles at mid-year
  pop <- GetPop.per.Year(res, trans, no.gens.in.year) #population at mid-year
  
  seroprev.age <- 1-(sus.pop/pop)
  seroprev.age[which(pop==0)] <- 1 #replace seroprev with 1 if no population in that age group
  
  return(seroprev.age)
}

#Function to get the total number of CRS pregnancies by age
#
#Parameters - 
#       a sim.result.MSIRV object
#       crs_rate_gestational_age_nbiweeks - vector - gestational ages (in biweeks) that coincides with crs_rate
#       crs_rate - vector - probability crs case given rubella infection in pregnancy and live birth that coincides with crs_rate_gestational_age_nbiweeks
#       births - vector - number of births per time period taken from the simulation
#       asfr.over.time - matrix age-specific fertility rate in every age class and year (nrow=7, ncol=length of sim in years)
#
#Returns -
#       an object with age classes set by fertility age classes
getCRScases.byage.per.stochastic.rate <- function(sim.res, crs_rate_gestational_age_nbiweeks=c(6,8), 
                                                  crs_rate=c(0.4, 0.15),
                                                  trans, no.gens.in.year, births, 
                                                  asfr.over.time, repro.age.sex.dist.over.time) {
  
  crs.age.time <- matrix(0,nrow=35,ncol=ncol(sim.res@.Data))
  for (a in 1:length(crs_rate)){
    
    if (a==1){
      # Get a vector providing risk for every age class of susceptibles
      CRSrisk <- getCumRiskTrimester(sim.res=sim.res, nbiweeks=crs_rate_gestational_age_nbiweeks[a])
      #plot(colSums(CRSrisk),type="l")
      
      # Get the number of individuals who become sick in nbiweeks
      SickInFollowingTrim <- sim.res[sim.res@s.inds,]*CRSrisk
      #plot(colSums(SickInFollowingTrim), type="l")
      
    } else { #if a>1
      # Get a vector providing risk for every age class of susceptibles
      CRSrisk.y <- getCumRiskTrimester(sim.res=sim.res, nbiweeks=crs_rate_gestational_age_nbiweeks[(a-1)])
      CRSrisk.o <- getCumRiskTrimester(sim.res=sim.res, nbiweeks=crs_rate_gestational_age_nbiweeks[a])
      CRSrisk <- CRSrisk.o-CRSrisk.y
      
      # Get the number of individuals who become sick in nbiweeks
      SickInFollowingTrim <- sim.res[sim.res@s.inds,]*CRSrisk
    }
    
    # Sum the 15-19 year age groups (in months) into one year age groups
    index.15.to.20 <- findInterval((sim.res@age.class-1)/12,c(15,16,17,18,19,20,101))+1
    tmp.15.to.20 <- matrix(NA, nrow=length(unique(index.15.to.20)), ncol=ncol(SickInFollowingTrim))
    for (i in 1:length(unique(index.15.to.20))){
      tmp.15.to.20[i,] <- apply(SickInFollowingTrim, 2, function(x) sum(x[index.15.to.20==i]))
    }
    SickInFollowingTrim.repro.ages <- rbind(tmp.15.to.20[-c(1,7),], SickInFollowingTrim[241:270,])
    
    #getting fertility indexing
    lower.age.boundary <- c(15,20,25,30,35,40,45,50) # xxamy - age.class change from lower.age.boundary <- c(0,15,20,25,30,35,40,45,49,100) # for asfr
    repro.ages <- 15:49
    fertility.category <- findInterval(repro.ages,lower.age.boundary)
    
    # Multiply to get number of pregnancies with CRS
    crs.age.time.tmp <- crs_rate[a]*(1/no.gens.in.year)*asfr.over.time[fertility.category,]*
      SickInFollowingTrim.repro.ages*repro.age.sex.dist.over.time[fertility.category,]
    
    #Adding to the number of CRS cases based on each gestational age group associated with a rate
    crs.age.time <- crs.age.time+crs.age.time.tmp
    
  }
  
  # Adjust the CRS totals given that births based on CBR not ASFR in the simulation
  #getting population size of women of childbearing age
  pop.age <- matrix(NA, (max(trans@age.class))/12, ncol(sim.res)) #xxamy - age.class change from pop.age <- matrix(NA, (max(trans@age.class)-1)/12, ncol(sim.res))
  for (t in 1:ncol(sim.res)){
    state <- sim.res[,t]
    vec0 <- GetNumber.per.AgeGroup(trans, state, tval=1)
    pop.age[,t] <- GetNumber.per.AgeYear(vec=vec0, age.classes=trans@age.class) 
  }
  row.names(pop.age) <- 0:100
  repro.ages.index <- findInterval(0:100, seq(15,50,5))+1
  repro.age.pop <- matrix(NA, nrow=length(unique(repro.ages.index)), ncol=ncol(pop.age))
  for (i in 1:length(unique(repro.ages.index))){
    repro.age.pop[i,] <- apply(pop.age, 2, function(x) sum(x[repro.ages.index==i]))
  }
  repro.age.pop <- repro.age.pop[-c(1,9),] #narrowing down to only ages 15-49
  #getting births based on age-specific fertility rates
  births.asfr <- colSums(repro.age.pop*asfr.over.time/no.gens.in.year*repro.age.sex.dist.over.time)
  #getting scalar (births based on sim divided by births based on asfr)
  scale <- births/births.asfr
  #adjusting crs given scalar
  crs.adj.age.time <- crs.age.time*scale
  
  return(crs.adj.age.time)
}

#Function to get the VIMC output needed from each result object
#
#Parameters - 
#       sim.res - a sim.result.MSIRV object
#       t.max - numeric - the t.max used to create sim.res
#       setup - country specific object used to create sim.res
#       year - numeric - used to create sim.res
#       crs_rate_gestational_age_nbiweeks - vector - gestational ages (in biweeks) that coincides with crs_rate_agespecific
#       crs_rate_agespecific - vector - probability crs case given rubella infection in pregnancy and live birth that coincides with crs_rate_gestational_age_nbiweeks
#       crs_rate_overall - numeric - average CRS rate over age groups per crs_rate_gestational_age_nbiweeks
#       fetal_death_rate - numric - prob fetal death given rubella infection in pregnancy
#       infant_death_rate - numeric - prob infant death given rubella infection in pregnancy
#
#Returns -
#       a large list of matrices by age and time in biweeks (population, rubella cases, rubella seroprevalence, CRS cases, CRS deaths, CRS dalys, CRS rate)
getOutput_201910gavi_v4 <- function(sim.res, t.max, setup, year, 
                                    crs_rate_gestational_age_nbiweeks, 
                                    crs_rate_agespecific, 
                                    crs_rate_overall, 
                                    fetal_death_rate, 
                                    infant_death_rate,
                                    scenario, run.number, R0){
  
  tmpr <- sim.res
  epi.state <- tmpr@experiment.def@trans@epi.class
  no.gens.in.year <- 1/tmpr@experiment.def@step.size
  
  #Population size by age cohort at mid-year (week 13 in a year)
  pop.1yage.1ytime <- GetPop.per.Year(res=tmpr@result@.Data, trans=tmpr@experiment.def@trans, no.gens.in.year)
  
  #Rubella cases by age cohort and year (sum over each)
  rubella.cases.1yage.1ytime <- GetRubellaCases.per.Year(res=tmpr@result@.Data, trans=tmpr@experiment.def@trans, 
                                                         epi.state=epi.state, no.gens.in.year)
  
  #Seroprevalence by age cohort and year (mid-year seroprev)
  rubella.seropos.1yage.1ytime <- GetRubellaSeroprevalence.per.Year(res=tmpr@result@.Data, trans=tmpr@experiment.def@trans, 
                                                                    epi.state=epi.state, no.gens.in.year)
  
  #CRS cases by womens repro age (15-49) and year
  ##setting up inputs to estimates CRS
  asfr.index <- which(colnames(setup$asfr.1950.2100)==as.character(year)):which(colnames(setup$asfr.1950.2100)==as.character(year+t.max-1))
  asfr.over.time <- matrix(apply(setup$asfr.1950.2100[,asfr.index], 1, function(x) rep(x, each=no.gens.in.year)), ncol=(no.gens.in.year*length(asfr.index)), byrow=T)
  asfr.over.time <- cbind(asfr.over.time[,1], asfr.over.time) #repeat first time point
  repro.age.sex.dist.over.time <- matrix(apply(setup$repro.age.sex.dist.1950.2100[,asfr.index], 1, function(x) rep(x, each=no.gens.in.year)), ncol=(no.gens.in.year*length(asfr.index)), byrow=T)
  repro.age.sex.dist.over.time <- cbind(repro.age.sex.dist.over.time[,1], repro.age.sex.dist.over.time) #repeat first row
  ##getting crs cases by biweek
  crs.1yage.biweek <- as.matrix(getCRScases.byage.per.stochastic.rate(sim.res=tmpr@result, 
                                                                      crs_rate_gestational_age_nbiweeks=crs_rate_gestational_age_nbiweeks, 
                                                                      crs_rate=crs_rate_agespecific, 
                                                                      trans=tmpr@experiment.def@trans, no.gens.in.year, 
                                                                      births=tmpr@result@births.each.timestep, asfr.over.time=asfr.over.time, 
                                                                      repro.age.sex.dist.over.time=repro.age.sex.dist.over.time))
  
  ##getting crs cases by year
  crs.1yage.1ytime <- GetSumInYear(mat=crs.1yage.biweek, no.gens.in.year=no.gens.in.year)
  
  #checking for NAs because it means something went wrong - very important because I change all NAs to zero below
  if(sum(is.na(rubella.cases.1yage.1ytime)) > 0 | sum(is.na(rubella.seropos.1yage.1ytime)) > 0 |
     sum(is.na(crs.1yage.1ytime)) > 0) {
    stop(paste("NAs produced by the model \n", 
               setup$iso3code, "no. rubella cases NAs:", sum(is.na(rubella.cases.1yage.1ytime)),
               "\n", setup$iso3code, "no. rubella seroprev NAs", sum(is.na(rubella.seropos.1yage.1ytime)), 
               "\n",setup$iso3code, "no. CRS cases NAs", sum(is.na(crs.1yage.1ytime)),
               "\n", scenario, ", run number:", run.number, ", R0:", R0))
  }
  
  #Averaging over two year time points because biannual dynamics 
  #Rubella Cases
  last.col <- ncol(rubella.cases.1yage.1ytime)
  rubella.cases.1yage.1ytime.orig <- rubella.cases.1yage.1ytime
  rubella.year.totals <- colSums(rubella.cases.1yage.1ytime)
  rubella.smooth.year.totals <- sapply(seq(2,(length(rubella.year.totals))), function(x) mean(rubella.year.totals[x:(x-1)]))
  rubella.cases.1yage.1ytime <- t(apply(rubella.cases.1yage.1ytime[,-1], 1, function(x) x/(rubella.year.totals[-1])*rubella.smooth.year.totals))
  rubella.cases.1yage.1ytime[is.na(rubella.cases.1yage.1ytime)] <- 0
  rubella.cases.1yage.1ytime <- cbind(rubella.cases.1yage.1ytime.orig[,1], rubella.cases.1yage.1ytime) #add back col 1
  #CRS Cases
  crs.1yage.1ytime.orig <- crs.1yage.1ytime
  crs.year.totals <- colSums(crs.1yage.1ytime)
  crs.smooth.year.totals <- sapply(seq(2,(length(crs.year.totals))), function(x) mean(crs.year.totals[x:(x-1)]))
  crs.1yage.1ytime <- t(apply(crs.1yage.1ytime[,-1], 1, function(x) x/(crs.year.totals[-1])*crs.smooth.year.totals))
  crs.1yage.1ytime[is.na(crs.1yage.1ytime)] <- 0
  crs.1yage.1ytime <- cbind(crs.1yage.1ytime.orig[,1], crs.1yage.1ytime) #add back col 1
  #Rubella Seroprevalence 
  rubella.seropos.1yage.1ytime.orig <- rubella.seropos.1yage.1ytime
  rubella.year.totals <- colSums(rubella.seropos.1yage.1ytime)
  rubella.smooth.year.totals <- sapply(seq(2,(length(rubella.year.totals))), function(x) mean(rubella.year.totals[x:(x-1)]))
  rubella.seropos.1yage.1ytime <- t(apply(rubella.seropos.1yage.1ytime[,-1], 1, function(x) x/(rubella.year.totals[-1])*rubella.smooth.year.totals))
  rubella.seropos.1yage.1ytime[is.na(rubella.seropos.1yage.1ytime)] <- 0
  rubella.seropos.1yage.1ytime <- cbind(rubella.seropos.1yage.1ytime.orig[,1], rubella.seropos.1yage.1ytime) #add back col 1
  #IF vaccine introduced - only average pre-vaccine, keep the post-vaccine raw output
  if(length(which(colSums(tmpr@result@.Data[tmpr@result@v.inds,])!=0))!=0) {
    # don't smooth after vaccine introduction
    index.vaccine.intro <- min(which(colSums(tmpr@result@.Data[tmpr@result@v.inds,])!=0))
    year.vaccine.intro <- ceiling(index.vaccine.intro/no.gens.in.year)
    rubella.cases.1yage.1ytime <- cbind(rubella.cases.1yage.1ytime[,1:(year.vaccine.intro-1)], rubella.cases.1yage.1ytime.orig[,year.vaccine.intro:last.col])
    crs.1yage.1ytime <- cbind(crs.1yage.1ytime[,1:(year.vaccine.intro-1)], crs.1yage.1ytime.orig[,year.vaccine.intro:last.col])
    rubella.seropos.1yage.1ytime <- cbind(rubella.seropos.1yage.1ytime[,1:(year.vaccine.intro-1)], rubella.seropos.1yage.1ytime.orig[,year.vaccine.intro:last.col])
  }
  
  #crs incidence rates
  births <- tmpr@result@births.each.timestep
  sidx = seq.int(from=1, to=length(births), by=no.gens.in.year)
  eidx = c((sidx-1)[2:length(sidx)], length(births))
  births.year = sapply(1:length(sidx), function(i) sum(births[sidx[i]:eidx[i]]))[1:(t.max)]
  crs.rate <- 100000*colSums(crs.1yage.1ytime)/births.year
  
  #Deaths from CRS
  fetal.deaths <- crs.1yage.1ytime/crs_rate_overall*fetal_death_rate
  child.deaths <- crs.1yage.1ytime*infant_death_rate
  deaths.1yage.1ytime <- fetal.deaths+child.deaths
  
  #CRS DALYs
  le.data <- c(62, 66, 74, 79)
  gbd.2010 <- c(29.2, 27.8, 22.9, 19.2)
  le.1yrtime <- setup$e0.1950.2100[which(1950:2100==year):length(setup$e0.1950.2100)]
  daly.percase.2010 <- predict(smooth.spline(le.data, gbd.2010), le.1yrtime)
  daly.percase.2010$y[daly.percase.2010$x<le.data[1]] <- gbd.2010[1]
  daly.percase.2010$y[daly.percase.2010$x>le.data[4]] <- gbd.2010[4]
  dalys.1yage.1ytime <- t(apply(crs.1yage.1ytime,1,function(x) x*daly.percase.2010$y[1:(t.max)]))
  
  #Reducing to only data between 2000 and 2100
  pop.1yage.1ytime <- round(pop.1yage.1ytime[,21:t.max],3) #want 2000-2100
  rubella.cases.1yage.1ytime <- round(rubella.cases.1yage.1ytime[,21:t.max],3) #want 2000-2100
  crs.1yage.1ytime <- round(crs.1yage.1ytime[,21:t.max],3) #want 2000-2100
  deaths.1yage.1ytime <- round(deaths.1yage.1ytime[,21:t.max],3) #want 2000-2100
  dalys.1yage.1ytime <- round(dalys.1yage.1ytime[,21:t.max],3) #want 2000-2100
  rubella.seropos.1yage.1ytime <- round(rubella.seropos.1yage.1ytime[,1:41],3) #want 1980-2020
  crs.rate <- round(crs.rate[1:t.max],3)
  
  return(list(pop.1yage.1ytime=pop.1yage.1ytime,
              rubella.cases.1yage.1ytime=rubella.cases.1yage.1ytime,
              crs.1yage.1ytime=crs.1yage.1ytime,
              deaths.1yage.1ytime=deaths.1yage.1ytime,
              dalys.1yage.1ytime=dalys.1yage.1ytime,
              rubella.seropos.1yage.1ytime=rubella.seropos.1yage.1ytime,
              crs.rate=crs.rate))
}







