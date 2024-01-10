#' Make a parametric smooth WAIFW from Farrington, J. Amer. Stat. Assoc.; 2005, 100 p370
#'
#' @param age.class.boundries the upper age limit for each age class in years
#' @param mu age of highest contact increases with mu
#' @param sig width around equal age diagonal increases with gam
#' @param gam decreases strength in other diagonal (shrinks high trans int)
#' @param delta background homogeneous contact rate
#'
#' @return a smooth WAIFW matrix
#' @export
#'

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
