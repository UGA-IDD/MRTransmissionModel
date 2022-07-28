#' Make a WAIFW matrix based on Polymod (ignoring rescaling by population size that might be sensible - just use raw contacts)
#'
#' @param age.class.boundries the upper age limit for each age class in YEARS
#' @param country levels possible "all", "IT" "DE" "LU" "NL" "PL" "GB" "FI" "BE", default is "GB"
#' @param bandwidth desired smooth bandwidth - default=c(3,3)
#' @param do.touch do you want to only include contacts involving touching? default is FALSE
#'
#' @importFrom KernSmooth bkde2D
#'
#' @return a WAIFW matrix based on the Polymod results from chosen location with row and colnames indicating age classes
#' @export
#'

get.polymod.WAIFW <- function (age.class.boundries = (1:90),
                               country="GB",
                               bandwidth=c(3,3),
                               do.touch=FALSE) {
  #require(KernSmooth)

  #bring in the polymod data on contacts for UK
  polymod <- MRTransmissionModel::polymodRaw
  #polymod <- read.csv("data/polymodRaw.csv")

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
