#' Function to find the non-scaled MR1 CDF (uses most recent DHS data as of 2019)
#'
#' @param uncode numeric; country UN code
#' @param min.age numeric; age in months
#' @param max.age numeric; age in months
#'
#' @importFrom utils read.csv
#' @import survival
#' @importFrom stats smooth.spline
#'
#' @return CDF vaccinated across age in months
#' @export
#'

get.MR1cdf.survival <- function(uncode, min.age=1, max.age=24){

  data <- MRTransmissionModel::DHS_ageMCV1_Data
  #data <- read.csv("./data/DHS_ageMCV1_Data.csv")

  dat <- data[which(data$uncode==uncode),]
  ages <- min.age:max.age
  if (nrow(dat)>0){ #if DHS data for specific country - use country data
    obj <- survival::Surv(time=dat$age, time2=rep(max(dat$age), length(dat$age)), event=as.numeric(dat$censoring), type="interval")
    out <- summary(survfit(obj~1))
    sm <- smooth.spline(c(0,out$time), c(0,(1-out$surv)))
    cdf <- predict(sm, ages)$y
  } else { #otherwise use region data

    cc <- MRTransmissionModel::country_codes
    #cc <- read.csv("./data/country_codes.csv")

    reg <- cc$Region_Code[which(cc$uncode==uncode)]
    if (reg=="AFRO"){

      dat <- MRTransmissionModel::DHS_ageMCV1_Data_AFRO
      #dat <- read.csv("./data/DHS_ageMCV1_Data_AFRO.csv")

    } else if (reg=="EMRO"){

      dat <- MRTransmissionModel::DHS_ageMCV1_Data_EMRO
      #dat <- read.csv("./data/DHS_ageMCV1_Data_EMRO.csv")

    } else if (reg=="SEARO"){

      dat <- MRTransmissionModel::DHS_ageMCV1_Data_SEARO
      #dat <- read.csv("./data/DHS_ageMCV1_Data_SEARO.csv")

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
