#' Function to convert the SIA age range into something usable
#'
#' @param sia.age.range input as AgeGroup from the WHO SIA Data
#'
#' @return Returns age range
#'
#' @export
#'

convertAgeSIA <- function(sia.age.range){

  age.lower <- rep(NA,length(sia.age.range))
  age.lower.unit <- rep(NA,length(sia.age.range))
  age.upper <- rep(NA,length(sia.age.range))
  age.upper.unit <- rep(NA,length(sia.age.range))

  sia.age.range[sia.age.range=="unknown"] <- "9-59 M"               # assume the unkowns are the usual 9mo-5yrs
  sia.age.range[sia.age.range=="Children at elementary"] <- "5-11 Y"
  sia.age.range[sia.age.range=="School-age"] <- "5-11 Y"
  sia.age.range[sia.age.range=="eligible children"] <- "9-59 M"
  sia.age.range[sia.age.range=="6 M+"] <- "6 M-15 Y"
  sia.age.range[sia.age.range=="<15 Y"] <- "9 M-14 Y"
  sia.age.range[sia.age.range=="<5 Y"] <- "9-59 M"
  sia.age.range[sia.age.range=="<4 Y"] <- "9-47 M"
  sia.age.range[sia.age.range=="5 Y"] <- "9-59 M"
  sia.age.range[sia.age.range=="14 Y"] <- "9 M-14 Y"
  sia.age.range[sia.age.range=="13 Y"] <- "9 M-13 Y"
  sia.age.range[sia.age.range=="<1 Y"] <- "6-12 M"
  sia.age.range[sia.age.range=="12 M"] <- "6-12 M"
  sia.age.range[sia.age.range=="1 Y"] <- "6-12 M"
  sia.age.range[sia.age.range=="4 Y"] <- "9-47 M"
  sia.age.range[sia.age.range=="9-5 Y"] <- "9-59 M"
  sia.age.range[sia.age.range=="1 Y school"] <- "5-6 Y"
  sia.age.range[sia.age.range=="1st year primary school"] <- "5-6 Y"

  sia.age.tmp <- rep(NA, length(sia.age.range))
  sia.age.tmp[grepl("-", sia.age.range)] <- as.vector(sia.age.range[grepl("-", sia.age.range)])

  for (i in 1:length(sia.age.tmp)){
    age.mat <- matrix(unlist(strsplit(sia.age.tmp[i], "-")), ncol=2, byrow=TRUE)
    # Define Lower Age
    if (grepl("M", age.mat[,1])){
      age.lower.unit[i] <- "M"
      age.lower[i] <- gsub("M", "", age.mat[,1])
      age.lower[i] <- gsub(" ", "", age.lower[i])
    } else if (grepl("Y", age.mat[,1])){
      age.lower.unit[i] <- "Y"
      age.lower[i] <- gsub("Y", "", age.mat[,1])
      age.lower[i] <- gsub(" ", "", age.lower[i])
    }
    # Define Upper Age
    if (grepl("M", age.mat[,2])){
      age.upper.unit[i] <- "M"
      age.upper[i] <- gsub("M", "", age.mat[,2])
      age.upper[i] <- gsub(" ", "", age.upper[i])
      if (grepl("<", age.upper[i])){
        age.upper[i] <- gsub("<", "", age.upper[i])
        age.upper[i] <- as.numeric(age.upper[i]) - 1
      }
    } else if (grepl("Y", age.mat[,2])){
      age.upper.unit[i] <- "Y"
      age.upper[i] <- gsub("Y", "", age.mat[,2])
      age.upper[i] <- gsub(" ", "", age.upper[i])
      if (grepl("<", age.upper[i])){
        age.upper[i] <- gsub("<", "", age.upper[i])
        age.upper[i] <- as.numeric(age.upper[i]) - 1
      }
    }

    if (grepl("M", age.mat[,1])==FALSE & grepl("Y", age.mat[,1])==FALSE){
      age.lower[i] <- gsub(" ", "", age.mat[,1])
      age.lower.unit[i] <- age.upper.unit[i]
    }
  }

  age.lower.res <- as.numeric(age.lower)
  age.upper.res <- as.numeric(age.upper)
  age.lower.res[age.lower.unit=="Y" & !is.na(age.lower.unit)] <- 12 * age.lower.res[age.lower.unit=="Y" & !is.na(age.lower.unit)]
  age.upper.res[age.upper.unit=="Y" & !is.na(age.upper.unit)] <- 12 * age.upper.res[age.upper.unit=="Y" & !is.na(age.upper.unit)]

  age.range <- as.matrix(cbind(age.lower.res, age.upper.res))

  return(age.range)
}
