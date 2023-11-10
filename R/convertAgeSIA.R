#' Function to convert WHOs SIA age range, as text, into something usable
#'
#' @param sia.age.range input as AgeGroup from the WHO SIA Data
#' @param check boolean of whether to print diagnostics of the functions success, default is FALSE
#' @return Returns age range in months
#'
#' @export
#'

convertAgeSIA <- function(sia.age.range, check=FALSE){

  age.lower <- rep(NA,length(sia.age.range))
  age.lower.unit <- rep(NA,length(sia.age.range))
  age.upper <- rep(NA,length(sia.age.range))
  age.upper.unit <- rep(NA,length(sia.age.range))

  sia.age.range[sia.age.range=="unknown"] <- "9-59 M" # assume the unknowns are the usual 9mo-5yrs
  sia.age.range[sia.age.range=="Children at elementary"] <- "5-11 Y"
  sia.age.range[sia.age.range=="School-age"] <- "5-11 Y"
  sia.age.range[sia.age.range=="eligible children"] <- "9-59 M"
  sia.age.range[sia.age.range=="6 M+"] <- "6 M-15 Y"
  sia.age.range[sia.age.range=="<15 Y"] <- "9 M-14 Y"
  sia.age.range[sia.age.range=="<5 Y"] <- "9-59 M"
  sia.age.range[sia.age.range=="<4 Y"] <- "9-47 M"
  sia.age.range[sia.age.range=="5 Y"] <- "9-59 M"
  sia.age.range[sia.age.range=="14 Y"] <- "9 M-14 Y"
  sia.age.range[sia.age.range=="6 Y"] <- "9 M-6 Y"
  sia.age.range[sia.age.range=="17 Y"] <- "9 M-17 Y"
  sia.age.range[sia.age.range=="18 Y"] <- "9 M-18 Y"
  sia.age.range[sia.age.range=="13 Y"] <- "9 M-13 Y"
  sia.age.range[sia.age.range=="7 Y"] <- "9 M-7 Y"
  sia.age.range[sia.age.range=="22 Y"] <- "9 M-22 Y"
  sia.age.range[sia.age.range=="<1 Y"] <- "6-12 M"
  sia.age.range[sia.age.range=="12 M"] <- "6-12 M"
  sia.age.range[sia.age.range=="1 Y"] <- "6-12 M"
  sia.age.range[sia.age.range=="4 Y"] <- "9-47 M"
  sia.age.range[sia.age.range=="9-5 Y"] <- "9-59 M"
  sia.age.range[sia.age.range=="1 Y school"] <- "5-6 Y"
  sia.age.range[sia.age.range=="1st year primary school"] <- "5-6 Y"
  sia.age.range[sia.age.range=="2-7 or 1-6 Y"] <- "9-83 M"
  sia.age.range[sia.age.range=="+38 Y+"] <- "38-99 Y"
  sia.age.range[sia.age.range=="1 -<5 Y"] <- "9-59 M"
  sia.age.range[sia.age.range=="12 M-45 Y+"] <- "12 M-45 Y"
  sia.age.range[sia.age.range=="6 M-< 5 Y"] <- "6-59 M"
  sia.age.range[sia.age.range=="6 M-<10 Y"] <- "6-119 M"
  sia.age.range[sia.age.range=="6 M-<15 Y"] <- "6-179 M"
  sia.age.range[sia.age.range=="6 M-<6 Y"] <- "6-71 M"
  sia.age.range[sia.age.range=="6 M-<7 Y"] <- "6-84 M"
  sia.age.range[sia.age.range=="6 M->15 Y"] <- "6 M-15 Y"
  sia.age.range[sia.age.range=="6/12-50 Y"] <- "6 M-50 Y"
  sia.age.range[sia.age.range=="9 M-5 Y 6 M"] <- "9-66 M"
  sia.age.range[sia.age.range=="9 M-<10 Y"] <- "9-119 M"
  sia.age.range[sia.age.range=="9 M-<15 Y"] <- "9-179 M"
  sia.age.range[sia.age.range=="9 M-<5 Y"] <- "9-59 M"
  sia.age.range[sia.age.range=="9-5 Y 6 M"] <- "9-66 M"
  sia.age.range[sia.age.range=="<2 Y"] <- "9-23 M"
  sia.age.range[sia.age.range=="<25 Y"] <- "9-299 M"
  sia.age.range[sia.age.range=="<41 Y"] <- "9-492 M"
  sia.age.range[sia.age.range==">12 Y or outside target vaccinated"] <- "12-99 Y"
  sia.age.range[sia.age.range==">15 Y"] <- "15-99 Y"
  sia.age.range[sia.age.range=="school children 7-12 Y"] <- "7-12 Y"
  sia.age.range[sia.age.range=="school children 1-7 Y"] <- "1-7 Y"
  sia.age.range[sia.age.range=="Adults"] <- "15-99 Y"
  sia.age.range[sia.age.range=="All ages"] <- "9 M-99 Y"
  sia.age.range[sia.age.range=="Grade 9"] <- "14-15 Y"
  sia.age.range[sia.age.range=="children and adults"] <- "9 M-99 Y"
  sia.age.range[sia.age.range=="1-<7 Y"] <- "9-83 M"
  sia.age.range[sia.age.range=="2.5 - 7Y"] <- "30 M-7 Y"
  sia.age.range[sia.age.range=="Unvaccinated adults 18-45 Y"] <- "18-45 Y"
  sia.age.range[sia.age.range=="unvaccinated children 2-17 Y"] <- "2-17 Y"

  #Dropping - nothing I can do with these
  sia.age.range[sia.age.range=="risk groups"] <- NA
  sia.age.range[sia.age.range=="Women and Children"] <- NA
  sia.age.range[sia.age.range=="Staff work in airports, ports and piers"] <- NA
  sia.age.range[sia.age.range=="Aborigines in Sg. Berua Village, Hulu Terengganu district, Terengganu"] <- NA
  sia.age.range[sia.age.range=="port heath workers"] <- NA
  sia.age.range[sia.age.range=="Health Care Workers"] <- NA
  sia.age.range[sia.age.range=="refugees"] <- NA
  sia.age.range[sia.age.range=="Military (1980-91 cohorts)"] <- NA
  sia.age.range[sia.age.range=="HCW(1980-91 cohorts)"] <- NA
  sia.age.range[sia.age.range=="Tourism workers"] <- NA
  sia.age.range[sia.age.range=="Airport workers"] <- NA
  sia.age.range[sia.age.range=="CBAW"] <- NA

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

  if(check){
    tmp <- data.frame(age.range)
    tmp$sia.age.range <- sia.age.range
    tmp <- tmp[is.na(tmp$sia.age.range),]
    print(paste(nrow(tmp), "SIAs did not have age specific target information and will not be used"))
    tmp <- tmp[is.na(tmp$age.lower.res) & !is.na(tmp$sia.age.range),]
    print(tmp)
    print("if 0 rows, then no sia.age.range was not accounted for or successfully changed to lower and upper age")
  }

  return(age.range)
}
