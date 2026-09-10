#' Function to convert WHOs SIA age range into something usable
#'
#' @param sia.age.range input as AgeGroup from the WHO SIA Data
#' @param check logical of whether to print diagnostics of the functions success, default is FALSE
#'
#' @return Returns age range in months
ConvertAgeSIA_v2 <- function(sia.age.range, check=TRUE){

  age.lower <- rep(NA,length(sia.age.range))
  age.lower.unit <- rep(NA,length(sia.age.range))
  age.upper <- rep(NA,length(sia.age.range))
  age.upper.unit <- rep(NA,length(sia.age.range))

  #normalizing the text to lower case and without spaces
  sia.age.range <- stringr::str_squish(sia.age.range)
  sia.age.range <- stringr::str_to_lower(sia.age.range)

  sia.age.range[sia.age.range=="unknown"] <- "9-59 M" # assume the unknowns are the usual 9mo-5yrs
  sia.age.range[sia.age.range=="9 m+"] <- "9-59 M"
  sia.age.range[sia.age.range=="children at elementary"] <- "5-11 Y"
  sia.age.range[sia.age.range=="school-age"] <- "5-11 Y"
  sia.age.range[sia.age.range=="eligible children"] <- "9-59 M"
  sia.age.range[sia.age.range=="6 m+"] <- "6 M-15 Y"
  sia.age.range[sia.age.range=="<15 y"] <- "9 M-14 Y"
  sia.age.range[sia.age.range=="<5 y"] <- "9-59 M"
  sia.age.range[sia.age.range=="<4 y"] <- "9-47 M"
  sia.age.range[sia.age.range=="5 y"] <- "9-59 M"
  sia.age.range[sia.age.range=="14 y"] <- "9 M-14 Y"
  sia.age.range[sia.age.range=="6 y"] <- "9 M-6 Y"
  sia.age.range[sia.age.range=="17 y"] <- "9 M-17 Y"
  sia.age.range[sia.age.range=="18 y"] <- "9 M-18 Y"
  sia.age.range[sia.age.range=="12 y"] <- "9 M-12 Y"
  sia.age.range[sia.age.range=="13 y"] <- "9 M-13 Y"
  sia.age.range[sia.age.range=="7 y"] <- "9 M-7 Y"
  sia.age.range[sia.age.range=="22 y"] <- "9 M-22 Y"
  sia.age.range[sia.age.range=="<1 y"] <- "6-12 M"
  sia.age.range[sia.age.range=="12 m"] <- "6-12 M"
  sia.age.range[sia.age.range=="18 m"] <- "9-18 M"
  sia.age.range[sia.age.range=="9 m"] <- "6-9 M"
  sia.age.range[sia.age.range=="1 y"] <- "6-12 M"
  sia.age.range[sia.age.range=="4 y"] <- "9-47 M"
  sia.age.range[sia.age.range=="9-5 y"] <- "9-59 M"
  sia.age.range[sia.age.range=="1 y school"] <- "5-6 Y"
  sia.age.range[sia.age.range=="1st year primary school"] <- "5-6 Y"
  sia.age.range[sia.age.range=="2-7 or 1-6 y"] <- "9-83 M"
  sia.age.range[sia.age.range=="+38 y+"] <- "38-99 Y"
  sia.age.range[sia.age.range=="1 -<5 y"] <- "9-59 M"
  sia.age.range[sia.age.range=="12 m-45 y+"] <- "12 M-45 Y"
  sia.age.range[sia.age.range=="6 m-< 5 y"] <- "6-59 M"
  sia.age.range[sia.age.range=="6 m-<10 y"] <- "6-119 M"
  sia.age.range[sia.age.range=="6 m-<15 y"] <- "6-179 M"
  sia.age.range[sia.age.range=="6 m-<6 y"] <- "6-71 M"
  sia.age.range[sia.age.range=="6 m-<7 y"] <- "6-84 M"
  sia.age.range[sia.age.range=="6 m->15 y"] <- "6 M-15 Y"
  sia.age.range[sia.age.range=="6/12-50 y"] <- "6 M-50 Y"
  sia.age.range[sia.age.range=="9 m-5 y 6 m"] <- "9-66 M"
  sia.age.range[sia.age.range=="9 m-<10 y"] <- "9-119 M"
  sia.age.range[sia.age.range=="9 m-<15 y"] <- "9-179 M"
  sia.age.range[sia.age.range=="9 m-<5 y"] <- "9-59 M"
  sia.age.range[sia.age.range=="9-5 y 6 m"] <- "9-66 M"
  sia.age.range[sia.age.range=="<2 y"] <- "9-23 M"
  sia.age.range[sia.age.range=="<7 y"] <- "9-84 M"
  sia.age.range[sia.age.range=="<25 y"] <- "9-299 M"
  sia.age.range[sia.age.range=="<41 y"] <- "9-492 M"
  sia.age.range[sia.age.range=="<17 y"] <- "9-204 M"
  sia.age.range[sia.age.range==">12 y or outside target vaccinated"] <- "12-99 Y"
  sia.age.range[sia.age.range==">15 y"] <- "15-99 Y"
  sia.age.range[sia.age.range=="18 y+"] <- "18-99 Y"
  sia.age.range[sia.age.range=="15 y+"] <- "15-99 Y"
  sia.age.range[sia.age.range=="school children 7-12 y"] <- "7-12 Y"
  sia.age.range[sia.age.range=="school children 1-7 y"] <- "1-7 Y"
  sia.age.range[sia.age.range=="adults"] <- "15-99 Y"
  sia.age.range[sia.age.range=="all ages"] <- "9 M-99 Y"
  sia.age.range[sia.age.range=="grade 9"] <- "14-15 Y"
  sia.age.range[sia.age.range=="children and adults"] <- "9 M-99 Y"
  sia.age.range[sia.age.range=="1-<7 y"] <- "9-83 M"
  sia.age.range[sia.age.range=="2.5 - 7y"] <- "30 M-7 Y"
  sia.age.range[sia.age.range=="unvaccinated adults 18-45 y"] <- "18-45 Y"
  sia.age.range[sia.age.range=="unvaccinated children 2-17 y"] <- "2-17 Y"
  sia.age.range[sia.age.range=="2-7 y or 1-6 y"] <- "1-7 Y"

  #Dropping - nothing I can do with these
  sia.age.range[sia.age.range=="risk groups"] <- NA
  sia.age.range[sia.age.range=="13 m-8 y from vulnerable groups"] <- NA
  sia.age.range[sia.age.range=="women and children"] <- NA
  sia.age.range[sia.age.range=="staff work in airports, ports and piers"] <- NA
  sia.age.range[sia.age.range=="aborigines in sg. berua village, hulu terengganu district, terengganu"] <- NA
  sia.age.range[sia.age.range=="port heath workers"] <- NA
  sia.age.range[sia.age.range=="health care workers"] <- NA
  sia.age.range[sia.age.range=="refugees"] <- NA
  sia.age.range[sia.age.range=="military (1980-91 cohorts)"] <- NA
  sia.age.range[sia.age.range=="hcw(1980-91 cohorts)"] <- NA
  sia.age.range[sia.age.range=="tourism workers"] <- NA
  sia.age.range[sia.age.range=="airport workers"] <- NA
  sia.age.range[sia.age.range=="cbaw"] <- NA
  sia.age.range[sia.age.range=="other"] <- NA
  sia.age.range[stringr::str_detect(sia.age.range,"risk group|rick group|hcw|health care worker|poe worker|traveller|airport worker|tourism worker|military|refugee|port heath|port health")] <- NA

  #turn any uppercase into lower case
  sia.age.range <- stringr::str_to_lower(sia.age.range)

  # prep loop
  sia.age.tmp <- rep(NA,length(sia.age.range))
  has.range <- !is.na(sia.age.range) & grepl("-",sia.age.range)
  sia.age.tmp[has.range] <- sia.age.range[has.range]

  for (i in seq_along(sia.age.tmp)){
    if(is.na(sia.age.tmp[i])) next
    age.mat <- matrix(unlist(strsplit(sia.age.tmp[i], "-")), ncol=2, byrow=TRUE)
    # Define Lower Age
    if (grepl("m", age.mat[,1])){
      age.lower.unit[i] <- "M"
      age.lower[i] <- gsub("m", "", age.mat[,1])
      age.lower[i] <- gsub(" ", "", age.lower[i])
    } else if (grepl("y", age.mat[,1])){
      age.lower.unit[i] <- "Y"
      age.lower[i] <- gsub("y", "", age.mat[,1])
      age.lower[i] <- gsub(" ", "", age.lower[i])
    }
    # Define Upper Age
    if (grepl("m", age.mat[,2])){
      age.upper.unit[i] <- "M"
      age.upper[i] <- gsub("m", "", age.mat[,2])
      age.upper[i] <- gsub(" ", "", age.upper[i])
      if (grepl("<", age.upper[i])){
        age.upper[i] <- gsub("<", "", age.upper[i])
        age.upper[i] <- as.numeric(age.upper[i]) - 1
      }
    } else if (grepl("y", age.mat[,2])){
      age.upper.unit[i] <- "Y"
      age.upper[i] <- gsub("y", "", age.mat[,2])
      age.upper[i] <- gsub(" ", "", age.upper[i])
      if (grepl("<", age.upper[i])){
        age.upper[i] <- gsub("<", "", age.upper[i])
        age.upper[i] <- as.numeric(age.upper[i]) - 1
      }
    }

    if (grepl("m", age.mat[,1])==FALSE & grepl("y", age.mat[,1])==FALSE){
      age.lower[i] <- gsub(" ", "", age.mat[,1])
      age.lower.unit[i] <- age.upper.unit[i]
    }
  }

  age.lower.res <- as.numeric(age.lower)
  age.upper.res <- as.numeric(age.upper)
  age.lower.res[age.lower.unit=="Y" & !is.na(age.lower.unit)] <- 12 * age.lower.res[age.lower.unit=="Y" & !is.na(age.lower.unit)]
  age.upper.res[age.upper.unit=="Y" & !is.na(age.upper.unit)] <- 12 * age.upper.res[age.upper.unit=="Y" & !is.na(age.upper.unit)]

  age.range <- as.matrix(cbind(age.lower.res, age.upper.res))
  colnames(age.range) <- c("age.lower.res","age.upper.res")

  if(check){
    tmp <- data.frame(age.range)
    tmp$sia.age.range <- sia.age.range
    tmp.drop <- tmp[is.na(tmp$sia.age.range),]
    print(paste(nrow(tmp.drop),"SIAs did not have age specific target information and will not be used"))
    tmp.unparsed <- tmp[is.na(tmp$age.lower.res) & !is.na(tmp$sia.age.range),]
    print(tmp.unparsed)
    print("if 0 rows, then all non-dropped sia.age.range values were accounted for")
  }

  return(age.range)
}
