#' Function to convert measles SIAs into something usable
#'
#' @param iso3code iso3code
#' @param uncode UN country code
#' @param upper.limit upper coverage limit as proportion
#' @param demog demography from setupCountry.Nov.2023()
#'
#' @importFrom readxl read_excel
#' @import dplyr
#'
#' @return age range, year, and coverage of each SIA
#' @export
#'

getWHO.measles.SIAs.Nov2023 <- function(iso3code, uncode, upper.limit=0.97, demog){

  filep <- system.file("/extdata/SIAs_WHO_Downloaded8Nov2023.xlsx",
                       package = "MRTransmissionModel")
  #filep <- "inst/extdata/SIAs_WHO_Downloaded8Nov2023.xlsx" #use prior to rebuilding package
  df <- readxl::read_excel(path = filep, sheet = "MR_summary_October2023")

  df <- df %>%
    dplyr::filter(Intervention !="Rubella" | Intervention !="mumps" | Intervention !="Mumps") %>%
    dplyr::filter(Activity=="Campaign" | Activity=="CatchUp" | Activity=="FollowUp" | Activity=="SIA" | Activity=="FollowUp - Phased") %>%
    dplyr::filter(!is.na(Extent)) %>%
    dplyr::filter(`Implementation status`=="Done") %>%
    dplyr::filter(!is.na(Year),
                  !is.na(`Age group`))

  mtch <- which(df$`ISO3 code`==iso3code,arr.ind=TRUE)
  df <- df[mtch,]

  if (nrow(df)!=0){

    age.range.new <- convertAgeSIA(sia.age.range = df$`Age group`)
    years.sia.new <- df$Year

    coverage.new <- age.range.new.min <- age.range.new.max <- rep(NA, length(years.sia.new))
    for (s in 1:length(years.sia.new)){
      age.min.months <- age.range.new[s,1]
      age.min <- ceiling(age.range.new[s,1]/12)
      age.max <- ceiling(age.range.new[s,2]/12)
      y <- years.sia.new[s]

      #in convertAgeSIA() NAs could be generated if age groups were not feasible to use - therefore if any NAs then give them 0 coverage
      if (is.na(age.min) | is.na(age.max)){
        coverage.new[s] <- 0
        age.range.new.min[s] <- 12
        age.range.new.max[s] <- 48
      } else {

        if (df$Extent[s]=="national" | df$Extent[s]=="National"){
          if (!is.na(df$`Survey results`[s])){
            coverage.new[s] <- (df$`Survey results`[s])/100
          } else {
            coverage.new[s] <- (df$`% Reached`[s])/100
          }
        } else { #Extent is Rollover-National or Subnational
          #get national population of the target age
          uncode <- countrycode(df$`ISO3 code`[s], "iso3c", "un")
          out <- demog
          if (age.min==1){
            tot.pop <- sum(out$pop.age.byageclasses.1950.2100[2:age.max,(y-1950)])*1000
            tot.pop <- as.numeric(tot.pop+(age.min.months/12*out$pop.age.byageclasses.1950.2100[1,(y-1950)]*1000))
          } else {
            tot.pop <- as.numeric(sum(out$pop.age.byageclasses.1950.2100[age.min:age.max,(y-1950)])*1000)
          }
          #multiply coverage by the targeted age then divide by national target age group to get national coverage
          if (!is.na(df$`Survey results`[s])){
            coverage.new[s] <- (df$`Survey results`[s])/100 * (df$`Target population`[s]) / tot.pop
          } else {
            coverage.new[s] <- (df$`Reached population`[s]) / tot.pop
          }
        }

      coverage.new[s] <- (pmin(coverage.new[s], upper.limit))
      age.range.new.min[s] <- age.range.new[s,1]
      age.range.new.max[s] <- age.range.new[s,2]
      }
    }

    #Remove SIA if we get NA for coverage
    if (any(is.na(coverage.new))){
      print(paste("SIA Coverage is NA in year:", years.sia.new[which(is.na(coverage.new))]))
      coverage.new <- coverage.new[!is.na(coverage.new)]
      age.range.new.min <- age.range.new.min[!is.na(coverage.new)]
      age.range.new.max <- age.range.new.max[!is.na(coverage.new)]
      years.sia.new <- years.sia.new[!is.na(coverage.new)]
    }

    #This combines information for in which multiple SIAs were conducted
    years.sia.new2 <- unique(years.sia.new)
    coverage.new2 <- age.range.new.min2 <- age.range.new.max2 <- rep(NA, length(years.sia.new2))
    for (y in 1:length(years.sia.new2)){
      coverage.new2[y] <-  pmin(sum(coverage.new[years.sia.new %in% years.sia.new2[y]], na.rm=T), upper.limit)
      age.range.new.min2[y] <- age.range.new.min[years.sia.new %in% years.sia.new2[y]][1]
      age.range.new.max2[y] <- age.range.new.max[years.sia.new %in% years.sia.new2[y]][1]
    }

  } else {
    coverage.new2 <- 0
    age.range.new.min2 <- 12
    age.range.new.max2 <- 48
    years.sia.new2 <- 2000
  }

  return(list(age.min=age.range.new.min2,
              age.max=age.range.new.max2,
              year=years.sia.new2,
              coverage=coverage.new2))
}
