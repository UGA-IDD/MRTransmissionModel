#' Function to convert rubella SIAs into something usable
#'
#' @param iso3code ido3code
#' @param uncode UN country code
#' @param upper.limit upper coverage limit as proportion
#'
#' @importFrom readxl read_excel
#' @import dplyr
#'
#' @return age range, year, and coverage of each SIA
#' @export
#'

getWHO.rubella.SIAs.Nov2023 <- function(iso3code, uncode, upper.limit=0.97){

  filep <- system.file("/extdata/SIAs_WHO_Downloaded8Nov2023.xlsx",
                       package = "MRTransmissionModel")
  #filep <- "inst/extdata/SIAs_WHO_Downloaded8Nov2023.xlsx" #use prior to rebuilding package
  df <- readxl::read_excel(path = filep, sheet = "MR_summary_October2023")

  df <- df %>%
    dplyr::filter(Intervention !="Measles" | Intervention !="MEASLES" | Intervention !="mumps" | Intervention !="Mumps") %>%
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

    coverage.new <- rep(NA, length(years.sia.new))
    for (s in 1:length(years.sia.new)){
      print(s)
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
          out <- getDemography.wpp2022(uncode=uncode)
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

        coverage.new <- (pmin(coverage.new, upper.limit))
        age.range.new.min <- age.range.new[,1]
        age.range.new.max <- age.range.new[,2]
      }
    }

  } else {
    coverage.new <- 0
    age.range.new.min <- 12
    age.range.new.max <- 48
    years.sia.new <- 2000
  }

  return(list(age.min=age.range.new.min,
              age.max=age.range.new.max,
              year=years.sia.new,
              coverage=coverage.new))
}
