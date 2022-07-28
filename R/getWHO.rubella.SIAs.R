#' Function to convert rubella SIAs into something usable
#'
#' @param iso3code ido3code
#' @param uncode UN country code
#' @param upper.limit upper coverage limit as proportion
#'
#' @importFrom readxl read_excel
#' @importFrom dplyr filter
#'
#' @return age range, year, and coverage of each SIA
#' @export
#'

getWHO.rubella.SIAs <- function(iso3code, uncode, upper.limit=0.97){

  filep <- system.file("/extdata/SIAs_WHO_DownloadedDec2021.xls",
                       package = "MRTransmissionModel")
  df <- read_excel(path = filep, sheet = "SIAs_2000_2022")

  #df <- read_excel("./data/SIAs_WHO_DownloadedDec2021.xls", sheet="SIAs_2000_2022")
  df <- df %>%
    filter(Intervention !="MCV", Intervention !="measles", Intervention !="Measles") %>%
    filter(Activity=="Campaign" | Activity=="CatchUp" | Activity=="FollowUp" | Activity=="SIA") %>%
    filter(Extent!="unknown") %>%
    filter(`Implementation status` =="done") %>%
    filter(!is.na(Year),
           !is.na(`Age group`))

  mtch <- which(df$`ISO3 code`==iso3code,arr.ind=TRUE)
  df <- df[mtch,]

  if (nrow(df)!=0){

    age.range.new <- convertAgeSIA(df$`Age group`)
    years.sia.new <- df$Year

    coverage.new <- rep(NA, length(years.sia.new))
    for (s in 1:length(years.sia.new)){
      age.min.months <- age.range.new[s,1]
      age.min <- ceiling(age.range.new[s,1]/12)
      age.max <- ceiling(age.range.new[s,2]/12)
      y <- years.sia.new[s]

      if (df$Extent[s]=="national" | df$Extent[s]=="National"){
        if (!is.na(df$`Survey results (%)`[s])){
          coverage.new[s] <- (df$`Survey results (%)`[s])/100
        } else {
          coverage.new[s] <- (df$`% Reached`[s])/100
        }
      } else { #Extent is Rollover-National, Sub-national, or Sub-National
        #get national population of the target age
        out <- getDemography.wpp2019(uncode=uncode)
        if (age.min==1){
          tot.pop <- sum(out$pop.age.byageclasses.1950.2100[2:age.max,(y-1950)])*1000
          tot.pop <- as.numeric(tot.pop+(age.min.months/12*out$pop.age.byageclasses.1950.2100[1,(y-1950)]*1000))
        } else {
          tot.pop <- as.numeric(sum(out$pop.age.byageclasses.1950.2100[age.min:age.max,(y-1950)])*1000)
        }
        #multiply coverage by the targeted age then divide by national target age group to get national coverage
        if (!is.na(df$`Survey results (%)`[s])){
          coverage.new[s] <- (df$`Survey results (%)`[s])/100 * (df$`Target population`[s]) / tot.pop
        } else {
          coverage.new[s] <- (df$`Reached population`[s]) / tot.pop
        }
      }
    }

    coverage.new <- (pmin(coverage.new, upper.limit))
    age.range.new.min <- age.range.new[,1]
    age.range.new.max <- age.range.new[,2]

  } else {
    coverage.new <- 0
    age.range.new.min <- 12
    age.range.new.max <- 48
    years.sia.new <- 2000
  }

  return(list(age.min=age.range.new.min, age.max=age.range.new.max, year=years.sia.new, coverage=coverage.new))
}
