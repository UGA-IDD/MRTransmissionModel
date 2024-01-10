#' Make a WAIFW matrix based on Prem et al. 2021 using pakistan for afghanistan
#'
#' @param age.class.boundries the upper age limit for each age class in YEARS
#' @param uncode country code
#' @param other.contact.matrix xxx
#' @param bandwidth desired smooth bandwidth - default=c(3,3)
#' @param adjustment_start_year for 202110gavi_v3 given that Shaun's R0 estimates are adjusted to 1980 pop structure
#' @param year xxx
#'
#' @import countrycode
#'
#' @return a WAIFW matrix based on the Polymod results from chosen location with row and colnames indicating age classes
#' @export
#'

get.prem.WAIFW <- function (age.class.boundries = (1:90),
                            uncode, other.contact.matrix=FALSE,
                            bandwidth=c(3,3), adjustment_start_year=FALSE, year=1980) {

  # Missing Prem Age Contacts
  #these equate to the following iso3:
  # "AND" "ATG" "AUS" "DMA" "GRD" "HTI" "JPN" "KIR" "LBN" "LIE" "MHL" "FSM" "MCO" "NRU" "PLW" "KNA" "SMR" "SYC" "SOM" "TUV", "XK"
  uncode.missing <- c(20, 28, 36, 212, 308, 332, 392, 296, 422, 438, 584, 583, 492, 520, 585, 659, 674, 690, 706, 798, 999)

  # If missing contact matrix - default is GBR
  if (uncode %in% c(20, 28, 36, 212, 308, 392, 422, 438, 583, 492, 520, 585, 659, 674, 690)){
    iso3code <- "GBR"
    print(paste("Prem Contact Matrix Missing for uncode:", uncode, " - GBR set as default"))
  }
  if (uncode==332) iso3code <- "DOM" #HTI - Haiti becomes Dominican Republic
  if (uncode==296) iso3code <- "FJI" #KIR - Kiribati becomes Fuji
  if (uncode==584) iso3code <- "SLB" #MHL - Marshall Islands becomes Solomon Islands
  if (uncode==706) iso3code <- "ETH" #SOM - Somalia becomes Ethiopia
  if (uncode==798) iso3code <- "TON" #TUV - Tuvalu becomes Tonga
  if (uncode==999) iso3code <- "ALB" #XK - Kosovo becomes Albania
  if (!uncode %in% uncode.missing) {
    iso3code <- countrycode::countrycode(uncode, origin="un", destination="iso3c")
  }

  contact_all <- MRTransmissionModel::contact_all
  #contact_all <- readRDS("./data/prem_contact_matrices_2021/contact_all.RDS")

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
