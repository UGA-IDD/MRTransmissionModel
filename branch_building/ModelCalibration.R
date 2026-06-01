
load("/Users/winter/Library/CloudStorage/GoogleDrive-amykwinter@gmail.com/My Drive/Indrajit_and_Amy/Seroprevalence data modelling framework/Get_Zambia_Data/simulated_zambia_data.RData")
serodata <- data.frame(
  survey.time.point = c(rep((2023-1980)*24,length(seroprev.1980to2025$age_years)),
                        rep((2019-1980)*24,length(seroprev.1980to2025$age_years)),
                        rep((2015-1980)*24,length(seroprev.1980to2025$age_years))),
  age.bin.lower = rep(seroprev.1980to2025$age_years,3),
  age.bin.upper = rep(seroprev.1980to2025$age_years+1,3),
  n_tested = rep(100, length(seroprev.1980to2025$year_2023)*3),
  n_positive = c(round(100*seroprev.1980to2025$year_2023),
                 round(100*seroprev.1980to2025$year_2019),
                 round(100*seroprev.1980to2025$year_2015))
)


library(MRTransmissionModel)
setup <- setupCountry.Nov2023(country="Zambia")
year <- 1980
t.max <- 45

source("R/GetPredictedMeaslesSerology.R")
source("R/GetMeaslesSeroprevalence.per.TimePoint.R")
source("R/AggregateMeaslesSeroprevalenceByAgeBins.R")
source("R/LogLikMeaslesSerology.R")

LogLikMeaslesSerology(R0 = 16, sia.scale = 0.8, serodata = serodata,
                      setup = setup, year = 1980, t.max = t.max,
                      age.classes = c(1:60, seq(72, 1212, 12)))

#this took 36 secs

fit <- FitMeaslesSerology(
  serodata = serodata,
  setup = setup,
  year = 1980,
  t.max = t.max,
  par.init = c(log(16), qlogis(0.8)),
  age.classes = c(1:60, seq(72, 1212, 12))
)
format(Sys.time(), "%Y-%m-%d %H:%M:%S")
#this took less than 2 hours

fit$par.natural
head(fit$predictions)
