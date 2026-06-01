
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

# --- Package loading ---
# devtools::load_all(".") crashes R on this machine -- do not use.
# Workflow that works:
#   devtools::document()           # regenerate NAMESPACE (run once after @export changes)
#   devtools::build()              # creates .tar.gz (does NOT install)
#   devtools::install(".", quick=TRUE)  # run when C code changes or to bake new R files into library()
library(MRTransmissionModel)
source("R/GetPredictedMeaslesSerology.R")
source("R/GetMeaslesSeroprevalence.per.TimePoint.R")
source("R/AggregateMeaslesSeroprevalenceByAgeBins.R")
source("R/LogLikMeaslesSerology.R")
source("R/TransformMeaslesSerologyParameters.R")
source("R/LogLikMeaslesSerology.transformed.R")
source("R/NegLogLikMeaslesSerology.transformed.R")
source("R/FitMeaslesSerology.R")

setup <- setupCountry.Nov2023(country="Zambia")
year  <- 1980
t.max <- 45


# Single likelihood evaluation at R0=16, rho=0 (independence)
LogLikMeaslesSerology(
  R0          = 16,
  rho         = 0,
  serodata    = serodata,
  setup       = setup,
  year        = 1980,
  t.max       = t.max,
  age.classes = c(1:60, seq(72, 1212, 12))
)


# Full MLE fit — estimates R0 and rho (RI-SIA dose correlation)
# Completed in < 45 minutes. Result: R0 = 18.07, rho = -0.27
format(Sys.time(), "%Y-%m-%d %H:%M:%S")
fit <- FitMeaslesSerology(
  serodata    = serodata,
  setup       = setup,
  year        = 1980,
  t.max       = t.max,
  par.init    = c(log(16), atanh(0)),
  age.classes = c(1:60, seq(72, 1212, 12))
)
format(Sys.time(), "%Y-%m-%d %H:%M:%S")

fit$par.natural   # R0 = 18.07, rho = -0.27
head(fit$predictions)
