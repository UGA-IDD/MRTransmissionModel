
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
library(MRTransmissionModel)
devtools::document()           # regenerate NAMESPACE (run once after @export changes)
devtools::build()              # creates .tar.gz (does NOT install)
source("R/GetPredictedMeaslesSerology.R")
source("R/GetMeaslesSeroprevalence.per.TimePoint.R")
source("R/AggregateMeaslesSeroprevalenceByAgeBins.R")
source("R/LogLikMeaslesSerology.R")
source("R/TransformMeaslesSerologyParameters.R")
source("R/LogLikMeaslesSerology.transformed.R")
source("R/NegLogLikMeaslesSerology.transformed.R")
source("R/FitMeaslesSerology.R")
source("R/EX.Country.part2.R")
#source("R/setClasses.R")
source("R/run.R")

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
  age.classes = c(1:60, seq(72, 1212, 12)),
)


# Full MLE fit — estimates R0 and rho (RI-SIA dose correlation)
# Completed in < 45 minutes. Result: R0 = 18.07, rho = -0.27
format(Sys.time(), "%Y-%m-%d %H:%M:%S")
fit <- FitMeaslesSerology(
  serodata    = serodata,
  setup       = setup,
  year        = 1980,
  t.max       = t.max,
  par.init    = c(log(16), atanh(0), qlogis(0.05)), #c(log(16), atanh(0)),
  age.classes = c(1:60, seq(72, 1212, 12)),
  hessian     = TRUE
)
format(Sys.time(), "%Y-%m-%d %H:%M:%S")

fit$par.natural   # R0 = 18.07, rho = -0.27
head(fit$predictions)

# Normal (Laplace) approximation to the posterior
# fit$cov.transformed is the covariance of (log(R0), atanh(rho)) at the MLE,
# computed by inverting the Hessian of the negative log-likelihood.
# Sampling from MVN(mean=MLE, Sigma=cov.transformed) gives an approximate
# posterior on the transformed scale; back-transforming gives (R0, rho).

n.samples <- 2000
set.seed(42)
samples.transformed <- MASS::mvrnorm(
  n     = n.samples,
  mu    = fit$par.transformed,
  Sigma = fit$cov.transformed
)

# Back-transform to natural scale
posterior.samples <- data.frame(
  R0  = exp(samples.transformed[, "log_R0"]),
  rho = tanh(samples.transformed[, "atanh_rho"])
)

# Posterior summaries: mean, SD, 95% credible interval
posterior.summary <- t(apply(posterior.samples, 2, function(x) c(
  mean    = round(mean(x), 3),
  sd      = round(sd(x), 3),
  lower95 = round(quantile(x, 0.025), 3),
  upper95 = round(quantile(x, 0.975), 3)
)))
print(posterior.summary)

# Marginal posterior plots (red line = MLE point estimate)
par(mfrow = c(1, 2))
hist(posterior.samples$R0,  breaks = 40, main = "Posterior: R0",
     xlab = "R0", col = "lightblue", border = "white")
abline(v = fit$par.natural$R0,  col = "red", lwd = 2)

hist(posterior.samples$rho, breaks = 40, main = "Posterior: rho",
     xlab = "rho", col = "lightblue", border = "white")
abline(v = fit$par.natural$rho, col = "red", lwd = 2)
par(mfrow = c(1, 1))

# Joint posterior scatter plot
plot(posterior.samples$R0, posterior.samples$rho,
     xlab = "R0", ylab = "rho", pch = 16, cex = 0.4, col = "steelblue",
     main = "Joint posterior: R0 vs rho")
points(fit$par.natural$R0, fit$par.natural$rho, col = "red", pch = 4, cex = 2, lwd = 2)


# --- Actual Zambia data ---
load("/Users/winter/Library/CloudStorage/GoogleDrive-amykwinter@gmail.com/My Drive/Indrajit_and_Amy/Seroprevalence data modelling framework/Get_Zambia_Data/all_zambia_data.RData")
serodata <- data.frame(
  survey.time.point = c(rep((2016-1980)*24,nrow(seroprev.2016.raw)), #last time step 2015
                        rep((2024-1980)*24,nrow(seroprev.2024.raw))), #last time step 2023
  age.bin.lower = c(seroprev.2016.raw$age.integer, seroprev.2024.raw$age.integer),
  age.bin.upper = c(seroprev.2016.raw$age.integer, seroprev.2024.raw$age.integer)+1,
  n_tested = c(seroprev.2016.raw$ntest,seroprev.2024.raw$ntest),
  n_positive = c(seroprev.2016.raw$npos,seroprev.2024.raw$npos)
)

LogLikMeaslesSerology(
  R0               = 16,
  rho              = 0,
  serodata         = serodata,
  setup            = setup,
  year             = 1980,
  t.max            = t.max,
  age.classes      = c(1:60, seq(72, 1212, 12)),
  age0is6to11monly = TRUE
)

format(Sys.time(), "%Y-%m-%d %H:%M:%S")
fit <- FitMeaslesSerology(
  serodata         = serodata,
  setup            = setup,
  year             = 1980,
  t.max            = t.max,
  par.init         = c(log(10), qlogis(0.2)),  # c(log(16), atanh(0), qlogis(0.7))
  fix.rho          = 1,
  age.classes      = c(1:60, seq(72, 1212, 12)),
  age0is6to11monly = TRUE,
  generation.time  = 0.5,
  seasonal.amp     = 0.15,
  hessian          = TRUE
)
format(Sys.time(), "%Y-%m-%d %H:%M:%S")
save(fit, file="fit-tmp.RData")
# took 35 minutes with two parameters, R0 = 7.19, rho = 0.9999996
# took 1h 15 minutes with three parameters (), R0 = 6.16, rho = 0.9999959, prop.inacc=0.4011502
# two parameters, R0 = 7.3, scale.sia = 0.18 (hard coded mcv1*0.85 and mcv2*0.90 and rho=1)
head(fit$predictions)
fit$par.natural
head(fit$predictions)
fit$predictions %>%
  filter(survey.time.point==864) %>%
  mutate(seroprev = n_positive/n_tested) %>%
ggplot() +
  geom_line(aes(x=age.bin.lower, y=seroprev)) +
  geom_line(aes(x=age.bin.lower, y=pred.seroprev), color="red")
fit$predictions %>%
  filter(survey.time.point==1056) %>%
  mutate(seroprev = n_positive/n_tested) %>%
  ggplot() +
  geom_line(aes(x=age.bin.lower, y=seroprev)) +
  geom_line(aes(x=age.bin.lower, y=pred.seroprev), color="red")
#FitMeaslesSerology()
#→ optim() calls NegLogLikMeaslesSerology.transformed()
#→ LogLikMeaslesSerology()
#→ GetPredictedMeaslesSerology()
#→ EX.Country.part1() + EX.Country.part2()


# What if i only fit it to 2024 data?
serodata <- data.frame(
  survey.time.point = c(rep((2024-1980)*24,nrow(seroprev.2024.raw))),
  age.bin.lower = c(seroprev.2024.raw$age.integer),
  age.bin.upper = c(seroprev.2024.raw$age.integer)+1,
  n_tested = c(seroprev.2024.raw$ntest),
  n_positive = c(seroprev.2024.raw$npos)
)
format(Sys.time(), "%Y-%m-%d %H:%M:%S")
fit <- FitMeaslesSerology(
  serodata         = serodata,
  setup            = setup,
  year             = 1980,
  t.max            = t.max,
  par.init         = c(log(16),qlogis(0.05)),
  age.classes      = c(1:60, seq(72, 1212, 12)),
  age0is6to11monly = TRUE,
  hessian          = TRUE
)
format(Sys.time(), "%Y-%m-%d %H:%M:%S")


#only fit R0 on older data ----
seroprev.2016.raw %>%
  filter(age.integer>=25) -> seroprev.2016.sub
serodata.subset <- data.frame(
  survey.time.point = c(rep((2016-1980)*24,nrow(seroprev.2016.sub))),
  age.bin.lower = c(seroprev.2016.sub$age.integer),
  age.bin.upper = c(seroprev.2016.sub$age.integer)+1,
  n_tested = c(seroprev.2016.sub$ntest),
  n_positive = c(seroprev.2016.sub$npos)
)
.cache <- new.env(parent = emptyenv())
iter   <- 0

fit.R0 <- optimize(
  f = function(log_R0) {
    iter <<- iter + 1
    R0 <- exp(log_R0)
    cat(sprintf("Iter %d | R0 = %.3f\n", iter, R0))
    -LogLikMeaslesSerology(
      R0               = R0,
      scale.sia        = 0,
      serodata         = serodata.subset,
      setup            = setup,
      year             = 1980,
      t.max            = t.max,
      age.classes      = c(1:60, seq(72, 1212, 12)),
      age0is6to11monly = TRUE,
      .cache           = .cache
    )
  },
  interval = c(log(2), log(50))
)

list(R0 = exp(fit.R0$minimum), logLik = -fit.R0$objective)


pred.full <- GetPredictedMeaslesSerology(
  serodata         = serodata,   # full dataset, both surveys
  setup            = setup,
  year             = 1980,
  t.max            = t.max,
  R0               = exp(fit.R0$minimum),
  scale.sia        = 0,
  age.classes      = c(1:60, seq(72, 1212, 12)),
  age0is6to11monly = TRUE
)

library(ggplot2)
pred.full %>%
  mutate(seroprev = n_positive / n_tested,
         year     = 1980 + survey.time.point / 24) %>%
  ggplot(aes(x = age.bin.lower)) +
  geom_point(aes(y = seroprev)) +
  geom_line(aes(y = pred.seroprev), color = "red") +
  facet_wrap(~year) +
  labs(x = "Age (years)", y = "Seroprevalence")

