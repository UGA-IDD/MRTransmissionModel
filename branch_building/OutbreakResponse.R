
# ============================================================
# Outbreak Response (OBR) — working example for Zambia
# ============================================================
#
# Compares scenarios using Zambia demography and vaccination history.
# All scenarios use the observational model (Se/Sp/background required).
#
# Trigger: sum of obs.TP (IgM-confirmed measles) in children 0-14 years
# over the trigger window reaches or.n.confirmations.target = 5.
#
# Campaign age targeting: CDF of accumulated case ages at or.vacc.agedist.percentile.
#   stop_testing_at_trigger = FALSE -> CDF built from obs.TP (confirmed cases)
#   stop_testing_at_trigger = TRUE  -> CDF built from all suspected cases
#
# R0 = 18 from MAP calibration to Zambia serology data.
# Step size: 0.5 months (24 steps/year).
# ============================================================

library(MRTransmissionModel)
devtools::document()
devtools::build()

setup <- setupCountry.Nov2023(country = "Zambia")
year  <- 1980
t.max <- 45

# --- Step 1: shared transient spin-up (part1) ---

EXt0 <- EX.Country.part1(
  uncode                         = setup$uncode,
  generation.time                = 0.5,
  age.classes                    = c(1:240, seq(252, 1212, 12)),
  maternal.decay.rt              = 0.45,
  exponent                       = 0.97,
  frequency.dep                  = TRUE,
  yr.births.per.1000.acrossyears = rep(setup$cbr.1950.2100[year - 1950 + 1], 20),
  intro.rate                     = 1 / 24 / 320,
  tot.pop                        = NULL,
  targeted.intro                 = FALSE,
  R0                             = 18,
  t.max                          = 20,
  get.births                     = setup$get.births.here,
  seasonal.amp                   = 0.15,
  flat.WAIFW                     = FALSE,
  country.specific.WAIFW         = TRUE,
  vynnycky.waifw                 = FALSE,
  vynnycky.waifw.betas           = c(2, 1),
  asdr.object                    = setup$asdr.object,
  year                           = year,
  use_montagu_demog              = FALSE
)

# --- Shared vaccination inputs ---

mr1cov <- setup$MCV1.coverage.1980to2100[(year - 1980 + 1):t.max] * 0.5
mr2cov <- setup$MCV2.coverage.1980to2100[(year - 1980 + 1):t.max] * 0.5
siacov <- setup$measlesSIA.coverage.1980to2100[(year - 1980 + 1):t.max] * 0.5

# --- Shared observational model inputs ---
# Background rash: exponentially declining with age, ~20% higher in rainy season (Nov-Apr).
# or.reporting.rate = 0.10: ~10% of true infections present as suspected cases (pre-testing).

age.classes.ex <- c(1:240, seq(252, 1212, 12))
n.age.ex       <- length(age.classes.ex)

age.bg.rate <- 0.4 * exp(-age.classes.ex / 48)   # cases per step per age class
month.mult  <- c(1.15, 1.15, 1.10, 1.05, 1.00, 0.90,   # Jan-Jun
                 0.85, 0.85, 0.90, 0.95, 1.05, 1.15)   # Jul-Dec
non.meas.bg <- outer(age.bg.rate, month.mult)            # n.age x 12

shared.obs <- list(
  or.reporting.rate              = 0.10,
  or.non.meas.cases.by.age.month = non.meas.bg,
  or.Se                          = 0.95,
  or.Sp                          = 0.97
)


# --- Scenario A: Routine + SIA only (no OBR) ---
# or.vacc.coverage = 0 disables the response campaign.

result.no.obr <- EX.Country.part2.OR(
  uncode                         = setup$uncode,
  generation.time                = 0.5,
  pop.rescale                    = setup$pop.total.1950.2100[
                                     seq(1990, (year + t.max - 10), 10) - 1950 + 1],
  pop.time                       = seq(1990, (year + t.max - 10), 10) - year,
  is.stochastic                  = FALSE,
  t.max                          = t.max,
  rescale.WAIFW                  = FALSE,
  yr.births.per.1000.acrossyears = setup$cbr.1950.2100[
                                     (year - 1950 + 1):((year - 1950 + 1) + t.max)],
  asdr.object                    = setup$asdr.object,
  year                           = year,
  EXt0                           = EXt0,
  time.specific.MR1cov           = mr1cov,
  time.specific.MR2cov           = mr2cov,
  time.specific.SIAcov           = siacov,
  time.specific.min.age.MR1      = rep(9,  t.max),
  time.specific.max.age.MR1      = rep(24, t.max),
  time.specific.min.age.MR2      = rep(25, t.max),
  time.specific.max.age.MR2      = rep(36, t.max),
  time.specific.min.age.SIA      = setup$age.min.sia.measles[(year - 1980 + 1):t.max],
  time.specific.max.age.SIA      = setup$age.max.sia.measles[(year - 1980 + 1):t.max],
  obj.vcdf.MR1                   = get.MR1cdf.survival(uncode = setup$uncode, 1, 24),
  obj.vcdf.MR2                   = get.vcdf.normal(25, 36),
  obj.prob.vsucc                 = pvacsuccess(c(1:240, seq(252, 1212, 12)),
                                               get.boulianne.vsucc()),
  sia.timing.in.year             = 3 / 12,
  MR1MR2correlation              = TRUE,
  MR1SIAcorrelation              = FALSE,
  MR2SIAcorrelation              = FALSE,
  intro.rate                     = 1 / 24 / 320,
  # --- Observational model ---
  or.reporting.rate              = shared.obs$or.reporting.rate,
  or.non.meas.cases.by.age.month = shared.obs$or.non.meas.cases.by.age.month,
  or.Se                          = shared.obs$or.Se,
  or.Sp                          = shared.obs$or.Sp,
  # --- OBR trigger ---
  or.n.confirmations.target      = 5,
  or.confirmation.delay          = 2,          # 1-month confirmation lag
  or.response.delay              = 2,          # 1-month mobilization lag
  or.trigger.window              = 4,          # 2-month lookback window
  or.trigger.window.age.lower    = 0,          # trigger: all ages 0-14 years
  or.trigger.window.age.upper    = 168,
  # --- OBR response (disabled) ---
  or.vacc.coverage               = 0,          # DISABLED
  or.vacc.agedist.percentile     = 0.95,
  stop_testing_at_trigger        = FALSE,
  or.min.interval                = 48,
  or.start.timestep              = (2020 - year) * 24 + 1
)


# --- Scenario B: Routine + SIA + RDT OBR (age targets from confirmed cases) ---
# stop_testing_at_trigger = FALSE: age targeting uses obs.TP accumulated
# through the response window (diagnostic confirmation continues).
# or.confirmation.delay = 0

result.obr.tp <- EX.Country.part2.OR(
  uncode                         = setup$uncode,
  generation.time                = 0.5,
  pop.rescale                    = setup$pop.total.1950.2100[
                                     seq(1990, (year + t.max - 10), 10) - 1950 + 1],
  pop.time                       = seq(1990, (year + t.max - 10), 10) - year,
  is.stochastic                  = FALSE,
  t.max                          = t.max,
  rescale.WAIFW                  = FALSE,
  yr.births.per.1000.acrossyears = setup$cbr.1950.2100[
                                     (year - 1950 + 1):((year - 1950 + 1) + t.max)],
  asdr.object                    = setup$asdr.object,
  year                           = year,
  EXt0                           = EXt0,
  time.specific.MR1cov           = mr1cov,
  time.specific.MR2cov           = mr2cov,
  time.specific.SIAcov           = siacov,
  time.specific.min.age.MR1      = rep(9,  t.max),
  time.specific.max.age.MR1      = rep(24, t.max),
  time.specific.min.age.MR2      = rep(25, t.max),
  time.specific.max.age.MR2      = rep(36, t.max),
  time.specific.min.age.SIA      = setup$age.min.sia.measles[(year - 1980 + 1):t.max],
  time.specific.max.age.SIA      = setup$age.max.sia.measles[(year - 1980 + 1):t.max],
  obj.vcdf.MR1                   = get.MR1cdf.survival(uncode = setup$uncode, 1, 24),
  obj.vcdf.MR2                   = get.vcdf.normal(25, 36),
  obj.prob.vsucc                 = pvacsuccess(c(1:240, seq(252, 1212, 12)),
                                               get.boulianne.vsucc()),
  sia.timing.in.year             = 3 / 12,
  MR1MR2correlation              = TRUE,
  MR1SIAcorrelation              = FALSE,
  MR2SIAcorrelation              = FALSE,
  intro.rate                     = 1 / 24 / 320,
  # --- Observational model ---
  or.reporting.rate              = shared.obs$or.reporting.rate,
  or.non.meas.cases.by.age.month = shared.obs$or.non.meas.cases.by.age.month,
  or.Se                          = shared.obs$or.Se,
  or.Sp                          = shared.obs$or.Sp,
  # --- OBR trigger ---
  or.n.confirmations.target      = 5,
  or.confirmation.delay          = 1,
  or.response.delay              = 6,
  or.trigger.window              = 4,
  or.trigger.window.age.lower    = 0,
  or.trigger.window.age.upper    = 168,
  # --- OBR response (age targets from TP) ---
  or.vacc.coverage               = 0.85,
  or.vacc.agedist.percentile     = 0.95,       # vaccinate up to 95th pctile of case age CDF
  stop_testing_at_trigger        = FALSE,       # use confirmed cases (obs.TP) for age targeting
  or.min.interval                = 48,
  or.start.timestep              = (2020 - year) * 24 + 1
)


# --- Scenario C: Routine + SIA + EIA OBR (age targets from suspected cases) ---
# stop_testing_at_trigger = TRUE: age targeting uses all suspected cases
# (true.reported + non-measles background; reporting process only).
# or.confirmation.delay = 2

result.obr.susp <- EX.Country.part2.OR(
  uncode                         = setup$uncode,
  generation.time                = 0.5,
  pop.rescale                    = setup$pop.total.1950.2100[
                                     seq(1990, (year + t.max - 10), 10) - 1950 + 1],
  pop.time                       = seq(1990, (year + t.max - 10), 10) - year,
  is.stochastic                  = FALSE,
  t.max                          = t.max,
  rescale.WAIFW                  = FALSE,
  yr.births.per.1000.acrossyears = setup$cbr.1950.2100[
                                     (year - 1950 + 1):((year - 1950 + 1) + t.max)],
  asdr.object                    = setup$asdr.object,
  year                           = year,
  EXt0                           = EXt0,
  time.specific.MR1cov           = mr1cov,
  time.specific.MR2cov           = mr2cov,
  time.specific.SIAcov           = siacov,
  time.specific.min.age.MR1      = rep(9,  t.max),
  time.specific.max.age.MR1      = rep(24, t.max),
  time.specific.min.age.MR2      = rep(25, t.max),
  time.specific.max.age.MR2      = rep(36, t.max),
  time.specific.min.age.SIA      = setup$age.min.sia.measles[(year - 1980 + 1):t.max],
  time.specific.max.age.SIA      = setup$age.max.sia.measles[(year - 1980 + 1):t.max],
  obj.vcdf.MR1                   = get.MR1cdf.survival(uncode = setup$uncode, 1, 24),
  obj.vcdf.MR2                   = get.vcdf.normal(25, 36),
  obj.prob.vsucc                 = pvacsuccess(c(1:240, seq(252, 1212, 12)),
                                               get.boulianne.vsucc()),
  sia.timing.in.year             = 3 / 12,
  MR1MR2correlation              = TRUE,
  MR1SIAcorrelation              = FALSE,
  MR2SIAcorrelation              = FALSE,
  intro.rate                     = 1 / 24 / 320,
  # --- Observational model ---
  or.reporting.rate              = shared.obs$or.reporting.rate,
  or.non.meas.cases.by.age.month = shared.obs$or.non.meas.cases.by.age.month,
  or.Se                          = shared.obs$or.Se,
  or.Sp                          = shared.obs$or.Sp,
  # --- OBR trigger ---
  or.n.confirmations.target      = 5,
  or.confirmation.delay          = 2,
  or.response.delay              = 6,
  or.trigger.window              = 4,
  or.trigger.window.age.lower    = 0,
  or.trigger.window.age.upper    = 168,
  # --- OBR response (age targets from suspected cases) ---
  or.vacc.coverage               = 0.85,
  or.vacc.agedist.percentile     = 0.95,
  stop_testing_at_trigger        = TRUE,        # use suspected cases for age targeting
  or.min.interval                = 48,
  or.start.timestep              = (2020 - year) * 24 + 1
)


# --- Inspect results ---

# When did OBR campaigns fire?
obr.tp.steps   <- which(result.obr.tp@result@or.times == 1)
obr.susp.steps <- which(result.obr.susp@result@or.times == 1)
cat("OBR campaigns (TP targeting) at years:        ",
    round(year + (obr.tp.steps   - 1) * 0.5 / 12, 2), "\n")
cat("OBR campaigns (suspected targeting) at years: ",
    round(year + (obr.susp.steps - 1) * 0.5 / 12, 2), "\n")

# Compare cumulative incidence
cum.I <- function(res) sum(res@result@.Data[res@result@i.inds, ])
cat("Cumulative I — no OBR:              ", round(cum.I(result.no.obr)),   "\n")
cat("Cumulative I — OBR (TP targeting):  ", round(cum.I(result.obr.tp)),   "\n")
cat("Cumulative I — OBR (susp targeting):", round(cum.I(result.obr.susp)), "\n")
cat("Reduction (TP):   ",
    round((1 - cum.I(result.obr.tp)   / cum.I(result.no.obr)) * 100, 1), "%\n")
cat("Reduction (susp): ",
    round((1 - cum.I(result.obr.susp) / cum.I(result.no.obr)) * 100, 1), "%\n")

# Inspect observational model output for TP scenario (post-2020)
t.start.2020 <- (2020 - year) * 24 + 1
rc.tp        <- result.obr.tp@result
total.TP     <- sum(rc.tp@obs.TP[,     t.start.2020:ncol(rc.tp@obs.TP)])
total.FP     <- sum(rc.tp@obs.FP_test[ t.start.2020:ncol(rc.tp@obs.FP_test)])
cat("Post-2020 IgM-confirmed TP (tested): ", round(total.TP), "\n")
cat("Post-2020 IgM-positive  FP (tested): ", round(total.FP), "\n")
cat("Estimated IgM+ proportion among tested:",
    round(total.TP / (total.TP + total.FP), 3), "\n")
