
# ============================================================
# Outbreak Response (OBR) — working example for Zambia
# ============================================================
#
# Compares two scenarios using Zambia demography and vaccination history:
#   A) Routine + SIA only (no OBR):  or.vacc.coverage = 0
#   B) Routine + SIA + OBR:          or.vacc.coverage = 0.85
#
# Both scenarios use identical routine and SIA vaccination so that any
# difference in outcomes is attributable solely to the OBR campaign.
#
# R0 = 18 from MAP calibration to Zambia serology data.
# Step size: 0.5 months (24 steps/year).
# ============================================================

library(MRTransmissionModel)
devtools::document()           # regenerate NAMESPACE (run once after @export changes)
devtools::build()              # creates .tar.gz (does NOT install)
#source("R/GetPredictedMeaslesSerology.R")

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
# Used identically in both scenarios so OBR impact is cleanly isolated.

mr1cov <- setup$MCV1.coverage.1980to2100[(year - 1980 + 1):t.max]*0.5
mr2cov <- setup$MCV2.coverage.1980to2100[(year - 1980 + 1):t.max]*0.5
siacov <- setup$measlesSIA.coverage.1980to2100[(year - 1980 + 1):t.max]*0.5


# --- Scenario A: Routine + SIA only (no OBR) ---
# or.vacc.coverage = 0 disables the response campaign while keeping all
# other OBR arguments syntactically present for easy comparison.

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
  # --- OBR parameters ---
  or.trigger.mode                = "I_scaled",   # use true I x reporting rate
  or.reporting.rate              = 1,            # no scaling (raw I)
  or.n.confirmations.target      = 5,
  or.confirmation.delay          = 2,            # 1-month lag to case confirmation
  or.response.delay              = 2,            # 1-month lag from trigger to campaign
  or.trigger.window              = 48,           # 2-year lookback window
  or.trigger.age.lower           = 0,            # surveillance: all ages 0-14 years
  or.trigger.age.upper           = 168,
  or.vacc.age.lower              = 6,            # campaign: 6 months - 14 years
  or.vacc.age.upper              = 168,          # fixed upper bound (overrides CDF)
  or.vacc.agedist.percentile     = NA_real_,     # not used when or.vacc.age.upper is set
  or.vacc.coverage               = 0,            # DISABLED — no OBR response
  or.min.interval                = 48,           # minimum 2 years between campaigns
  or.start.timestep              = (2020 - year) * 24 + 1  # OBR available from 2020
)


# --- Scenario B: Routine + SIA + OBR ---
#
# OBR trigger: sum of I x reporting.rate in children 0-14 years over the
# past 2 months (4 half-month steps), observed with a 1-month confirmation
# lag, reaches or.n.confirmations.target = 5. Campaign fires 1 month later
# (or.response.delay = 2 steps), targeting 6 months - 14 years at 85%
# coverage. Campaigns cannot re-trigger within 2 years (48 steps).

result.obr <- EX.Country.part2.OR(
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
  # --- OBR parameters ---
  or.trigger.mode                = "I_scaled",
  or.reporting.rate              = 1,
  or.n.confirmations.target      = 5,
  or.confirmation.delay          = 2,
  or.response.delay              = 2,
  or.trigger.window              = 4,
  or.trigger.age.lower           = 0,
  or.trigger.age.upper           = 168,
  or.vacc.age.lower              = 6,
  or.vacc.age.upper              = 168,
  or.vacc.agedist.percentile     = NA_real_,
  or.vacc.coverage               = 0.85,
  or.min.interval                = 48,
  or.start.timestep              = (2020 - year) * 24 + 1
)


# --- Scenario C: Routine + SIA + OBR (confirmed trigger mode) ---
#
# Uses the full observational model pipeline to trigger OBR:
#   1. True infections are scaled by or.reporting.rate to give suspected cases
#      that reach the surveillance system (care-seeking x clinical recognition,
#      before any testing).
#   2. Non-measles rash-febrile illness (or.non.meas.cases.by.age.month)
#      is added to get total suspected cases.
#   3. Tests are allocated across age groups proportionally. Testing stops once
#      or.n.confirmations.target / pos.rate tests have been done.
#   4. Confirmed cases (obs.TP + obs.FP_test + obs.TP_clinical + obs.FP_clinical)
#      are summed; trigger fires when the or.trigger.window sum >= or.n.confirmations.target.
#
# Se = 0.95, Sp = 0.97 are standard measles IgM values from the literature.
# Background rash: exponentially declining with age (peak in infants),
# ~20% higher in Zambia's rainy season (Nov-Apr = months 11, 12, 1-4).
# or.reporting.rate = 0.10: ~10% of true infections present as suspected cases
# to the health system (care-seeking x clinical recognition, before testing).

age.classes.ex <- c(1:240, seq(252, 1212, 12))
n.age.ex       <- length(age.classes.ex)

age.bg.rate  <- 0.4 * exp(-age.classes.ex / 48)   # cases per step per age class; peak in infants
month.mult   <- c(1.15, 1.15, 1.10, 1.05, 1.00, 0.90,  # Jan-Jun
                  0.85, 0.85, 0.90, 0.95, 1.05, 1.15)  # Jul-Dec
non.meas.bg  <- outer(age.bg.rate, month.mult)     # n.age x 12 matrix

result.obr.confirmed <- EX.Country.part2.OR(
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
  # --- OBR parameters (confirmed trigger mode) ---
  or.trigger.mode                = "confirmed",
  or.reporting.rate              = 0.10,           # care-seeking x clinical recognition (pre-testing)
  or.non.meas.cases.by.age.month = non.meas.bg,   # non-measles background (n.age x 12 matrix)
  or.Se                          = 0.95,           # IgM test sensitivity
  or.Sp                          = 0.97,           # IgM test specificity
  or.n.confirmations.target      = 5,
  or.confirmation.delay          = 2,
  or.response.delay              = 2,
  or.trigger.window              = 4,
  or.trigger.age.lower           = 0,
  or.trigger.age.upper           = 168,
  or.vacc.age.lower              = 6,
  or.vacc.age.upper              = 168,
  or.vacc.agedist.percentile     = NA_real_,
  or.vacc.coverage               = 0.85,
  or.min.interval                = 48,
  or.start.timestep              = (2020 - year) * 24 + 1
)


# --- Inspect results ---

# When did OBR campaigns fire?
obr.timesteps <- which(result.obr@result@or.times == 1)
obr.years     <- year + (obr.timesteps - 1) * 0.5 / 12
cat("OBR campaigns fired at years:", round(obr.years, 2), "\n")

# Pre-planned SIA years (same in both scenarios)
sia.timesteps <- which(result.obr@result@sia.times == 1)
sia.years     <- year + (sia.timesteps - 1) * 0.5 / 12
cat("SIA campaigns at years:", round(sia.years, 2), "\n")

# Compare cumulative incidence across all three scenarios
cum.I.no.obr   <- sum(result.no.obr@result@.Data[result.no.obr@result@i.inds, ])
cum.I.obr      <- sum(result.obr@result@.Data[result.obr@result@i.inds, ])
cum.I.confirmed <- sum(result.obr.confirmed@result@.Data[result.obr.confirmed@result@i.inds, ])
cat("Cumulative I — no OBR:           ", round(cum.I.no.obr), "\n")
cat("Cumulative I — OBR (I_scaled):   ", round(cum.I.obr), "\n")
cat("Cumulative I — OBR (confirmed):  ", round(cum.I.confirmed), "\n")
cat("Reduction (I_scaled):  ", round((1 - cum.I.obr       / cum.I.no.obr) * 100, 1), "%\n")
cat("Reduction (confirmed): ", round((1 - cum.I.confirmed / cum.I.no.obr) * 100, 1), "%\n")

# When did confirmed-mode OBR campaigns fire?
obr.conf.timesteps <- which(result.obr.confirmed@result@or.times == 1)
obr.conf.years     <- year + (obr.conf.timesteps - 1) * 0.5 / 12
cat("Confirmed-mode OBR campaigns fired at years:", round(obr.conf.years, 2), "\n")

# Inspect observational model output for the confirmed scenario.
# Total confirmed cases (true + false positive) summed over all ages and the
# post-2020 period. Useful for checking whether the obs model is firing at
# a plausible rate.
t.start.2020 <- (2020 - year) * 24 + 1
rc.conf      <- result.obr.confirmed@result
total.TP     <- sum(rc.conf@obs.TP[,      t.start.2020:ncol(rc.conf@obs.TP)])
total.FP     <- sum(rc.conf@obs.FP_test[, t.start.2020:ncol(rc.conf@obs.FP_test)])
total.TP.cli <- sum(rc.conf@obs.TP_clinical[, t.start.2020:ncol(rc.conf@obs.TP_clinical)])
cat("Post-2020 confirmed true positives (tested):     ", round(total.TP), "\n")
cat("Post-2020 confirmed false positives (tested):    ", round(total.FP), "\n")
cat("Post-2020 true positives (clinical/epi-linked):  ", round(total.TP.cli), "\n")
cat("Estimated IgM+ proportion among tested cases:",
    round(total.TP / (total.TP + total.FP), 3), "\n")
