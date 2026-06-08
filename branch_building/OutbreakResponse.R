
# ============================================================
# Outbreak Response (OBR) — working example for Zambia
# ============================================================
#
# This script demonstrates EX.Country.part2.OR() using Zambia demography
# and vaccination history. It compares two scenarios:
#   A) Routine + SIA vaccination only (no OBR)
#   B) Routine + SIA + OBR triggered when cumulative I exceeds a threshold
#
# Uses R0 = 18 from the MAP calibration to serology data.
# ============================================================

library(MRTransmissionModel)
source("R/setClasses.R")
source("R/run.R")
source("R/EX.Country.part2.OR.R")

setup <- setupCountry.Nov2023(country = "Zambia")
year  <- 1980
t.max <- 45

# --- Step 1: EX.Country.part1 (shared transient run) ---

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


# --- Scenario A: Routine + SIA only (no OBR) ---
# Uses EX.Country.part2.OR with or.vacc.coverage = 0 to disable the response.

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
  time.specific.MR1cov           = setup$MCV1.coverage.1980to2100[(year - 1980 + 1):t.max]*0.5,
  time.specific.MR2cov           = setup$MCV2.coverage.1980to2100[(year - 1980 + 1):t.max]*0.5,
  time.specific.SIAcov           = setup$measlesSIA.coverage.1980to2100[(year - 1980 + 1):t.max],
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
  or.total.delay                 = 4,
  or.trigger.window              = 48,
  or.threshold.type              = "count",
  or.threshold.value             = 500,
  or.trigger.age.lower           = 0,
  or.trigger.age.upper           = 168,
  or.vacc.age.lower              = 6,
  or.vacc.age.upper              = 168,
  or.vacc.coverage               = 0,           # disabled — no OBR response
  or.min.interval                = 48,
  or.start.timestep              = (2020 - year) * 24 + 1
)


# --- Scenario B: Routine + SIA + OBR ---
#
# OBR trigger: if the sum of I in children 0-14 years (0-168 months)
# over the past 2 years (48 half-month steps), as observed with a
# 2-month total lag (4 steps), exceeds 500 infectious person-steps,
# fire a campaign targeting 6 months - 14 years at 85% coverage.
# Campaigns cannot re-trigger within 2 years (48 steps).

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
  time.specific.MR1cov           = setup$MCV1.coverage.1980to2100[(year - 1980 + 1):t.max],
  time.specific.MR2cov           = setup$MCV2.coverage.1980to2100[(year - 1980 + 1):t.max],
  time.specific.SIAcov           = setup$measlesSIA.coverage.1980to2100[(year - 1980 + 1):t.max],
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
  or.total.delay                 = 4,        # 2-month detection + response lag (4 half-month steps)
  or.trigger.window              = 48,       # sum I over past 2 years (48 half-month steps)
  or.threshold.type              = "count",  # absolute infectious person-steps
  or.threshold.value             = 500,      # fire OBR if sum(I) in window >= 500
  or.trigger.age.lower           = 0,        # surveillance: 0-14 years
  or.trigger.age.upper           = 168,      # (168 months = 14 years)
  or.vacc.age.lower              = 6,        # campaign targets 6 months - 14 years
  or.vacc.age.upper              = 168,
  or.vacc.coverage               = 0.85,
  or.min.interval                = 48,       # minimum 2 years between campaigns
  or.start.timestep              = (2020 - year) * 24 + 1  # OBR available from 2020 onward
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

# Compare cumulative incidence (total I summed across all time steps)
cum.I.no.obr <- sum(result.no.obr@result@.Data[result.no.obr@result@i.inds, ])
cum.I.obr    <- sum(result.obr@result@.Data[result.obr@result@i.inds, ])
cat("Cumulative I — no OBR:", round(cum.I.no.obr), "\n")
cat("Cumulative I — with OBR:", round(cum.I.obr), "\n")
cat("Reduction:", round((1 - cum.I.obr / cum.I.no.obr) * 100, 1), "%\n")

