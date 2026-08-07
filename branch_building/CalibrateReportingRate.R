
# ============================================================
# Calibrate or.reporting.rate (r) for Zambia — 5-rule + Rt > 1 anchor
# ============================================================
#
# Goal: estimate or.reporting.rate so that when the surveillance system
# accumulates 5 IgM-confirmed cases (obs.TP) over the trigger window, the
# underlying true infections are at an epidemiologically meaningful level —
# i.e., Rt has just crossed above 1 (epidemic onset).
#
# Method (from Building_Observational_Model.md):
#
#   SIMPLE.FORMULA = TRUE  (r absorbs everything):
#     r = 5 / n_infections_at_Rt_crossing
#
#   SIMPLE.FORMULA = FALSE (model-consistent; Se and window applied separately):
#     r = 5 / (n_infections_at_Rt_crossing × Se × or.trigger.window)
#     because the obs model computes obs.TP = I × r × Se and the trigger
#     sums obs.TP over or.trigger.window steps.
#
# Procedure:
#   1. Spin up with EX.Country.part1() (deterministic, shared).
#   2. Run n.reps stochastic replicates with EX.Country.part2() — no OBR
#      machinery needed, just raw transmission + vaccination dynamics.
#   3. For each run, find epidemic onset moments in the post-eval.year window:
#      first time I grows for n.consec consecutive steps after a trough.
#   4. Record total true I at each onset, pool across all runs.
#   5. Apply the 5-rule.
# ============================================================

library(MRTransmissionModel)

# ============================================================
# User settings
# ============================================================

n.reps         <- 200      # number of stochastic replicates
year           <- 1980     # simulation start year (must match EXt0)
t.max          <- 45       # years to simulate (1980–2025)
hist.year      <- 2015     # use actual coverage up to (but not including) this year
eval.year      <- 2015     # earliest year to look for epidemic onsets (10-year window: 2015–2025)
n.consec       <- 3        # consecutive growing steps required to call Rt > 1 (3 steps = 6 weeks)
or.Se          <- 0.95     # IgM sensitivity (only used when SIMPLE.FORMULA = FALSE)
or.trigger.win <- 4        # trigger window steps (only used when SIMPLE.FORMULA = FALSE)
n.target       <- 5        # confirmations target (or.n.confirmations.target)
SIMPLE.FORMULA <- TRUE     # TRUE:  r = n.target / mean(I_at_onset)
                           # FALSE: r = n.target / (mean(I_at_onset) × Se × trigger.window)

set.seed(42)

# ============================================================
# Step 1: Setup and deterministic spin-up
# ============================================================

setup <- setupCountry.Nov2023(country = "Zambia")

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

# Coverage strategy: actual historical estimates 1980–2014 to build the
# realistic 2015 susceptibility profile, then zero all vaccination 2015–2025
# so that outbreaks occur in the eval window for onset detection.
hist.idx <- seq_len(hist.year - year)          # indices for 1980–2014 (years 1–35)
eval.idx <- (hist.year - year + 1):t.max       # indices for 2015–2024 (years 36–45)

mr1cov <- setup$MCV1.coverage.1980to2100[(year - 1980 + 1):t.max]
mr2cov <- setup$MCV2.coverage.1980to2100[(year - 1980 + 1):t.max]
siacov <- setup$measlesSIA.coverage.1980to2100[(year - 1980 + 1):t.max]

mr1cov[eval.idx] <- 0
mr2cov[eval.idx] <- 0
siacov[eval.idx] <- 0

t.eval.start <- (eval.year - year) * 24 + 1

# ============================================================
# Step 2: Helper — find epidemic onset timesteps
# ============================================================
#
# Returns indices (relative to the input I.total vector) of epidemic onsets.
# An onset = first step of a run of >= n.consec consecutive steps where
# I(t) > I(t-1), not already inside a flagged epidemic.
# Epidemic state resets once I drops below reset.fraction * I_at_onset.

find_epidemic_onsets <- function(I.total, n.consec = 3, reset.fraction = 0.5) {
  n       <- length(I.total)
  growing <- c(FALSE, diff(I.total) > 0)

  onsets  <- integer(0)
  in.epi  <- FALSE
  onset.I <- NA_real_

  for (t in n.consec:n) {
    sustained <- all(growing[(t - n.consec + 1):t])

    if (sustained && !in.epi) {
      onset.t <- t - n.consec + 1
      onsets  <- c(onsets, onset.t)
      in.epi  <- TRUE
      onset.I <- I.total[onset.t]
    }

    if (in.epi && !is.na(onset.I)) {
      if (I.total[t] < reset.fraction * onset.I || I.total[t] < 1) {
        in.epi  <- FALSE
        onset.I <- NA_real_
      }
    }
  }

  return(onsets)
}

# ============================================================
# Step 3: Run stochastic replicates with EX.Country.part2()
# ============================================================

all.onset.I <- numeric(0)

cat("Running", n.reps, "stochastic replicates...\n")

for (rep in seq_len(n.reps)) {

  if (rep %% 20 == 0) cat("  replicate", rep, "/", n.reps, "\n")

  result <- EX.Country.part2(
    uncode                         = setup$uncode,
    generation.time                = 0.5,
    pop.rescale                    = setup$pop.total.1950.2100[
                                       seq(1990, (year + t.max - 10), 10) - 1950 + 1],
    pop.time                       = seq(1990, (year + t.max - 10), 10) - year,
    is.stochastic                  = TRUE,
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
    intro.rate                     = 1 / 24 / 320
  )

  I.mat   <- result@result@.Data[result@result@i.inds, ]
  I.total <- colSums(I.mat)
  I.eval  <- I.total[t.eval.start:length(I.total)]

  onsets <- find_epidemic_onsets(I.eval, n.consec = n.consec)
  if (length(onsets) > 0) {
    all.onset.I <- c(all.onset.I, I.eval[onsets])
  }
}

# ============================================================
# Step 4: Apply the 5-rule
# ============================================================

cat("\n--- Results ---\n")
cat("Total epidemic onsets detected:", length(all.onset.I), "\n")
cat("Mean I at onset:   ", round(mean(all.onset.I), 1), "\n")
cat("Median I at onset: ", round(median(all.onset.I), 1), "\n")
cat("SD I at onset:     ", round(sd(all.onset.I), 1), "\n")

if (SIMPLE.FORMULA) {
  r.hat <- n.target / mean(all.onset.I)
  cat("\nFormula: r = n.target / mean(I_at_onset)\n")
} else {
  r.hat <- n.target / (mean(all.onset.I) * or.Se * or.trigger.win)
  cat("\nFormula: r = n.target / (mean(I_at_onset) × Se × trigger.window)\n")
}

cat("Calibrated or.reporting.rate:", round(r.hat, 4), "\n")
cat("As a percentage:", round(r.hat * 100, 2), "%\n")

# ============================================================
# Step 5: Diagnostics
# ============================================================

par(mfrow = c(1, 2))

hist(all.onset.I,
     breaks = 30, col = "steelblue", border = "white",
     main = "True I at epidemic onset",
     xlab = "Total I at Rt > 1 crossing",
     ylab = "Count across stochastic runs")
abline(v = mean(all.onset.I),   col = "red",    lwd = 2, lty = 1)
abline(v = median(all.onset.I), col = "orange", lwd = 2, lty = 2)
legend("topright",
       legend = c(paste("mean =",   round(mean(all.onset.I),   1)),
                  paste("median =", round(median(all.onset.I), 1))),
       col = c("red", "orange"), lty = c(1, 2), lwd = 2, bty = "n")

if (SIMPLE.FORMULA) {
  implied.confirmed <- all.onset.I * r.hat
  xlab.str <- "I × r (target = 5)"
} else {
  implied.confirmed <- all.onset.I * r.hat * or.Se * or.trigger.win
  xlab.str <- "I × r × Se × window (target = 5)"
}
hist(implied.confirmed,
     breaks = 30, col = "darkorange", border = "white",
     main = "Implied confirmed cases at onset\n(should be ~5)",
     xlab = xlab.str, ylab = "Count")
abline(v = n.target, col = "red", lwd = 2, lty = 2)

par(mfrow = c(1, 1))
