# AI Session Log — MRTransmissionModel

## Branch: model-calibration
## Goal: Build measles serology calibration

---

## Background

This package (`MRTransmissionModel`) is an MSIRV (Maternal/Susceptible/Infected/Recovered/Vaccinated) compartmental measles and rubella transmission model. A branch called `model-calibration` was started to add functionality for calibrating the model to observed serological data.

The calibration fits two parameters by default:
- `R0` — basic reproduction number (optimized on log scale: `theta[1] = log(R0)`)
- `scale.sia` — multiplicative scaling factor on SIA coverage, in (0, 1) (optimized on logit scale: `theta[2] = qlogis(scale.sia)`)

`rho` (RI-SIA dose correlation) can optionally be estimated as a third parameter by setting `fix.rho = NA`. By default it is fixed at 1. See design decisions below.

The entry point for fitting is `FitMeaslesSerology()`. For single evaluations or diagnostic fitting, use `LogLikMeaslesSerology()` directly or `optimize()` for 1-parameter R0-only fits.

---

## Files added/modified on this branch

### New R files
| File | Purpose |
|------|---------|
| `R/GetPredictedMeaslesSerology.R` | Runs Part 1 + Part 2 simulation, returns predicted seroprevalence aligned to observed data rows |
| `R/GetMeaslesSeroprevalence.per.TimePoint.R` | Extracts immune + total population by age group at requested time points only |
| `R/AggregateMeaslesSeroprevalenceByAgeBins.R` | Aggregates age-year seroprevalence into user-specified age bins |
| `R/LogLikMeaslesSerology.R` | `LogLikMeaslesSerology()` — binomial log-likelihood |
| `R/TransformMeaslesSerologyParameters.R` | `TransformMeaslesSerologyParameters()` — maps theta → natural scale (R0, scale.sia, rho) |
| `R/LogLikMeaslesSerology.transformed.R` | `LogLikMeaslesSerology.transformed()` — convenience wrapper on transformed scale |
| `R/NegLogLikMeaslesSerology.transformed.R` | `NegLogLikMeaslesSerology.transformed()` — negative log-lik for use with `optim()` |
| `R/FitMeaslesSerology.R` | `FitMeaslesSerology()` — full MLE fit |

### Modified R files
| File | What changed |
|------|-------------|
| `R/GetNumber.per.AgeYear.R` | Vectorized the inner double for-loop with `colSums(matrix(...))` — was 240 R iterations, now 1 call |
| `R/GetPredictedMeaslesSerology.R` | Added `scale.sia`, `rho`, `fix.rho` parameters; custom VE and vaccination age windows; MCV1/MCV2 scaled by 0.85/0.90 |
| `R/EX.Country.part2.R` | Routing: `SIAinacc=TRUE` goes to `vaccinationlimitations`; `correlation/MR1MR2correlation` goes to `vaccinationcorrelation` |
| `R/run.R` | New `setMethod("run", "experiment.updatedemog.vaccinationchange.vaccinationcorrelation", ...)` added; intro rate indexing fixed |
| `R/setClasses.R` | New S4 class `experiment.updatedemog.vaccinationchange.vaccinationcorrelation` added |
| `DESCRIPTION` | Added new R files to `Collate`; removed phantom `NegLogLikMeaslesSerology.R` entry |

---

## The rho and dose correlation parameters

`rho` is the Pearson correlation between RI (routine immunization) and SIA doses received by a child. The combined vaccination probability at each SIA time step uses:
```
pvacc = r + s - r*s - rho * sqrt(r*(1-r)*s*(1-s))
```
where `r` = routine coverage and `s` = SIA coverage.

- `rho = 0`: independence — combined coverage = r + s - r*s
- `rho = 1`: perfect dependence — children reached by RI are the same ones reached by SIAs; combined coverage ≈ max(r, s)
- `MR1MR2correlation = TRUE`: handled separately in `get.routine.time.age.specific()`; accounts for correlation between MR1 and MR2 routine doses

**In calibration, `rho` defaults to 1** (`MR1SIAcorrelation = 1`, `MR2SIAcorrelation = 1`, `MR1MR2correlation = TRUE`). See the Zambia calibration findings below for why.

### fix.rho — estimating vs. fixing rho

All calibration functions accept a `fix.rho` argument:
- `fix.rho = 1` (default): rho fixed at 1, 2-parameter model (R0, scale.sia)
- `fix.rho = NA`: rho is estimated as a third parameter (R0, scale.sia, rho)

`TransformMeaslesSerologyParameters(theta, fix.rho)` always returns a list with all three fields: `R0`, `scale.sia`, `rho`. When `fix.rho` is a number, `theta` must have length 2 and `rho = fix.rho` is passed through. When `fix.rho = NA`, `theta` must have length 3 and `rho = tanh(theta[3])`.

### S4 class: `experiment.updatedemog.vaccinationchange.vaccinationcorrelation`

Defined in `setClasses.R`. Inherits from `experiment.updatedemog.vaccinationchange`. Slots:
- `MR1MR2correlation` — logical
- `MR1SIAcorrelation` — numeric (the `rho` value for MR1-SIA correlation)
- `MR2SIAcorrelation` — numeric (the `rho` value for MR2-SIA correlation)

`EX.Country.part2()` routes to this class when `MR1MR2correlation=TRUE` or either SIA correlation is non-zero. `SIAinacc=TRUE` routes to `vaccinationlimitations` (not this class).

---

## The scale.sia parameter

`scale.sia` is a multiplicative factor applied to all SIA coverage values before they enter the simulation:
```r
sia.cov <- setup$measlesSIA.coverage.1980to2100[...] * scale.sia
```

- `scale.sia = 1`: SIA coverage at reported levels (default)
- `scale.sia = 0`: SIAs contribute nothing
- `scale.sia = 0.5`: SIAs reach 50% of their reported coverage

**Transform**: `scale.sia = plogis(theta[2])`, constraining to (0, 1).

This parameter was introduced after finding that reported SIA coverage levels (0.89–0.97 per campaign, some covering ages 6 months to 15 years) caused near-perfect predicted vaccination for most birth cohorts when combined with high MCV1 coverage — systematically over-predicting seroprevalence.

---

## Vaccine effectiveness and age windows

Current settings in `GetPredictedMeaslesSerology.R`:

```r
vacc_succ_obj <- pvacsuccess(age.classes, get.boulianne.vsucc())
vacc_succ_obj@prob.vsucc <- rep(1, length(vacc_succ_obj@ages))
```
VE = 1 for all age classes. This is valid because no doses are given before 9 months (MR1 window starts at 10 months), so the VE value for younger ages is never applied.

MCV1 and MCV2 coverage are scaled before entering the simulation:
```r
time.specific.MR1cov = setup$MCV1.coverage.1980to2100[...] * 0.85
time.specific.MR2cov = setup$MCV2.coverage.1980to2100[...] * 0.90
```

Vaccination age windows:
- MR1: 10–12 months (uniform vcdf), `get.vcdf.uniform(10, 12)`
- MR2: 19–25 months (uniform vcdf), `get.vcdf.uniform(19, 25)`

---

## Zambia calibration findings

### Coverage inputs
- **MCV1**: 0% before 1983, then 57–97% from 1983 onward (~85–97% since 1990)
- **MCV2**: 0% before 2014, rising from 33% in 2014 to ~81% by 2021
- **SIAs**: 7 campaigns total (2002, 2003, 2007, 2010, 2012, 2016, 2020); coverage 0.69–0.97 per campaign; age ranges vary — 2002/2003 were very wide (6–180 months = up to 15 years)

### Why rho was fixed at 1 and scale.sia introduced

When fitting with `R0`, `rho`, and `scale.sia` all free, both `rho → 1` and `scale.sia → 0` simultaneously. These parameters are **degenerate** when `scale.sia ≈ 0`: the rho term `rho * sqrt(r*(1-r)*s*(1-s)) ≈ 0` when `s ≈ 0`, so rho has no effect. The optimizer can push both to their boundaries freely.

Diagnostic fits confirmed:
- Fitting with `scale.mcv1` free and SIAs zeroed: `scale.mcv1 → 1` (routine coverage does not need scaling)
- The over-prediction comes specifically from SIA campaigns piling on top of already-high MCV1 for cohorts born 1988–2010

**Conclusion**: Fix `rho = 1` (maximum correlation between all dose types, conservative assumption) and estimate only `R0` and `scale.sia`.

### Current fit (simulated Zambia data, 3 survey years 2015/2019/2023)

With MCV1 × 0.85, MCV2 × 0.90, VE = 1, MR1 window 10–12 months, MR2 window 19–25 months, the optimizer is converging toward R0 ≈ 7.3, scale.sia ≈ 0.18. This is a meaningful non-zero scale.sia (unlike earlier fits), suggesting the MCV coverage scaling breaks the degeneracy enough for SIAs to contribute.

### Concern about 2016 data sensitivity

The 2016 serological survey may use a less sensitive assay than the 2024 survey, causing systematically depressed observed seroprevalence in ages 1–30. This could explain why the optimizer consistently eliminates SIA contributions when fit to the full real Zambia dataset — the model is trying to match artificially low 2016 seroprevalence.

### Recommended calibration strategy (real Zambia data)

Fit to a subset that avoids the potentially low-sensitivity 2016 ages 1–25 cohorts:
- **2016 ages 25+** (born ≤ 1991): immunity is primarily from natural infection; ages 30+ (born ≤ 1986) had no SIA exposure at all — purely R0-driven
- **2024 ages 0–4**: recent cohorts, mainly MCV1/MCV2 driven; no SIA exposure (2020 SIA targeted ages 9–59 months, so 0–4 year olds in 2024 were not reached)

```r
serodata.2016.older <- serodata[
  serodata$survey.time.point == (2016 - 1980) * 24 &
  serodata$age.bin.lower >= 25, ]

serodata.2024 <- data.frame(
  survey.time.point = rep((2024 - 1980) * 24, nrow(seroprev.2024.raw)),
  age.bin.lower     = seroprev.2024.raw$age.integer,
  age.bin.upper     = seroprev.2024.raw$age.integer + 1,
  n_tested          = seroprev.2024.raw$ntest,
  n_positive        = seroprev.2024.raw$npos
)

serodata.subset <- rbind(serodata.2016.older, serodata.2024)
```

Since 2016 ages 30+ have no SIA exposure, `scale.sia` is weakly identified from this subset. Consider fixing `scale.sia = 0` and fitting R0 only using `optimize()`.

### R0-only fit result (real Zambia data)

Fitting to 2016 ages 25+ and 2024 ages 0–4 with `scale.sia = 0` fixed:
```
R0 ≈ 7.9
```
Converged in 11 iterations using Brent's method via `optimize()`.

R0 ≈ 7–8 is consistently lower than typical measles estimates for sub-Saharan Africa (12–18). This may reflect: (1) genuine lower transmission in Zambia's population structure, (2) serological data inconsistency with vaccination history, or (3) assay sensitivity issues in the 2016 survey.

---

## Standard calibration workflow

### 2-parameter fit (R0 + scale.sia, rho fixed at 1)

```r
fit <- FitMeaslesSerology(
  serodata         = serodata,
  setup            = setup,
  year             = 1980,
  t.max            = t.max,
  fix.rho          = 1,
  age.classes      = c(1:60, seq(72, 1212, 12)),
  age0is6to11monly = TRUE,
  hessian          = TRUE
)
fit$par.natural   # $R0, $scale.sia, $rho (= 1)
```

### 3-parameter fit (R0 + scale.sia + rho estimated)

```r
fit <- FitMeaslesSerology(
  serodata         = serodata,
  setup            = setup,
  year             = 1980,
  t.max            = t.max,
  fix.rho          = NA,
  par.init         = c(log(16), qlogis(0.5), atanh(0)),
  age.classes      = c(1:60, seq(72, 1212, 12)),
  age0is6to11monly = TRUE,
  hessian          = TRUE
)
fit$par.natural   # $R0, $scale.sia, $rho
```

### 1-parameter fit (R0 only, scale.sia fixed at 0)

```r
iter <- 0

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
      age0is6to11monly = TRUE
    )
  },
  interval = c(log(2), log(50))
)

list(R0 = exp(fit.R0$minimum), logLik = -fit.R0$objective)
```

### Profile likelihood (fix rho at several values, fit R0 + scale.sia)

To check whether rho is identifiable from the data:

```r
rho.values <- c(0, 0.3, 0.6, 0.9)

profile.results <- lapply(rho.values, function(rho.fixed) {
  fit <- FitMeaslesSerology(
    serodata         = serodata,
    setup            = setup,
    year             = 1980,
    t.max            = t.max,
    fix.rho          = rho.fixed,
    age.classes      = c(1:60, seq(72, 1212, 12)),
    age0is6to11monly = TRUE
  )
  list(rho = rho.fixed, R0 = fit$par.natural$R0,
       scale.sia = fit$par.natural$scale.sia, logLik = fit$logLik)
})
```

---

## Speedup fixes implemented

### Fix 1 — Post-processing only at observed time points
**File**: `GetMeaslesSeroprevalence.per.TimePoint.R`

Only processes the unique time points in `serodata$survey.time.point` (typically 2–5), not all ~1,081 simulation time steps. Roughly 360× faster for that function.

### Fix 2 — Vectorized age aggregation
**File**: `GetNumber.per.AgeYear.R`

Replaced a double `for` loop (240 R iterations) with:
```r
pop.per.youngage.year <- colSums(matrix(vec[seq_len(top.one.month.age)], nrow = 12))
```

---

## Age class coarsening

Default: `c(1:240, seq(252, 1212, 12))` = 321 classes → state vector 1,605 rows
Coarse:  `c(1:60, seq(72, 1212, 12))` = 156 classes → state vector 780 rows (~4× speedup)

Pass `age.classes` directly to `FitMeaslesSerology` and it propagates all the way through to `EX.Country.part1` / `EX.Country.part2`.

---

## Package loading workflow

**`devtools::load_all(".")` crashes R on this machine — do not use it.**

```r
devtools::document()           # regenerates NAMESPACE + man/
devtools::build()              # creates .tar.gz (does NOT install)
devtools::install(".", quick=TRUE)  # bakes new R files into library()
library(MRTransmissionModel)

# Source calibration files not yet in installed package:
source("R/GetPredictedMeaslesSerology.R")
source("R/GetMeaslesSeroprevalence.per.TimePoint.R")
source("R/AggregateMeaslesSeroprevalenceByAgeBins.R")
source("R/LogLikMeaslesSerology.R")
source("R/TransformMeaslesSerologyParameters.R")
source("R/LogLikMeaslesSerology.transformed.R")
source("R/NegLogLikMeaslesSerology.transformed.R")
source("R/FitMeaslesSerology.R")
source("R/EX.Country.part2.R")
source("R/run.R")
```

---

## Posterior uncertainty — Normal (Laplace) approximation

Pass `hessian = TRUE` to `FitMeaslesSerology`. The return object includes `fit$cov.transformed` — covariance matrix on the transformed scale (`log_R0`, `logit_scale.sia`).

```r
samples.transformed <- MASS::mvrnorm(
  n = 2000, mu = fit$par.transformed, Sigma = fit$cov.transformed
)
posterior.samples <- data.frame(
  R0        = exp(samples.transformed[, "log_R0"]),
  scale.sia = plogis(samples.transformed[, "logit_scale.sia"])
)
```

**Validity check**: if marginal histograms look skewed or any parameter is near its boundary, the normal approximation is unreliable — upgrade to MCMC.

---

## Future work: full Bayesian MCMC ⚠️ REMINDER

**Not yet done — worth doing when time allows.**

**Why deferred**: ~25–30 sec per likelihood evaluation → ~35 hours for 5,000 iterations. Feasible overnight but not interactive.

**How to implement**:
- MH on `(log(R0), logit(scale.sia))`, initialized at the MLE
- Proposal covariance: `2.38²/2 * fit$cov.transformed`
- `adaptMCMC` R package handles adaptive tuning automatically
- Run 3+ chains for convergence diagnostics (R-hat)

---

## Approaches tried and abandoned

### rho as free parameter
Initially `rho` (Pearson RI-SIA correlation) was estimated alongside R0. When fitting to full Zambia data (2016 + 2024), `rho → 1` consistently. When `scale.sia` was added as a third parameter, both `rho → 1` and `scale.sia → 0` simultaneously — a degenerate solution where rho is irrelevant. Rho is now fixed at 1 by default (`fix.rho = 1`), but can be freed with `fix.rho = NA`.

### prop.inacc as free parameter
A `prop.inacc` parameter (fraction of population permanently unreachable by all vaccination) was implemented and tested. Fitting gave `R0 = 6.16, rho ≈ 1, prop.inacc = 0.40` — all three parameters at implausible extremes. Diagnosis: coverage inputs and rho/prop.inacc are not jointly identifiable. Replaced by `scale.sia` (scales SIA coverage only), which is more targeted.

### scale.mcv1 as free parameter
Tested fitting with `scale.mcv1` (scaling routine MCV1+MCV2) and SIAs zeroed. The optimizer converged to `scale.mcv1 ≈ 1` — routine coverage does not need scaling. The problem is specifically SIA coverage being too high relative to what serology implies.

### Part 1 result cache
A `.cache` environment was added to `GetPredictedMeaslesSerology()` to skip the 20-year Part 1 burn-in when R0 was unchanged between calls. Removed because the cache required exact floating-point equality of R0 to hit, and Nelder-Mead changes R0 at every step — so the cache never fired during optimization. It would only have helped with an explicit grid search (also removed). The two speedups that do help are Fix 1 and Fix 2 above.

---

## Key design decisions

- `survey.time.point` in `serodata` must be a **time-step index** (not a calendar year). For `year=1980`, `generation.time=0.5` (24 steps/year): calendar year `y` → index `(y - 1980) * 24`.
- `age.bin.lower` and `age.bin.upper` are in **single-year units** (0, 1, 2, ...), not months.
- `rho` defaults to 1 in calibration (`MR1SIAcorrelation = 1`, `MR2SIAcorrelation = 1`, `MR1MR2correlation = TRUE`). All dose types are assumed maximally correlated — conservative assumption that minimizes the joint vaccination probability.
- `scale.sia` is applied to the SIA coverage vector before it enters `EX.Country.part2()`: `sia.cov * scale.sia`. It does not affect MCV1 or MCV2.
- `rho` is available as a free parameter via `fix.rho = NA` in all calibration functions, but is degenerate with `scale.sia` when fit to the full Zambia dataset.
- For 1-parameter fits (R0 only, scale.sia fixed), use `optimize()` with Brent's method — faster and cleaner than `optim()` for 1D problems.
- `age0is6to11monly = TRUE` should be passed for the Zambia real data (2024 survey ages 0–4 includes infants; age-0 seroprevalence is computed from age classes 7–12, i.e. months 6–11, only).
- No `...` in any function signature — all arguments are listed explicitly in every function in the call chain.
