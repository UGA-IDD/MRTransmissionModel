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
| `R/TransformMeaslesSerologyParameters.R` | `TransformMeaslesSerologyParameters()` — maps theta → natural scale (R0, scale.sia, rho, and optionally phi, kappa) |
| `R/LogLikMeaslesSerology.transformed.R` | `LogLikMeaslesSerology.transformed()` — convenience wrapper on transformed scale |
| `R/NegLogLikMeaslesSerology.transformed.R` | `NegLogLikMeaslesSerology.transformed()` — negative log-lik for use with `optim()` |
| `R/FitMeaslesSerology.R` | `FitMeaslesSerology()` — full MLE fit (serology only) |
| `R/GetAgeSpecificAnnualIncidence.R` | Extracts ΔR from result matrix to compute age-specific annual incidence |
| `R/GetPredictedMeaslesSerologyAndCases.R` | Runs simulation once, returns both serology and case incidence predictions |
| `R/LogLikMeaslesCases.R` | `LogLikMeaslesCases()` — NegBin log-likelihood for age-stratified case counts |
| `R/LogLikMeaslesSerologyAndCases.R` | `LogLikMeaslesSerologyAndCases()` — combined serology + case log-likelihood |
| `R/NegLogLikMeaslesSerologyAndCases.transformed.R` | Negative combined log-lik on transformed scale for use with `optim()` |
| `R/FitMeaslesSerologyAndCases.R` | `FitMeaslesSerologyAndCases()` — full joint MLE fit (serology + cases) |

### Modified R files
| File | What changed |
|------|-------------|
| `R/GetNumber.per.AgeYear.R` | Vectorized the inner double for-loop with `colSums(matrix(...))` — was 240 R iterations, now 1 call |
| `R/GetPredictedMeaslesSerology.R` | Added `scale.sia`, `rho`, `fix.rho` parameters; custom VE and vaccination age windows; MCV1/MCV2 scaled by 0.85/0.90 |
| `R/EX.Country.part2.R` | Routing: `SIAinacc=TRUE` goes to `vaccinationlimitations`; `correlation/MR1MR2correlation` goes to `vaccinationcorrelation`; `as.numeric()` wrapping on `MR1SIAcorrelation`/`MR2SIAcorrelation` assignment for backward compatibility |
| `R/run.R` | New `setMethod("run", "experiment.updatedemog.vaccinationchange.vaccinationcorrelation", ...)` added; intro rate indexing fixed |
| `R/setClasses.R` | New S4 class `experiment.updatedemog.vaccinationchange.vaccinationcorrelation` added |
| `R/FitMeaslesSerology.R` | Added `prior.R0.meanlog = log(14)` and `prior.R0.sdlog = 0.4` arguments for MAP consistency with joint fit |
| `R/NegLogLikMeaslesSerology.transformed.R` | Added R0 prior (`dlnorm`) applied when `fix.R0 = NA`; fixed misleading docstring |
| `DESCRIPTION` | Added new R files to `Collate`; removed phantom `NegLogLikMeaslesSerology.R` entry |
| `NAMESPACE` | Added 13 new calibration function exports, `exportClasses(experiment.updatedemog.vaccinationchange.vaccinationcorrelation)`, and `importFrom(stats, dlnorm)` |
| `.gitignore` | Added `src/*.o`, `src/*.so`, `src/*.dll` to prevent compiled binary files from being tracked |

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

### Candidate calibration strategy: subset by assay concern (exploratory only)

This approach was explored but **not adopted** as the current calibration. It is documented for reference.

The 2016 assay may have lower sensitivity than the 2024 assay, depressing observed seroprevalence in ages 1–25. One response is to fit only to age ranges less sensitive to this:
- **2016 ages 25+** (born ≤ 1991): primarily natural infection; ages 30+ (born ≤ 1986) had no SIA exposure
- **2024 ages 0–4**: recent cohorts, mainly MCV1/MCV2-driven

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

**Current approach**: all age groups from both the 2016 and 2024 surveys are used. The potential assay sensitivity difference between surveys is treated as a model limitation rather than a reason to subset.

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
source("R/GetAgeSpecificAnnualIncidence.R")
source("R/GetPredictedMeaslesSerologyAndCases.R")
source("R/LogLikMeaslesCases.R")
source("R/LogLikMeaslesSerologyAndCases.R")
source("R/NegLogLikMeaslesSerologyAndCases.transformed.R")
source("R/FitMeaslesSerologyAndCases.R")
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

## Future work: incorporating case data

Available case data:
- Suspected cases by age per year (2005–2025)
- Suspected cases by year and month, aggregated over age (2005–2025)
- Notable outbreaks: 2010/11 and 2024/25

### Structural limitation: deterministic model cannot reproduce outbreaks

The model runs with `is.stochastic = FALSE`. A deterministic model predicts the *expected* endemic trajectory — a smooth curve. Outbreaks in 2010/11 and 2024/25 are stochastic events requiring susceptible accumulation plus a spark and Rt crossing 1 transiently. No choice of R0 or scale.sia will make a deterministic model reproduce sharp outbreak peaks; an optimizer fitting to outbreak case counts would distort parameters trying to match averages.

To genuinely capture outbreak dynamics: `is.stochastic = TRUE` + many stochastic replicates + a particle filter (SMC/sequential Monte Carlo) for the likelihood. This is a major undertaking.

### What case data can add within the current deterministic framework

**Reporting fraction from annual case counts**

Add a reporting fraction φ (logit-transformed) and a negative-binomial log-likelihood on annual cases:
```
cases_obs[y] ~ NegBin(phi * incidence_pred[y], overdispersion)
```
Joint log-likelihood = serology LL + case LL. This lets case data inform average transmission intensity independently of serology. Negative-binomial (not Poisson) is necessary to absorb outbreak-year variance that the deterministic model structurally cannot explain. Exclude 2010/11 and 2024/25 from this likelihood, or treat them as outliers.

**Seasonal amplitude from monthly case data**

Monthly aggregate cases (no age breakdown needed) are well-suited for calibrating `seasonal.amp`, currently fixed at 0.15. Use non-outbreak years only (e.g., 2005–2009, 2012–2023). Fit the within-year shape while fixing year-level totals.

### Recommended approach

1. Keep serology calibration as-is for R0 and scale.sia
2. Add φ (reporting fraction) and NegBin overdispersion estimated from annual case totals in non-outbreak years
3. Use monthly cases in non-outbreak years to calibrate `seasonal.amp`
4. Treat 2010/11 and 2024/25 outbreak years as unexplained variance — do not try to fit them with the deterministic model

Note on age distribution: if distributing monthly cases by age is needed, assuming the year-average age structure is applied to all months in that year is workable but adds noise. Start with aggregated annual totals.

---

## R0 identification problem and strategies ⚠️ CURRENT PRIORITY

R0 is weakly identified by the available serology data. The core issue:

- **Young ages (<25)**: seroprevalence is high but vaccination explains it — R0 contributes almost nothing to the likelihood
- **Old ages (>25, 2016 survey)**: seroprevalence reflects natural infection, but at ~90% the curve is flat — many R0 values (8, 12, 16, 20) all produce ~90% cumulative immunity by age 30 in a pre-vaccination setting

The optimizer lands at R0 ≈ 7–8 not because that's the true value but because the serology surface is nearly flat with respect to R0. The three strategies below address this.

---

### Strategy 1 — Age-stratified case counts ✅ IMPLEMENTED

**Why it helps where serology doesn't**: Cases concentrate in the *currently susceptible* age groups. The age distribution of cases is a direct function of R0 and vaccination coverage:
- Low R0 (7–8): susceptibles accumulate slowly → cases spread across a wider age range, more in older children and adults
- High R0 (15–18): high force of infection → cases concentrate tightly in young unvaccinated children

This signal is independent of the serology flatness problem. Age-stratified case data (suspected cases by age per year, 2005–2025) are already available.

**Implementation**: The model predicts age-specific incidence. Annual incidence by age group can be extracted from the R compartment:

```r
# incidence in age group a, year y = change in R compartment over that year
incidence.pred[a, y] <- sum(R.compartment[age.rows.for.a, end.of.year.y]) -
                         sum(R.compartment[age.rows.for.a, start.of.year.y])
```

Add a NegBin log-likelihood on age-specific case counts:

```r
ll.cases.age <- sum(dnbinom(
  x    = cases.obs[a, y],       # observed suspected cases by age and year
  mu   = phi * incidence.pred[a, y],
  size = kappa,
  log  = TRUE
))

ll.total <- ll.serology + ll.cases.age
```

New parameters: `phi` (reporting fraction, logit scale) and `kappa` (NegBin overdispersion, log scale).

**Year exclusions** (from inspecting annual totals):
- **2006**: total = 0, almost certainly missing data — exclude
- **2010, 2011**: 10,779 and 12,726 cases (20–40× background) — major outbreak, exclude
- **2022–2025**: 1,180–1,520 at peak (3–5× background), likely COVID-19 vaccination disruption catchup — exclude or include and let κ absorb
- **2026**: partial year and beyond model range (t.max = 45 → simulation ends 2025) — exclude

Clean fitting years: **2003–2005, 2007–2009, 2012–2021** (14 years of stable endemic transmission). Also exclude **2009** — see pop.rescale artifact below.

**Fitting window comparison (Zambia real data)**:
- Including 2022–2024: converged to R0 ≈ 15.5, kappa ≈ 0.485 (high overdispersion — absorbing COVID-era disruption years)
- Excluding 2022–2024 (`fit2`, years 2003–2005, 2007–2009, 2012–2021): converged to R0 = 16.83, scale.sia = 0.60, **phi = 1 (boundary), kappa ≈ 0** — degenerate solution; see below

**pop.rescale artifact — 2009 spike** ⚠️

`GetPredictedMeaslesSerologyAndCases` originally passed `pop.rescale` to `EX.Country.part2`. The rescaling events at simulation years 10/20/30 (calendar 1990/2000/2010) correspond to time steps 240/480/720. Calendar year 2009 has `t.end = (2009-1980)*24 + 23 = 720` — exactly the rescaling step. The R compartment at t=720 is post-rescale while t.start is pre-rescale, so ΔR for 2009 absorbs the entire population jump → predicted incidence = 205,734 (vs ~58,000 for surrounding years).

**Fix**: `pop.rescale = NULL` hardcoded in `GetPredictedMeaslesSerologyAndCases`. Skipping rescaling is acceptable for calibration since we are fitting epidemic dynamics, not absolute population size. Exclude 2009 from `casedata` regardless as a precaution.

**phi/kappa degeneracy — `fit2` result** ⚠️

Result: R0 = 16.83, scale.sia = 0.60, phi = 1 (boundary), kappa ≈ 0. This is an optimization failure with two signals:

- **phi → 1 (boundary)**: the logit transform constrains phi to (0, 1). Hitting 1 means the optimizer wants phi > 1, which is impossible. This indicates the model predicts far more incidence than observed cases even at 100% reporting — the optimizer is trying to scale UP predictions but can't.
- **kappa → 0 (NegBin collapse)**: with kappa → 0, variance → ∞, making the NegBin so diffuse it barely penalizes any prediction. The optimizer escapes a bad case fit by rendering the case likelihood nearly uninformative.

**Root cause**: The model over-predicts incidence by ~100×. Observed/predicted ratios by year:

```
year   obs    pred     ratio
2003   399   68,463   0.006
2004   675   65,834   0.010
2007   878   58,866   0.015
...
2021   166   12,215   0.014
```

True phi ≈ 0.006–0.015 (0.6–1.5% reporting), not 0.1 (the default starting value). Starting from phi = 0.1, the optimizer moved in the wrong direction. This is a bad starting value problem, not a model problem.

**Fixes applied**:
1. `pop.rescale = NULL` in `GetPredictedMeaslesSerologyAndCases` (removes 2009 artifact)
2. kappa floor of 0.1 added in `TransformMeaslesSerologyParameters` (`pars$kappa <- max(exp(...), 0.1)`) — prevents NegBin collapse
3. Compute data-informed phi starting value before fitting:
```r
phi.init <- sum(casedata$cases) / sum(pred.from.trial$casedata$pred.incidence)
# use sum(), NOT median() — median gives Inf when pred.incidence = 0 for some rows
```
4. Use informed `par.init`:
```r
par.init = c(log(12), qlogis(0.5), qlogis(phi.init), log(1))
```

**Case data format**:
```r
# A tibble: rows × 3
#   age   year  cases
#   <dbl> <dbl> <dbl>
#       0  2003    36
#       0  2004    32
#       1  2003    ...
```
`age` is integer years (0, 1, 2, ...), `year` is calendar year, `cases` is suspected case count.

**Mapping case data to model output**:
- Case `age = a` (years) → model R-compartment rows for age classes `(a*12 + 1)` to `((a+1)*12)` (months)
- Case `year = y` → model time steps `(y - 1980)*24 + 1` to `(y - 1980)*24 + 24`
- Annual incidence for age `a` in year `y` ≈ ΔR for those rows over those columns:

```r
GetAgeSpecificAnnualIncidence <- function(res, r.inds, year.sim, years.case, ages.case) {
  # res      — full result matrix (rows = states, cols = time steps)
  # r.inds   — row indices of R compartment (from tmp.res@result@r.inds)
  # year.sim — simulation start year (1980)
  # years.case, ages.case — vectors of years and ages from casedata

  r.mat <- res[r.inds, , drop = FALSE]   # R compartment only

  n <- length(years.case)
  incidence <- numeric(n)

  for (i in seq_len(n)) {
    age.rows  <- (ages.case[i] * 12 + 1):((ages.case[i] + 1) * 12)
    t.start   <- (years.case[i] - year.sim) * 24 + 1
    t.end     <- t.start + 23

    r.start <- sum(r.mat[age.rows, t.start])
    r.end   <- sum(r.mat[age.rows, t.end])
    incidence[i] <- max(0, r.end - r.start)
  }

  return(incidence)
}
```

`max(0, ...)` guards against small negative values from population rescaling.

**NegBin log-likelihood on cases**:
```r
ll.cases <- sum(dnbinom(
  x    = casedata$cases,
  mu   = pmax(1e-6, phi * incidence.pred),
  size = kappa,
  log  = TRUE
))
```

**Implementation**: Complete. The following new files implement the joint fit:

- `GetAgeSpecificAnnualIncidence(res, trans, epi.state, year.sim, casedata)` — appends `pred.incidence` to casedata using ΔR from the R compartment
- `GetPredictedMeaslesSerologyAndCases(serodata, casedata, ...)` — runs one simulation, returns `list(serodata=..., casedata=...)` with both predictions
- `LogLikMeaslesCases(casedata, phi, kappa, eps)` — NegBin LL given `pred.incidence`
- `LogLikMeaslesSerologyAndCases(R0, rho, scale.sia, phi, kappa, serodata, casedata, ...)` — calls `GetPredictedMeaslesSerologyAndCases` once then computes both LLs
- `NegLogLikMeaslesSerologyAndCases.transformed(theta, serodata, casedata, ...)` — negative LL on transformed scale
- `FitMeaslesSerologyAndCases(serodata, casedata, ...)` — top-level joint MLE fit
- `TransformMeaslesSerologyParameters(theta, fix.rho, include.case.params)` — extended to handle case params; when `include.case.params = TRUE`: `phi = plogis(theta[n.sero+1])`, `kappa = exp(theta[n.sero+2])`

**Parameter layout** — all combinations of `fix.R0` and `fix.rho`:

| `fix.R0` | `fix.rho` | Free parameters | theta length (with case params) |
|---|---|---|---|
| `NA` (default) | `1` (default) | R0, scale.sia, phi, kappa | 4 |
| a number | `1` | scale.sia, phi, kappa | 3 |
| `NA` | `NA` | R0, scale.sia, rho, phi, kappa | 5 |
| a number | `NA` | scale.sia, rho, phi, kappa | 4 |

theta elements always appear in this order, with fixed parameters omitted: `log(R0)`, `qlogis(scale.sia)`, `atanh(rho)`, `qlogis(phi)`, `log(kappa)`.

**par.init defaults**: `R0 = 12`, `scale.sia = 0.5`, `phi = 0.1` (10% reporting), `kappa = 1` — **do not use phi = 0.1 for real Zambia data**; use `phi.init = sum(cases)/sum(pred.incidence)` ≈ 0.01 instead (see degeneracy note above)

---

**Usage** — 4-parameter joint fit (rho fixed at 1):

```r
fit <- FitMeaslesSerologyAndCases(
  serodata         = serodata,
  casedata         = casedata,   # tibble: age, year, cases; exclude 2006/2010/2011/2026
  setup            = setup,
  year             = 1980,
  t.max            = 45,
  fix.rho          = 1,
  age.classes      = c(1:60, seq(72, 1212, 12)),
  age0is6to11monly = TRUE,
  hessian          = TRUE
)
fit$par.natural      # $R0, $scale.sia, $rho (=1), $phi, $kappa
fit$predictions$casedata  # pred.incidence + pred.cases columns appended
```

The simulation runs once per optimizer iteration (no double-counting).
`fit$predictions$casedata$pred.cases = phi * pred.incidence`.

---

### Strategy 2 — Informative prior on R0 ✅ IMPLEMENTED

Sub-Saharan Africa measles literature consistently puts R0 at 12–18. If the Zambia data genuinely can't distinguish R0 = 8 from R0 = 15, place a log-normal prior centered at ~14. This is Bayesian regularization — borrowing from what's known elsewhere when local data are weakly informative. Requires MCMC to propagate uncertainty properly (MAP estimate with a prior is straightforward but underestimates uncertainty).

**Implementation**: `prior.R0.meanlog = log(14)` and `prior.R0.sdlog = 0.4` arguments added to both `FitMeaslesSerology()` and `NegLogLikMeaslesSerology.transformed()`. When `fix.R0 = NA`, the prior is applied:

```r
ll.prior <- if (is.na(fix.R0)) {
  dlnorm(pars$R0, meanlog = prior.R0.meanlog, sdlog = prior.R0.sdlog, log = TRUE)
} else {
  0
}
return(-(ll + ll.prior))
```

The prior drops out automatically when `fix.R0` is a number (profile likelihood), so profile fits are unaffected. `sdlog = 0.4` gives a 95% prior interval of roughly 7–28. This matches the prior already in `FitMeaslesSerologyAndCases()`, making both calibration approaches consistent — both are now MAP when R0 is free, both are pure conditional MLE when R0 is fixed.

---

### Strategy 3 — Profile R0 over joint serology + cases fit ✅ IMPLEMENTED

Fix R0 at a grid of plausible values and fit (scale.sia, phi, kappa) at each. If the profile likelihood is nearly flat, R0 is weakly identified and a biologically plausible value can be chosen. If there is a clear peak, it is found without assuming an answer.

`fix.R0` is now supported in `FitMeaslesSerologyAndCases` (and `TransformMeaslesSerologyParameters`). When `fix.R0` is a number, R0 drops out of theta and the R0 prior is not applied — giving a pure 3-parameter conditional fit.

**Current implementation**: serology-only profile (no case data). Each fit is 1-parameter (scale.sia only with fix.rho=1), so convergence is fast — typically under 30 iterations each. Note: `sink()` inside forked processes breaks R's connection table — do not use it inside `mclapply`.

```r
library(parallel)

R0.grid <- c(7, 9, 11, 13, 15, 18)

profile <- mclapply(R0.grid, function(r0) {
  fit <- FitMeaslesSerology(
    serodata         = serodata,
    setup            = setup,
    year             = 1980,
    t.max            = t.max,
    fix.rho          = 1,
    fix.R0           = r0,
    age.classes      = c(1:60, seq(72, 1212, 12)),
    age0is6to11monly = TRUE
  )
  list(R0 = r0, scale.sia = fit$par.natural$scale.sia,
       logLik = fit$logLik, converged = fit$convergence == 0)
}, mc.cores = min(length(R0.grid), parallel::detectCores() - 1))

do.call(rbind, lapply(profile, as.data.frame))
```

`mclapply` uses fork-based parallelism on Mac — all sourced functions and objects are inherited automatically, no export needed. Iteration output will not appear in the console during the run. If any element comes back as `try-error` or `NULL`, inspect with `profile[[i]]`.

The profile table (logLik vs R0) shows how much the serology can constrain R0. A flat profile confirms non-identifiability; a peak indicates where the data prefers.

### Profile results (real Zambia data, all ages 2016 + 2024)

```
 R0  scale.sia   logLik  converged
  7  0.164      -543.2   TRUE
  9  0.255      -580.3   TRUE
 11  0.279      -693.2   TRUE
 13  0.320      -833.2   TRUE
 15  0.372      -986.6   TRUE
 18  0.449     -1228.1   TRUE
```

The profile is **strongly monotonically decreasing** — the MLE is off the left edge of the grid, below R0=7. The serology data is informative about R0 but points to a biologically implausible value (R0<7 for measles in sub-Saharan Africa). This is strong evidence of systematic bias, most likely the 2016 assay sensitivity issue.

Note: scale.sia increases with R0 — at biologically plausible R0=12–15, scale.sia≈0.32–0.37 (SIAs achieving ~32–37% of reported coverage in terms of immunological effect).

### R0-fixed profile fits (R0 = 12, 14, 16) with visualization

Following the profile, three fixed-R0 fits were run using `FitMeaslesSerology` with `fix.R0` and `fix.rho = 1` (1-parameter: scale.sia only):

```
 R0   scale.sia   saved file
 12   ~0.28       fit.R0.12.RData
 14    0.346      fit.R0.14.RData
 16   ~0.41       fit.R0.16.RData
```

These are used both as alternative R0 assumptions and to construct a **R0 sensitivity band** on bar plots (see below).

---

## Visualization: predicted vs. observed seroprevalence bar plots

Bar plots compare observed seroprevalence (grey bars + 95% exact binomial CI) against model-predicted seroprevalence (red line/points), aggregated to 5-year age bins. Separate plots are made per survey year and per assumed R0.

### Aggregation to 5-year bins

```r
agg.pred5 <- function(p) {
  p %>%
    mutate(
      survey.year = round(1980 + survey.time.point / 24),
      age.grp     = floor(age.bin.lower / 5) * 5
    ) %>%
    group_by(survey.year, age.grp) %>%
    summarise(pred.prop = sum(pred.imm.pop) / sum(pred.pop), .groups = "drop")
}
```

Use `round()` on `survey.year` to avoid floating-point join failures.

### Uncertainty band: R0 sensitivity range

The Laplace (Hessian) CI on scale.sia is too narrow to be visible — once R0 is fixed at 14, uncertainty in scale.sia has almost no leverage on predicted seroprevalence (most of the immune fraction at older ages is determined by R0 via natural infection, not SIA coverage). The CI difference between lo/hi across age bins is ~0.01.

A more informative and interpretable uncertainty display is the **R0 sensitivity band**: predictions at R0=12 and R0=16 (each at their MLE scale.sia) bracket the R0=14 prediction. This is what appears as the red ribbon in the plots.

```r
pred.ci <- dplyr::left_join(
  agg.pred5(fit.R0.12$predictions) %>% dplyr::rename(pred.prop.lo = pred.prop),
  agg.pred5(fit.R0.16$predictions) %>% dplyr::rename(pred.prop.hi = pred.prop),
  by = c("survey.year", "age.grp")
)
```

### Plot function

```r
make.seroprev.plot <- function(data, title.str) {
  ggplot(data, aes(x = factor(age.label, levels = unique(age.label)))) +
    geom_col(aes(y = obs.prop), fill = "grey70", width = 0.6) +
    geom_errorbar(aes(ymin = ci.lo, ymax = ci.hi), width = 0.25, colour = "grey30") +
    geom_ribbon(aes(ymin = pred.prop.lo, ymax = pred.prop.hi, group = 1),
                fill = "red", alpha = 0.2) +
    geom_point(aes(y = pred.prop), colour = "red", size = 2.5) +
    geom_line(aes(y = pred.prop, group = 1), colour = "red", linewidth = 0.8) +
    scale_y_continuous(limits = c(0, 1),
                       labels = function(x) paste0(round(x * 100), "%")) +
    labs(x = "Age group (years)", y = "Seroprevalence", title = title.str) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}
```

The ribbon (`geom_ribbon`) uses `pred.prop.lo`/`pred.prop.hi` from the R0 sensitivity band. If these are NA (no ribbon requested), ggplot2 silently skips the ribbon layer.

### Plots produced

- **R0 = 14** (central estimate): 2016 survey, 2024 survey — red ribbon = R0 12–16 band
- **R0 = 12** (lower bound): 2016 survey, 2024 survey — no ribbon (single-scenario)
- **R0 = 16** (upper bound): 2016 survey, 2024 survey — no ribbon (single-scenario)

Six plots total; generated by `make.seroprev.plot()` + `dplyr::filter(pred.long, survey.year == YYYY)`.

---

### Assay sensitivity analysis: how much correction shifts R0 to ~14?

Since the profile points to R0<7, we asked: what 2016 assay sensitivity would make R0≈14 the MLE? The correction inflates only the 2016 n_positive values by 1/s (capped at n_tested), leaving the 2024 data unchanged. FitMeaslesSerology is then run with R0 free.

```r
sensitivity.grid <- c(0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
tp.2016  <- (2016 - 1980) * 24
idx.2016 <- serodata$survey.time.point == tp.2016

sens.results <- mclapply(sensitivity.grid, function(s) {
  serodata.adj <- serodata
  serodata.adj$n_positive[idx.2016] <- pmin(
    round(serodata$n_positive[idx.2016] / s),
    serodata$n_tested[idx.2016]
  )
  fit <- FitMeaslesSerology(
    serodata = serodata.adj, setup = setup,
    year = 1980, t.max = t.max, fix.rho = 1,
    age.classes = c(1:60, seq(72, 1212, 12)), age0is6to11monly = TRUE
  )
  list(sensitivity = s, R0 = fit$par.natural$R0,
       scale.sia = fit$par.natural$scale.sia,
       logLik = fit$logLik, converged = fit$convergence == 0)
}, mc.cores = min(length(sensitivity.grid), parallel::detectCores() - 1))
```

Results:

```
 sensitivity    R0    scale.sia   logLik   converged
 0.5          82.4    0.570      -438.98   TRUE
 0.6          81.2    0.547      -410.89   TRUE
 0.7          80.4    0.530      -398.83   TRUE
 0.8          72.4    0.660      -421.27   TRUE
 0.9          14.7    0.518      -593.39   TRUE
 1.0           7.4    0.182      -540.74   TRUE
```

**Key finding**: a sensitivity of ~0.9 (10% correction) shifts the MLE to R0≈14.7 — biologically plausible. Only a 10% reduction in 2016 assay sensitivity is needed; this is within normal assay-to-assay variability.

Below s=0.9, R0 blows up to 72–82: over-correction pushes many cells to 100% seroprevalence (n_positive capped at n_tested), making R0 unidentifiable. `FitMeaslesSerology` has no R0 prior, so it runs away.

**Important**: logLik values are not comparable across sensitivity rows — each row fits a different (adjusted) dataset. They cannot be used to select the best sensitivity value.

**Interpretation for reporting**: *"The low R0 estimate (≈7) is consistent with a 2016 assay sensitivity of ~90% relative to the 2024 assay — a difference that would require independent assay validation to confirm."*

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

## ~~⚠️ BEFORE MERGING TO MAIN — backward compatibility fix required~~ ✅ RESOLVED

The `experiment.updatedemog.vaccinationchange.vaccinationcorrelation` S4 class uses `MR1SIAcorrelation` and `MR2SIAcorrelation` as **numeric** slots (storing the rho value 0–1). The existing `vaccinationlimitations` class keeps them as **boolean** (unchanged). Backward compatibility is preserved by wrapping the assignment in `EX.Country.part2.R` with `as.numeric()`:

```r
EX@MR1SIAcorrelation <- as.numeric(MR1SIAcorrelation)
EX@MR2SIAcorrelation <- as.numeric(MR2SIAcorrelation)
```

`as.numeric(TRUE)` → `1`, `as.numeric(FALSE)` → `0`, numeric values pass through unchanged. Old code passing `TRUE`/`FALSE` continues to work.

**Architecture decision**: `MR1SIAcorrelation` and `MR2SIAcorrelation` remain in the `vaccinationlimitations` class as boolean guards (always `FALSE` in practice — if `TRUE` with SIAinacc or SIAinefficient, the run method stops immediately). The `vaccinationcorrelation` class owns the numeric rho values. Slots were not removed from `vaccinationlimitations` to preserve backward compatibility with any code that constructs that class directly.

---

## Calibration limitations (reviewer notes)

These limitations apply to the current joint serology + age-stratified cases MAP fit.

**1. Deterministic model cannot reproduce outbreak dynamics.**
`is.stochastic = FALSE` predicts a smooth endemic trajectory. Outbreak peaks (2010/11, 2024/25) cannot be reproduced by any choice of parameters — the optimizer distorts R0 trying to match averages of structurally incompatible data. Mitigated by excluding outbreak years from `casedata`.

**2. Only two cross-sectional seroprevalence surveys.**
Serology is from 2016 and 2024 only. Both surveys are used in full. There is no longitudinal cohort tracking, so the age profile of cumulative immunity is the only seroprevalence signal — there is no direct observation of acquisition rates between surveys or within birth cohorts.

**3. Potential assay sensitivity difference between surveys.**
The 2016 survey may use a less sensitive assay than the 2024 survey. If 2016 seroprevalence is systematically suppressed relative to true immunity, the optimizer will try to explain the low 2016 values with lower transmission or coverage, which biases R0 downward and scale.sia upward. This cannot be corrected without independent assay calibration data.

**4. Suspected cases, not confirmed.**
`casedata` contains suspected measles cases. Reporting fraction φ and overdispersion κ partially absorb misclassification, but if non-measles febrile illness has strong seasonal or age structure, it will appear as a structured residual that the model interprets as signal.

**5. NegBin overdispersion absorbs rather than explains variance.**
The NegBin likelihood handles outbreak-year and misclassification variance by inflating κ. This is statistically valid but means the case likelihood can become nearly uninformative (κ → small) if the model fit is poor — equivalent to the case data providing little constraint on parameters. A floor on κ (currently 0.1) prevents complete collapse but is an ad hoc fix, not a modelling solution.

**6. R0–φ non-identifiability.**
As R0 increases, predicted incidence increases, and φ decreases to compensate (φ × incidence stays approximately constant). The log-normal MAP prior on R0 (centered at log(14), sdlog = 0.4) is the only constraint breaking this ridge. Without it, the optimizer finds R0 = 1,000, φ = 0.00001 as equally good as R0 = 14, φ = 0.01. The prior width is judgement-based, not data-derived.

**7. Fixed parameters not estimated.**
MCV1 × 0.85, MCV2 × 0.90, VE = 1, generation time = 0.5 months, seasonal amplitude = 0.15, and all vaccination age windows are fixed. These have been informed by the serology-only fitting and biological priors but are not jointly estimated. Misspecification propagates silently to R0 and scale.sia.

**8. Local optima with Nelder-Mead.**
Nelder-Mead makes no global-optimality guarantees. The 4–5 dimensional joint likelihood surface has a shallow ridge (R0–φ) and a floor constraint (κ). Results should be cross-checked with multiple starting values or a grid search over R0 to confirm the reported optimum is not a local one.

**9. No formal uncertainty quantification.**
MAP point estimates are reported. The Laplace (Hessian) approximation is available (`hessian = TRUE`) but the R0–φ ridge and the κ floor make the likelihood non-Gaussian near the optimum. Intervals from the Hessian will understate uncertainty on φ and κ. Full Bayesian MCMC is the appropriate next step.

**10. Single country, single calibration.**
All parameter choices (prior center, coverage scaling, excluded years) are Zambia-specific. Applying this framework to another country requires re-diagnosing phi/kappa starting values, reviewing excluded years, and reconsidering the R0 prior.

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
