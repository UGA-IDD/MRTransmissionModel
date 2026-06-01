# AI Session Log — MRTransmissionModel

## Branch: model-calibration
## Goal: Build and speed up measles serology calibration

---

## Background

This package (`MRTransmissionModel`) is an MSIRV (Maternal/Susceptible/Infected/Recovered/Vaccinated) compartmental measles and rubella transmission model. A branch called `model-calibration` was started to add functionality for calibrating the model to observed serological data.

The calibration fits two parameters:
- `R0` — basic reproduction number (optimized on log scale, i.e. `theta[1] = log(R0)`)
- `rho` — Pearson correlation between RI (routine immunization) and SIA doses, in (-1, 1) (optimized on atanh scale, i.e. `theta[2] = atanh(rho)`)

The entry point for fitting is `FitMeaslesSerology()`. For Bayesian work, the key likelihood function is `LogLikMeaslesSerology()`.

---

## Files added/modified on this branch

### New R files
| File | Purpose |
|------|---------|
| `R/GetPredictedMeaslesSerology.R` | Runs Part 1 + Part 2 simulation, returns predicted seroprevalence aligned to observed data rows |
| `R/GetMeaslesSeroprevalence.per.TimePoint.R` | Extracts immune + total population by age group at requested time points only |
| `R/AggregateMeaslesSeroprevalenceByAgeBins.R` | Aggregates age-year seroprevalence into user-specified age bins |
| `R/LogLikMeaslesSerology.R` | `LogLikMeaslesSerology()` — binomial log-likelihood |
| `R/TransformMeaslesSerologyParameters.R` | `TransformMeaslesSerologyParameters()` — maps `(log(R0), atanh(rho))` → `(R0, rho)` |
| `R/LogLikMeaslesSerology.transformed.R` | `LogLikMeaslesSerology.transformed()` — convenience wrapper on transformed scale |
| `R/NegLogLikMeaslesSerology.transformed.R` | `NegLogLikMeaslesSerology.transformed()` — negative log-lik for use with `optim()` |
| `R/FitMeaslesSerology.R` | `FitMeaslesSerology()` — full MLE fit with optional grid search |

### Modified R files
| File | What changed |
|------|-------------|
| `R/GetNumber.per.AgeYear.R` | Vectorized the inner double for-loop with `colSums(matrix(...))` — was 240 R iterations, now 1 call |
| `R/EX.Country.part2.R` | Routing logic updated to use `vaccinationcorrelation` class when `MR1MR2correlation=TRUE` or `MR1SIAcorrelation != 0` or `MR2SIAcorrelation != 0` |
| `R/run.R` | New `setMethod("run", "experiment.updatedemog.vaccinationchange.vaccinationcorrelation", ...)` added; intro rate indexing fixed |
| `R/setClasses.R` | New S4 class `experiment.updatedemog.vaccinationchange.vaccinationcorrelation` added |
| `DESCRIPTION` | Added new R files to `Collate`; removed phantom `NegLogLikMeaslesSerology.R` entry |

---

## The rho parameter — what it is and how it works

`rho` is the Pearson correlation between the number of RI (routine immunization) doses and SIA doses a child receives. It replaces the old `sia.scale` multiplicative adjustment.

- `rho = 0`: independence (RI and SIA coverage uncorrelated)
- `rho > 0`: children who receive RI are more likely to also be reached by SIAs
- `rho < 0`: SIAs preferentially reach children missed by RI

The combined vaccination probability at each SIA time step uses the bivariate formula:
```
pvacc = r + s - r*s - rho * sqrt(r*(1-r)*s*(1-s))
```
where `r` = routine coverage and `s` = SIA coverage for that age class.

This is implemented in the `run` method for `experiment.updatedemog.vaccinationchange.vaccinationcorrelation` in `run.R`.

### New S4 class: `experiment.updatedemog.vaccinationchange.vaccinationcorrelation`

Defined in `setClasses.R`. Inherits from `experiment.updatedemog.vaccinationchange`. New slots:
- `MR1MR2correlation` — logical
- `MR1SIAcorrelation` — numeric (the `rho` value for MR1-SIA correlation)
- `MR2SIAcorrelation` — numeric (the `rho` value for MR2-SIA correlation)

`EX.Country.part2()` routes to this class when `MR1MR2correlation=TRUE` OR either SIA correlation is non-zero. With `rho=0`, the vaccinationcorrelation class is still used (because `MR1MR2correlation=TRUE` is always passed from `GetPredictedMeaslesSerology`), but the formula reduces to simple independent union: `r + s - r*s`.

---

## Speedup fixes implemented

### Fix 1 — Post-processing only at observed time points
**File**: `GetMeaslesSeroprevalence.per.TimePoint.R`

Only processes the unique time points present in `serodata$survey.time.point` (typically 2–5), not all ~1,081 simulation time steps. Roughly 360× faster for that function.

### Fix 2 — Vectorized age aggregation
**File**: `GetNumber.per.AgeYear.R`

Replaced a double `for` loop (240 R iterations) with:
```r
pop.per.youngage.year <- colSums(matrix(vec[seq_len(top.one.month.age)], nrow = 12))
```

### Fix 3 — Part 1 result cache
**File**: `GetPredictedMeaslesSerology.R`

Added `.cache = NULL` parameter (pass `new.env(parent = emptyenv())`). The 3 × 20-year Part 1 transient burn-in result is cached by R0. When only `rho` changes between calls (same R0), Part 1 is skipped entirely — most useful during grid search where R0 is held fixed across all `rho` values in a row.

### Fix 4 — Nelder-Mead optimizer + internal cache
**File**: `FitMeaslesSerology.R`

- Default optimizer: `"Nelder-Mead"` (derivative-free; avoids finite-difference gradient evaluations that would otherwise triple the number of likelihood calls)
- Cache created internally and shared across all evaluations within a single `FitMeaslesSerology` call

---

## Age class coarsening — tested and working

Default: `c(1:240, seq(252, 1212, 12))` = 321 classes → state vector 1,605 rows  
Coarse:  `c(1:60, seq(72, 1212, 12))` = 156 classes → state vector 780 rows (~4× speedup on matrix multiply)

**Tested**: `FitMeaslesSerology(..., age.classes=c(1:60, seq(72,1212,12)))` completed in **< 45 minutes** and returned:
```
R0  = 18.07
rho = -0.27
```
The negative `rho` indicates SIAs preferentially reach children missed by routine immunization, which is epidemiologically plausible.

Pass `age.classes` directly to `FitMeaslesSerology` and it propagates all the way through:

`FitMeaslesSerology` → `NegLogLikMeaslesSerology.transformed` → `LogLikMeaslesSerology` → `GetPredictedMeaslesSerology` → `EX.Country.part1` / `EX.Country.part2`

---

## Package loading workflow

**`devtools::load_all(".")` crashes R on this machine — do not use it.**

The workflow that works:

```r
devtools::document()   # regenerates NAMESPACE + man/ from @export tags
devtools::build()      # creates .tar.gz bundle (does NOT install)
library(MRTransmissionModel)   # loads the previously installed version

# Source any new R files not yet in the installed package:
source("R/GetPredictedMeaslesSerology.R")
source("R/GetMeaslesSeroprevalence.per.TimePoint.R")
source("R/AggregateMeaslesSeroprevalenceByAgeBins.R")
source("R/LogLikMeaslesSerology.R")
```

The `source()` calls are needed because `library()` loads the installed package (which was built from an earlier snapshot), while the new R files live only in the source tree. Once you run `devtools::install(".", quick=TRUE)` to rebuild the installed package with the new files, the `source()` calls are no longer needed.

**When to re-install** (i.e. run `devtools::install(".", quick=TRUE)`):
- After any change to C code in `src/` (requires recompilation)
- After adding new R files to the package so they appear in `library()` without manual `source()`

**When `document()` alone is enough**:
- After adding or changing `@export`, `@param`, or other roxygen tags (regenerates `NAMESPACE` and `man/`)

---

## Test data and script

Simulated Zambia data loaded from:
```
/Users/winter/Library/CloudStorage/GoogleDrive-amykwinter@gmail.com/My Drive/Indrajit_and_Amy/Seroprevalence data modelling framework/Get_Zambia_Data/simulated_zambia_data.RData
```

Working test script: `branch_building/ModelCalibration.R`

```r
devtools::document()
devtools::build()
library(MRTransmissionModel)
source("R/GetPredictedMeaslesSerology.R")
source("R/GetMeaslesSeroprevalence.per.TimePoint.R")
source("R/AggregateMeaslesSeroprevalenceByAgeBins.R")
source("R/LogLikMeaslesSerology.R")
setup <- setupCountry.Nov2023(country="Zambia")
year  <- 1980
t.max <- 45

# Single likelihood evaluation
LogLikMeaslesSerology(
  R0 = 16, rho = 0,
  serodata = serodata, setup = setup,
  year = 1980, t.max = t.max,
  age.classes = c(1:60, seq(72, 1212, 12))
)

# Full MLE fit
fit <- FitMeaslesSerology(
  serodata    = serodata,
  setup       = setup,
  year        = 1980,
  t.max       = t.max,
  par.init    = c(log(16), atanh(0)),
  age.classes = c(1:60, seq(72, 1212, 12))
)
fit$par.natural  # R0 = 18.07, rho = -0.27
```

---

## Posterior uncertainty — Normal (Laplace) approximation

The normal approximation gives a posterior distribution over `(R0, rho)` from the Hessian of the negative log-likelihood at the MLE — no additional model runs required.

Pass `hessian=TRUE` to `FitMeaslesSerology`. The return object includes:
- `fit$cov.transformed` — 2×2 covariance matrix on the `(log(R0), atanh(rho))` scale

Sample from the approximate posterior:
```r
samples.transformed <- MASS::mvrnorm(
  n = 2000, mu = fit$par.transformed, Sigma = fit$cov.transformed
)
posterior.samples <- data.frame(
  R0  = exp(samples.transformed[, "log_R0"]),
  rho = tanh(samples.transformed[, "atanh_rho"])
)
```
Full code including summaries and plots is in `branch_building/ModelCalibration.R`.

**Validity check**: if marginal histograms look skewed or the joint scatter is non-elliptical, upgrade to MCMC (see below).

---

## Future work: full Bayesian MCMC ⚠️ REMINDER

**Not yet done — worth doing when time allows.**

The normal approximation is convenient but approximate. Full Metropolis-Hastings MCMC would give exact posterior samples and catch any non-normality (skew, nonlinear R0-rho correlation, etc.).

**Why deferred**: ~25-30 sec per likelihood evaluation → ~35 hours for 5,000 iterations. Feasible overnight but not interactive.

**How to implement**:
- MH on `(log(R0), atanh(rho))`, initialized at the MLE
- Proposal covariance: `2.38²/2 * fit$cov.transformed` (standard tuning)
- Pass a single `.cache` env to every likelihood call throughout the chain
- `adaptMCMC` R package handles adaptive tuning automatically
- Run 3+ chains for convergence diagnostics (R-hat)

---

## Remaining potential improvements

1. **MCMC** — see reminder above.

2. **Parallelise the grid search** — `mapply` in `FitMeaslesSerology` is single-threaded. Could use `parallel::mcmapply` or `future.apply::future_mapply`.

3. **`NegLogLikMeaslesSerology.transformed` iteration counter** — defined inside `local()` so the counter persists across multiple `FitMeaslesSerology` calls in the same R session. Iteration numbers will be wrong on a second call.

4. **`age.classes` default inconsistency** — some function signatures still carry old defaults that differ from the coarse vector used in calibration. These defaults are never used during calibration (caller always passes explicitly) but are misleading.

---

## Key design decisions

- `survey.time.point` in `serodata` must be a **time-step index** (not a calendar year). For `year=1980`, `generation.time=0.5` (24 steps/year): calendar year `y` → index `(y - 1980) * 24`.
- `age.bin.lower` and `age.bin.upper` are in **single-year units** (0, 1, 2, ...), not months.
- `rho` is passed as `MR1SIAcorrelation = rho` and `MR2SIAcorrelation = rho` to `EX.Country.part2()` (same correlation assumed for MR1-SIA and MR2-SIA pairs).
