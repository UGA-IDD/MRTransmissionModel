# AI Session Log — MRTransmissionModel

## Session date: 2026-05-31
## Branch: model-calibration
## Goal: Build and speed up measles serology calibration

---

## Background

This package (`MRTransmissionModel`) is an MSIRV (Maternal/Susceptible/Infected/Recovered/Vaccinated) compartmental measles and rubella transmission model. A branch called `model-calibration` was started to add functionality for calibrating the model to observed serological data.

The calibration fits two parameters:
- `R0` — basic reproduction number (optimized on log scale)
- `sia.scale` — multiplicative adjustment to SIA (supplemental immunization activity) coverage (optimized on logit scale)

The entry point for fitting is `FitMeaslesSerology()`.

---

## Files added/modified this session

### New files (on `model-calibration` branch)
| File | Purpose |
|------|---------|
| `R/GetPredictedMeaslesSerology.R` | Runs Part 1 + Part 2 simulation, returns predicted seroprevalence aligned to observed data rows |
| `R/GetMeaslesSeroprevalence.per.TimePoint.R` | Extracts immune + total population by age group at requested time points only |
| `R/AggregateMeaslesSeroprevalenceByAgeBins.R` | Aggregates age-year seroprevalence into user-specified age bins |
| `R/LogLikMeaslesSerology.R` | Contains `LogLikMeaslesSerology`, `TransformMeaslesSerologyParameters`, `LogLikMeaslesSerology.transformed`, `NegLogLikMeaslesSerology.transformed`, `FitMeaslesSerology` |

### Modified files
| File | What changed |
|------|-------------|
| `R/GetNumber.per.AgeYear.R` | Vectorized the inner double for-loop with `colSums(matrix(...))` — was 240 R iterations, now 1 call |

---

## Speedup fixes implemented

### Fix 1 — Post-processing only at observed time points
**File**: `GetMeaslesSeroprevalence.per.TimePoint.R`

Previously iterated over all ~1,081 simulation time steps. Now only processes the unique time points present in `serodata$survey.time.point` (typically 2–5). Roughly 360× faster for that function.

### Fix 2 — Vectorized age aggregation
**File**: `GetNumber.per.AgeYear.R`

Replaced a double `for` loop (240 R iterations) with:
```r
pop.per.youngage.year <- colSums(matrix(vec[seq_len(top.one.month.age)], nrow = 12))
```

### Fix 3 — Part 1 result cache
**File**: `GetPredictedMeaslesSerology.R`

Added `.cache = NULL` parameter (pass `new.env(parent = emptyenv())`). The 3 × 20-year Part 1 transient burn-in result is cached by R0. When only `sia.scale` changes between calls (same R0), Part 1 is skipped entirely (~57% of compute avoided on cache hits).

### Fix 4 — FitMeaslesSerology rewrite
**File**: `LogLikMeaslesSerology.R`

- Default optimizer changed from `"BFGS"` to `"Nelder-Mead"` (derivative-free; avoids finite-difference gradient evaluations)
- Cache created internally at the start of each `FitMeaslesSerology` call and shared across all evaluations
- Optional `n.grid` parameter for coarse grid search before optimization (see below)

---

## Age class coarsening — tested and working

The default age class vector `c(1:240, seq(252, 1212, 12))` = 321 classes → state vector 1,605 rows. This makes the core matrix multiply (1,605 × 1,605 per time step) the main bottleneck.

A coarser vector `c(1:60, seq(72, 1212, 12))` = 156 classes → state vector 780 rows was tested and confirmed working end-to-end. This gives roughly 4× speedup on the matrix multiply.

**Tested**: `LogLikMeaslesSerology(R0=16, sia.scale=0.8, ..., age.classes=c(1:60, seq(72,1212,12)))` returned `-305.0634` in 36 seconds.

**You can now pass `age.classes` directly to `FitMeaslesSerology`**:
```r
fit <- FitMeaslesSerology(
  serodata    = serodata,
  setup       = setup,
  year        = 1980,
  t.max       = t.max,
  par.init    = c(log(16), qlogis(0.8)),
  age.classes = c(1:60, seq(72, 1212, 12))
)
```

---

## age.classes propagation — confirmed complete

A full audit was done of the call chain:

`FitMeaslesSerology` → `NegLogLikMeaslesSerology.transformed` → `LogLikMeaslesSerology` → `GetPredictedMeaslesSerology` → `EX.Country.part1` → `Get.CountryX.Starting.Pop.MSIRV` → sub-functions

**All call sites already pass `age.classes` explicitly** (no function relies on its default value during calibration). Downstream functions in the `run()` method use `exper@trans@age.class` from the experiment object, which is built from the correct `age.classes`. No edits were needed to achieve propagation — it was already in place.

---

## The Part 1 cache and Nelder-Mead

The Part 1 cache helps most when R0 stays the same between consecutive likelihood evaluations. This happens systematically during grid search (`n.grid > 1`) — sia.scale varies fastest so each R0 value generates multiple cache hits. During Nelder-Mead itself, R0 changes on almost every iteration (Nelder-Mead moves both parameters simultaneously), so cache hits are occasional but not reliable.

**Recommendation**: if calibration is too slow, use `n.grid = 4` or `5` to run a grid search first, then Nelder-Mead starts near the optimum and needs fewer iterations.

```r
fit <- FitMeaslesSerology(
  serodata         = serodata,
  setup            = setup,
  year             = 1980,
  t.max            = t.max,
  age.classes      = c(1:60, seq(72, 1212, 12)),
  n.grid           = 5,
  R0.range         = c(5, 25),
  sia.scale.range  = c(0.3, 1.0)
)
```

---

## Where we left off

A full `FitMeaslesSerology` run was in progress at session end:
```r
fit <- FitMeaslesSerology(
  serodata    = serodata,
  setup       = setup,
  year        = 1980,
  t.max       = t.max,
  par.init    = c(log(16), qlogis(0.8)),
  age.classes = c(1:60, seq(72, 1212, 12))
)
```
It was running and printing iterations (had reached iter 12 when session ended). Each evaluation takes ~36 seconds. Expected runtime: 30–90 minutes.

The simulated Zambia data used for testing is loaded from:
```
/Users/winter/Library/CloudStorage/GoogleDrive-amykwinter@gmail.com/My Drive/Indrajit_and_Amy/Seroprevalence data modelling framework/Get_Zambia_Data/simulated_zambia_data.RData
```
The test script is `Untitled.R` in the project root.

---

## Remaining potential improvements (not yet done)

1. **Parallelise the grid search** — `mapply` in the grid search is single-threaded. Could use `parallel::mcmapply` or `future.apply::future_mapply` to evaluate grid points in parallel.
2. **BFGS with cache** — now that the Part 1 cache makes gradient evaluations cheaper, BFGS might converge faster than Nelder-Mead for smooth likelihood surfaces. Worth benchmarking once the current run completes and you know the optimum.
3. **`age.classes` defaults inconsistency** — several function signatures still have old defaults (`c(1:240, seq(241,720,12))` or similar) that don't match the current standard. These defaults are never used during calibration (because the caller always passes explicitly), but they are misleading. Could be cleaned up.
4. **`NegLogLikMeaslesSerology.transformed` iteration counter** — defined inside `local()` so the counter persists across multiple `FitMeaslesSerology` calls in the same R session. Minor cosmetic issue but iteration numbers will be wrong on a second call.

---

## Key design decisions

- `survey.time.point` in `serodata` must be a **time-step index** (not a calendar year). For `year=1980`, `generation.time=0.5` (24 steps/year), calendar year `y` → index `(y - 1980) * 24`.
- `age.bin.lower` and `age.bin.upper` are in **single-year units** (0, 1, 2, ...), not months.
- The bug `age0is9to12monly` (wrong name) → `age0is9to11monly` was fixed in `LogLikMeaslesSerology`'s call to `GetPredictedMeaslesSerology`.
