# AI Session Log — MRTransmissionModel

## Branch: outbreak-response
## Goal: Revise model to conduct outbreak response

---

## ⚠️ Pre-existing bugs to fix (found during model-calibration code review)

Both bugs are in the `school.spatial` run method in `R/run.R`. They exist on `main` and on `model-calibration` — not introduced by either branch.

**Bug 1 — Wrong slot name (critical): `R/run.R` line 1171**
```r
obj.prob.vsucc = exper@obj.prob.vsucc.schoolvacc  # slot does not exist
```
Should be:
```r
obj.prob.vsucc = exper@obj.prob.vsucc  # inherited from school parent class
```
Will throw a runtime error if the `experiment.updatedemog.vaccinationchange.school.spatial` run method is ever executed.

**Bug 2 — Vector vs. element NA check (moderate): `R/run.R` lines 1246, 1256, 1263**
```r
is.na(index.school.vacc)    # checks whole vector — always returns a vector, not scalar
```
Should be:
```r
is.na(index.school.vacc[t]) # checks element t
```
Breaks the if-else branch logic controlling when school vaccination is applied each time step.

---

## Origial Thoughts from Amy
We want to modify the model to allow outbreak response "investigations" and "responses".  The goal is to conduct scenario modeling to evaluate the impact of outbreak response on measles case burden.  I don't have deaths in the model, so we can't go that far.  

The dimensions or characteristics of each scenario we want to be able to incorporate include:
1. delay in case detection
2. delay in response 
together 1 and 2 basically equate to total delay (k time steps). I envision in the model, this means that we will at the number of infections in the model k time steps back
3. outbreak response threshold modeled as either total number of infections or incidence (i.e., infections / population) that will trigger a response.  
so based on looking k time steps back, the threshold tells us whether or not we will conduct a response.
4. outbreak response vaccination campaign, that allows different age targets and coverage values.
so if an outbreak response is triggered per 1,2,and 3, then #4 tells of the response

What I am still thinking over is a scenario dimension on proportion of suspected cases tested, but this requires adding observational noise to infection that is not currently included in the model. If we included this dimention, we could add an additional option to outbreak response of futher investigation that conducted more testing. But, adding this observational process is not trivial.  Techincally we could easily model rubella with the same model, but Zambia has introduced rubella vaccine so it is unlikely to be resulting in a bunch of noise.  "Maculopapularskin rash may be caused by numerous viruses (e.g.rubella virus, parvovirus B19, human herpesviruses 6 and 7, entero-and adenoviruses, Epstein-Barr virus, cytomegalovirus, coxsack-ievirus, as well as, in travellers returning from endemic areas, dengue virus, Chikungunya-, Zika-, West Nile-, Ross River-, and Sindbis virus) and the bacterium Streptococcus pyogenes, but also by allergies and drugs".  I do have the proportion IgM positive over age and year and over month and year in Zambia, maybe we can do something with that?  But, know it is very noisy often b/c some months or ages have <5 suspected cases tested so the proportion test positive can jump from 0 and 1 easily.

---

## S4 Class Design: experiment.updatedemog.vaccinationchange.vaccinationcorrelation.outbreakresponse

The OBR experiment class extends `experiment.updatedemog.vaccinationchange.vaccinationcorrelation` with the following slots. **Note:** `or.total.delay` is not supplied directly by the user — it is computed inside `EX.Country.part2.OR()` as `or.confirmation.delay + or.response.delay` and then stored in this slot.

### Trigger slots
| Slot | Type | Description |
|------|------|-------------|
| `or.trigger.mode` | character | `"I_scaled"` (true infections × reporting rate) or `"confirmed"` (full observational model) |
| `or.reporting.rate` | numeric | Reporting rate r ∈ (0,1]; scalar or vector of length n.age. Used in both modes. Default 1 (no scaling). |
| `or.total.delay` | numeric | Total delay in time steps; set to `or.confirmation.delay + or.response.delay` by helper |
| `or.trigger.window` | numeric | Number of steps over which cumulative metric is summed |
| `or.n.confirmations.target` | numeric | Number of (confirmed) cases at which the trigger fires; default 5 |
| `or.trigger.age.lower` | numeric | Lower age bound (months) for trigger surveillance; `NA` = all ages |
| `or.trigger.age.upper` | numeric | Upper age bound (months) for trigger surveillance; `NA` = all ages |
| `or.start.timestep` | numeric | Earliest time step at which OBR can trigger (1-indexed); default 1 |

### Observational model slots (required for `or.trigger.mode = "confirmed"`)
| Slot | Type | Description |
|------|------|-------------|
| `or.non.meas.cases.by.age.month` | matrix | Non-measles suspected cases; rows = model age classes, columns = 12 calendar months |
| `or.Se` | numeric | Diagnostic test sensitivity ∈ [0,1] |
| `or.Sp` | numeric | Diagnostic test specificity ∈ [0,1] |

### Response slots
| Slot | Type | Description |
|------|------|-------------|
| `or.response.delay` | numeric | Time steps from trigger to campaign delivery; case age distribution is accumulated over this window |
| `or.vacc.age.lower` | numeric | Lower age bound (months) for response campaign; default 0 |
| `or.vacc.age.upper` | numeric | Upper age bound (months) for response campaign; `NA` = use `or.vacc.agedist.percentile` |
| `or.vacc.agedist.percentile` | numeric | Percentile (0–1) of cumulative case age CDF used as campaign upper age bound when `or.vacc.age.upper` is `NA` |
| `or.vacc.coverage` | numeric | Coverage of response campaign (0–1); set to 0 to disable OBR |
| `or.min.interval` | numeric | Minimum steps between successive OBR triggers |

### `EX.Country.part2.OR()` helper — user-facing delay arguments

The helper accepts `or.confirmation.delay` and `or.response.delay` separately (both required, no defaults):

```r
or.total.delay <- or.confirmation.delay + or.response.delay
```

`or.total.delay` is then stored in the experiment object. `or.confirmation.delay` is not stored as a slot — only the computed total and `or.response.delay` are kept.

### Result class

```r
setClass(
  "sim.results.MSIRV.update.demog.vaccine.change.outbreakresponse",
  contains = "sim.results.MSIRV.update.demog.vaccine.change",
  slots = list(
    or.times        = "numeric",  # length numTimeSteps; 1 when OBR campaign fired
    obs.TP          = "matrix",   # age x timestep; true measles cases tested & test-positive
    obs.FN_test     = "matrix",   # age x timestep; true measles cases tested & test-negative
    obs.FP_test     = "matrix",   # age x timestep; non-measles cases tested & test-positive
    obs.TN          = "matrix",   # age x timestep; non-measles cases tested & test-negative
    obs.TP_clinical = "matrix",   # age x timestep; true measles cases untested, clinical/epi-linked
    obs.FP_clinical = "matrix"    # age x timestep; non-measles cases untested, clinical/epi-linked
  )
)
```

The obs matrices are always allocated (even when no observational model is used); they contain zeros when `or.trigger.mode = "I_scaled"` and `or.non.meas.cases.by.age.month` is not supplied.

### Trigger logic

Two trigger modes:

**`"I_scaled"` (default):** metric = sum of `I * or.reporting.rate` over the surveillance window:
```r
t.end   <- t - exper@or.total.delay
t.start <- t.end - exper@or.trigger.window + 1
metric  <- sum(rc[trigger.i.inds, t.start:t.end] * or.rep.rate)
```
Trigger fires when `metric >= or.n.confirmations.target`.

**`"confirmed"`:** uses the full observational model (Se/Sp/non-measles background) to compute `confirmed.trigger[t]` each step. The same lookback window is then summed over `confirmed.trigger` rather than raw I.

The earliest time step at which OBR can trigger:
```r
or.min.t <- max(exper@or.total.delay + exper@or.trigger.window, exper@or.start.timestep)
```

### Dynamic campaign age targeting

When `or.vacc.age.upper` is `NA`, the campaign upper age is derived from the case age distribution accumulated during the `or.response.delay` window:
```r
cdf       <- cumsum(pending.case.by.age) / sum(pending.case.by.age)
pct.idx   <- which(cdf >= exper@or.vacc.agedist.percentile)[1]
vacc.upper <- age.classes[pct.idx]
```
In `"I_scaled"` mode the accumulation uses `I * r + non-measles background`; in `"confirmed"` mode it uses the confirmed case counts from the obs matrices.

### Integration into package ✅ COMPLETE

| File | Change |
|------|--------|
| `R/setClasses.R` | Added `experiment.updatedemog.vaccinationchange.vaccinationcorrelation.outbreakresponse` and `sim.results.MSIRV.update.demog.vaccine.change.outbreakresponse` class definitions |
| `R/run.R` | Added `setMethod("run", "experiment...outbreakresponse", ...)` after the `vaccinationcorrelation` run method |
| `R/EX.Country.part2.OR.R` | `EX.Country.part2.OR()` helper function |
| `NAMESPACE` | Exports for class and helper |
| `DESCRIPTION` | Collate entry for `EX.Country.part2.OR.R` |
| `branch_building/OutbreakResponse.R` | Scenario construction example only |

Notes on `EX.Country.part2.OR()`:
- `SIAinacc`, `SIAinefficient`, and `prop.inacc` are excluded — incompatible with the `vaccinationcorrelation` parent class
- `MR1MR2correlation`, `MR1SIAcorrelation`, `MR2SIAcorrelation` default to full correlation (`TRUE`, `1`, `1`)
- `MR1SIAcorrelation` and `MR2SIAcorrelation` are wrapped in `as.numeric()` for backward compatibility
- `or.vacc.coverage = 0` pattern for "no OBR" baseline — both scenarios use `EX.Country.part2.OR()` for consistency

---

## Observational Model: Suspected Cases and IgM Positivity

### The core problem

The model predicts true measles infections. Surveillance data are suspected cases — a rash-febrile illness that could be measles or any of a dozen other causes (rubella, parvovirus B19, human herpesviruses 6/7, enteroviruses, adenoviruses, EBV, CMV, coxsackievirus, Streptococcus pyogenes, allergies, drugs, and in travellers: dengue, Chikungunya, Zika, West Nile, etc.).

The key identity is:

```
P(IgM+ | suspected) = measles_incidence / (measles_incidence + background_rash)
```

The model only supplies the numerator.

**Zambia data context**: Proportion IgM+ is available by age/year and by month/year. It is very noisy because some months or age strata have <5 suspected cases tested, causing the proportion to jump between 0 and 1.

### Options (ranked by complexity)

**Option 1 — Use IgM-confirmed cases as calibration target (simplest)**

Multiply suspected cases × IgM+ proportion to get approximate confirmed counts, then aggregate to annual totals to smooth the noise. This sidesteps the background rash problem entirely. The reporting fraction `φ` then absorbs: (a) clinical presentation rate, (b) presenting as a suspected case, and (c) getting tested and confirmed. `φ` is smaller, but the signal is cleaner. The noisy monthly proportions are only a problem if used as a calibration target; annual aggregation removes most of the noise.

**Option 2 — Background rash as a nuisance parameter (moderate)**

Model `total_suspected[t] = measles[t] + B`, where `B` is a constant (or slowly varying) background rash rate. Then `P(IgM+) = measles[t] / (measles[t] + B)`. `B` can be estimated from inter-epidemic troughs where measles incidence is near zero. Adds one parameter and makes the IgM+ proportion directly modelable.

**Option 3 — Beta-binomial likelihood for monthly noise**

A beta-binomial likelihood absorbs over-dispersion from months with <5 tests without requiring aggregation. Most useful if monthly data is the calibration target. For case calibration, annual confirmed counts are probably simpler; reserve the monthly data for calibrating `seasonal.amp`.

**Option 4 — Time-varying φ for the "proportion tested" OBR dimension**

`φ` becomes a two-level step function: `φ_low` at baseline (routine surveillance), `φ_high` during an active OBR investigation. The switch timing is a scenario parameter tied to the OBR trigger. This is tractable without a full observation model — parameterize as `(φ_baseline, φ_investigation)` and add the switch as a scenario dimension.

### Current status

**The observational model (Option 2 / Option 3 hybrid) is implemented** as part of the `"confirmed"` trigger mode in `R/run.R`. At each time step it:

1. Computes total suspected cases = true reported measles (`I * r`) + non-measles background (from `or.non.meas.cases.by.age.month`)
2. Derives an expected IgM positivity rate using Se and Sp
3. Allocates tests proportionally across age classes up to `or.n.confirmations.target / pos.rate` tests
4. Fills `obs.TP`, `obs.FN_test`, `obs.FP_test`, `obs.TN`, `obs.TP_clinical`, `obs.FP_clinical` matrices

The `"I_scaled"` mode (default) skips the full obs model and simply scales true infections by `or.reporting.rate` — this is equivalent to Option 1.

**For calibration**: using IgM-confirmed counts (suspected × proportion IgM+, aggregated annually) as a calibration target and letting `φ` absorb under-detection remains the recommended starting point (Option 1 / `"I_scaled"` mode). The `"confirmed"` mode is available for scenarios that explicitly model the surveillance/testing process as a trigger mechanism.

The biggest contributors to background rash in Zambia are likely enteroviruses and parvovirus B19. Rubella noise is probably small given Zambia's vaccine introduction. **Background rash as a constant nuisance `B` is probably the right framing** for `or.non.meas.cases.by.age.month` rather than trying to model each pathogen separately.

---

## Time Step: Current Design and Considerations for Change

### Current step size

The model runs at `step.size = 1/24` years, i.e., **2-week (0.5-month) steps** with 24 steps per year. This is set in the experiment object and propagates into the WAIFW matrix and age-survival matrix.

### Why finer resolution might matter for OBR

With 2-week steps, all delay and window parameters (`or.total.delay`, `or.trigger.window`, `or.min.interval`) must be integer multiples of 2 weeks. A 3-week detection lag, for example, cannot be represented exactly — it would round to 2 or 4 weeks. Weekly steps would halve that rounding error and allow more realistic delay scenarios.

### Why changing the step size is non-trivial

The core blocker is in `src/MRfxns.c`. In the function `calc_phi_and_update_tran`, the I→I entry of the transition matrix is **hardcoded to zero**:

```c
/* line ~240 */
REAL(tran_matrix)[ind] = 0;  /* I -> I = 0: everyone recovers in one step */
```

This means regardless of step size, every infected individual recovers in exactly one time step. With the current 2-week step this is correct (measles infectious period ≈ 2 weeks). With a weekly step, the I→I entry should be approximately **0.5** (one-period survival probability for a 2-week infectious period).

### What would need to change for weekly steps

1. **`calc_phi_and_update_tran` in `src/MRfxns.c`**: Add a `recovery_prob` argument. Replace the hardcoded `0` with `1 - recovery_prob` for the I→I diagonal, and use `recovery_prob` for the I→R entries. For weekly steps, `recovery_prob = 0.5`.

2. **WAIFW scaling**: The WAIFW matrix encodes per-step transmission. When halving the step size, each entry must be scaled by `new_step / old_step = 0.5` to preserve the same per-period force of infection. This can be done in R before the matrix is passed to C.

3. **age.surv.matrix**: The aging and survival entries are also per-step. These must be recomputed or rescaled for the new step size. Currently this is handled in `GetWAIFW` or the experiment setup functions.

### Deterministic workaround (no C change required)

The original R implementation of the transition is preserved in commented-out code in `R/SIRIDTran.R` (lines ~32–41). For the **deterministic model only**, it would be possible to intercept the transition matrix after the C call and manually reset the I→I diagonal to the correct value before applying `tran.matrix %*% state`. This avoids modifying C but only works for the deterministic path; the stochastic path calls a separate C function (`do_ID_transition_SIR_stochastic_moves_cl`) and cannot be patched this way.

### Pragmatic assessment

**2-week resolution is likely sufficient for OBR scenario modeling.** Most public health response decisions (detection, mobilization, campaign implementation) do not operate at sub-2-week precision in practice. If we need to compare scenarios that differ by a 1-week delay, we can approximate by comparing 2-week vs. 4-week delay scenarios as bounds. Changing the step size would require C code modification, rebuilding the package, and re-validating the model — a non-trivial investment for marginal gain in temporal precision.