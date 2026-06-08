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

The OBR experiment class extends `experiment.updatedemog.vaccinationchange.vaccinationcorrelation` with two groups of new slots.

### Trigger slots
| Slot | Type | Description |
|------|------|-------------|
| `or.total.delay` | numeric | Total delay in time steps (detection + response lag) |
| `or.trigger.window` | numeric | Number of steps over which cumulative I is summed |
| `or.threshold.value` | numeric | Threshold value (count or rate) |
| `or.threshold.type` | character | `"count"` or `"incidence"` |
| `or.trigger.age.lower` | numeric | Lower age bound (months) for trigger surveillance; `NA` = all ages |
| `or.trigger.age.upper` | numeric | Upper age bound (months) for trigger surveillance; `NA` = all ages |

### Response slots
| Slot | Type | Description |
|------|------|-------------|
| `or.vacc.age.lower` | numeric | Lower age bound (months) for response campaign |
| `or.vacc.age.upper` | numeric | Upper age bound (months) for response campaign |
| `or.vacc.coverage` | numeric | Coverage of response campaign (0–1); set to 0 to disable OBR |
| `or.min.interval` | numeric | Minimum steps between successive OBR triggers |
| `or.start.timestep` | numeric | Earliest time step at which OBR can trigger (1-indexed); default 1 |

### Result class

```r
setClass(
  "sim.results.MSIRV.update.demog.vaccine.change.outbreakresponse",
  contains = "sim.results.MSIRV.update.demog.vaccine.change",
  slots = list(or.times = "numeric")
)
```
`or.times` is a vector of length `t.max` with 1 at each step an OBR was triggered, 0 otherwise.

### Trigger metric: cumulative I (not ΔR)

The trigger looks backward in the simulation record `rc`:

```r
t.end   <- t - exper@or.total.delay
t.start <- t.end - exper@or.trigger.window + 1
cum.I   <- sum(rc[trigger.i.inds, t.start:t.end])
```

We use the sum of I compartments (not change in R) because the I→R transition takes one full time step, which would add an additional implicit delay if ΔR were used as the metric.

The earliest time step at which OBR can trigger is:
```r
or.min.t <- max(or.total.delay + or.trigger.window, or.start.timestep)
```
This ensures there is always enough history to look back without indexing out of bounds, and additionally respects the user-specified start point (e.g. `or.start.timestep = (2020 - year) * 24 + 1` to prevent OBR before 2020).

### Integration into package ✅ COMPLETE

All pieces extracted from `branch_building/OutbreakResponse.R` and moved into the package:

| File | Change |
|------|--------|
| `R/setClasses.R` | Added `experiment.updatedemog.vaccinationchange.vaccinationcorrelation.outbreakresponse` and `sim.results.MSIRV.update.demog.vaccine.change.outbreakresponse` class definitions |
| `R/run.R` | Added `setMethod("run", "experiment...outbreakresponse", ...)` after the `vaccinationcorrelation` run method |
| `R/EX.Country.part2.OR.R` | New file — `EX.Country.part2.OR()` helper function (mirrors `EX.Country.part2()`, always creates an OBR experiment object) |
| `NAMESPACE` | Added `export(EX.Country.part2.OR)`, `exportClasses(experiment...outbreakresponse)`, `exportClasses(sim.results...outbreakresponse)` |
| `DESCRIPTION` | Added `'EX.Country.part2.OR.R'` to Collate after `'EX.Country.part2.R'` |
| `branch_building/OutbreakResponse.R` | Stripped to scenario construction example only (all other code now in package) |

Post-integration additions:
- Added `or.start.timestep` slot (class, run method, helper function) — lets user prevent OBR from triggering before a given time step
- `or.min.t` in run method is now `max(or.total.delay + or.trigger.window, or.start.timestep)`
- `or.vacc.coverage = 0` pattern established for "no OBR" baseline — both scenarios use `EX.Country.part2.OR()` for consistency
- Working Zambia example in `branch_building/OutbreakResponse.R` with two scenarios (A: `or.vacc.coverage = 0`, B: `or.vacc.coverage = 0.85`)

Notes on `EX.Country.part2.OR()`:
- `SIAinacc`, `SIAinefficient`, and `prop.inacc` are excluded — incompatible with the `vaccinationcorrelation` parent class
- `MR1MR2correlation`, `MR1SIAcorrelation`, `MR2SIAcorrelation` default to full correlation (`TRUE`, `1`, `1`) — appropriate for OBR scenarios where all campaigns are correlated
- `MR1SIAcorrelation` and `MR2SIAcorrelation` are wrapped in `as.numeric()` for backward compatibility with boolean input

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

### Current status and recommendation

No code has been written for the observational model. The discussion is design-level only.

**Recommended starting point**: Option 1. Use IgM-confirmed counts (suspected × proportion IgM+, aggregated annually) as the calibration target, and let `φ` absorb all sources of under-detection. Background rash as a nuisance parameter (Option 2) is the right longer-term framing, but Option 1 is simpler and avoids the need to estimate `B`.

The biggest contributors to background rash in Zambia are likely enteroviruses and parvovirus B19 (no vaccine, no model). Rubella noise is probably small given Zambia's vaccine introduction, and the Zambia-specific Chikungunya/dengue/Zika burden is low. **Background rash as a constant nuisance `B` is probably the right framing** rather than trying to model each pathogen separately.

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