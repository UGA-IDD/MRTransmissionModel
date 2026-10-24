# Measles Surveillance Observational Model
## From Simulated True Infections to Reported and Confirmed Cases

---

## Overview

A measles transmission model produces **true infections**. But a surveillance system does not observe true infections — it observes **suspected reported cases**, some of which are confirmed by diagnostic testing. To evaluate surveillance programs, we need an observational model that bridges simulated truth and observed surveillance output.

The pipeline has two conceptually distinct parts:

1. **True infections → Suspected reported cases** (the reporting process)
2. **Suspected reported cases → Confirmed vs. Discarded** (the diagnostic process)

---

## Part 1: True Infections → Suspected Reported Cases

### The Key Distinction

Suspected cases are **not** simply a fraction of true infections. They come from two completely separate processes:

```
suspected_reported(t) = true_reported_infections(t) + non-measles reported cases(t)
```

- `true_reported_infections(t)` — real measles cases that entered the surveillance system
- `non-measles reported cases(t)` — non-measles fever-rash illness (rubella, roseola, dengue, etc.) that clinicians could not rule out and reported as suspected measles


The non-measles reported cases are background surveillance noise — they exist independently of whether measles is circulating.

---

### Step 1: Estimating True Reported Infections

```
true_reported_infections(t) = true_infections(t) × r
```

Where true_infections(t) comes directly from the simulations, and **r** is the overall reporting rate (i.e., the probability that a true infection becomes a true reported case). See below for more information about estimating **r**.

---

### Step 2: Estimating non-measles reported cases

Fnon-measles reported cases are derived directly from empirical surveillance data using the month-specific IgM positivity rate:

```
non-measles reported cases(m) = number_suspected_reported(m) × [1 - positivity_rate(m)]
```

Where `m` indexes calendar month (1–12), giving a seasonal profile of the false suspect background.

The IgM positivity rate varies by month — it is **low** when non-measles fever-rash illness dominates (many false suspects, few true cases) and **high** when measles is circulating strongly. This seasonal variation reflects the **ratio of true to false suspects** in the surveillance system, not just absolute measles burden.



### Why Empirical Suspected Cases and Simulated True Reported Infections Are on the Same Scale

A natural concern: are empirically-derived non-measles reported cases counts on the same scale as the true reported infections from the simulated data?

**The answer is yes, and here is why:**

- Empirical suspected cases are **reported suspected cases** — already filtered through the surveillance system
- Simulated true reported infections are `true_infections × r` — also already filtered through the same surveillance system by construction
- The 5-rule calibration (see below) forced **r** to make those two things live on the same observational scale

This means the false suspect background derived from empirical data and the true reported infections derived from simulation are directly compatible in magnitude. No further rescaling is needed.

---

## Part 2: Suspected Reported Cases → Confirmed vs. Discarded

### The Role of Se and Sp

The IgM diagnostic test has known sensitivity (Se) and specificity (Sp). These are applied separately to the two components of suspected cases.

**Apply Se to true reported infections:**

```
TP(t) = true_reported_infections(t) × Se
FN(t) = true_reported_infections(t) × (1 - Se)
```

- **TP**: real measles cases correctly confirmed by the test
- **FN**: real measles cases missed by the test and discarded

**Apply Sp to non-measles reported cases:**

```
TN(t) = non-measles reported cases(t) × Sp
FP(t) = non-measles reported cases(t) × (1 - Sp)
```

- **TN**: non-measles cases correctly ruled out by the test
- **FP**: non-measles cases incorrectly confirmed as measles

---

### What the Surveillance System Observes

```
confirmed(t) = TP(t) + FP(t)
discarded(t) = TN(t) + FN(t)
```

The full 2×2 table {TP, FP, TN, FN} constitutes a complete observational model sitting on top of the transmission model. These quantities can be sampled from stochastically in simulation.

---

## Calibrating the Reporting Rate r

### The 5-Rule

The reporting rate **r** (`or.reporting.rate` in the model) is calibrated using a program decision rule: OBR campaigns are triggered when IgM-confirmed cases (obs.TP) accumulate to **5** over the surveillance window. This threshold is operationally meaningful — it represents the point at which the surveillance system signals a genuine outbreak requiring response.

Note on model architecture: `or.reporting.rate` is the *pre-testing* fraction (care-seeking × clinical recognition). The obs model then applies Se on top, and the trigger sums obs.TP over `or.trigger.window` steps. Two calibration formulas are therefore possible:

```
Simple (r absorbs everything):
  r = 5 / number_infections_at_outbreak

Model-consistent (Se and window applied separately in obs model):
  r = 5 / (number_infections_at_outbreak × Se × or.trigger.window)
```

`number_infections_at_outbreak` is the mean number of true infections at epidemic onset (Rt first crossing above 1), estimated from stochastic simulation runs.

---

### Estimating number_infections_at_outbreak via Rt > 1

The key quantity is estimated **purely from transmission model output** — no observed data required.

Rt > 1 crossings are detected operationally in the stochastic I time series: since every infected individual recovers in exactly one time step, I(t)/I(t−1) ≈ Rt. A sustained crossing is defined as `n.consec = 3` consecutive steps (6 weeks) of growing I, following a trough. After each onset is flagged, epidemic state resets once I drops below `reset.fraction × I_at_onset` to avoid re-counting the same epidemic.

**Procedure:**

1. Run `n.reps` stochastic replicates using `EX.Country.part2()` (no OBR machinery needed — just raw transmission + vaccination dynamics).
2. For each replicate, restrict to the **evaluation window** (2015–2025) and find all epidemic onset timesteps.
3. Record total true I at each onset; pool across all runs and replicates.
4. Average:

```
number_infections_at_outbreak = mean(total I at all detected onsets)
```

5. Apply the chosen formula to get `r`.

---

### Vaccination coverage strategy for simulation

The choice of coverage inputs matters because it determines the susceptibility structure at epidemic onset — and therefore the value of `number_infections_at_outbreak`.

- **Too low coverage (e.g., vaccine-free)**: I at onset is artificially high; r would be underestimated for realistic Zambia conditions.
- **Actual current coverage throughout**: few outbreaks occur in the eval window (Zambia approaches herd immunity), yielding too few onset events for reliable estimation.

**Implemented approach:**

- **1980–2014**: actual historical Zambia MCV1, MCV2, and SIA coverage estimates. This builds the correct 2015 susceptibility profile — the immunity earned by the real vaccination programme.
- **2015–2024**: all coverage zeroed. This removes vaccination as a competing force during the evaluation window, ensuring outbreaks occur and providing a sufficient sample of onset events across stochastic runs.

The resulting onsets occur in a population with a realistic 2015 immune profile but without ongoing vaccination — a deliberate design choice to get enough epidemic events for calibration. The I_at_onset values therefore reflect what an outbreak looks like given the historical immunity background.

---

### Allowing r to Vary by Year

The reporting rate may not be stable over time — surveillance capacity, health system investment, and reporting infrastructure can improve or degrade across years. To capture this:

1. **Estimate r(year)** separately for each year using the Rt > 1 anchor applied to simulation runs calibrated to that year's conditions, solving:

```
r(year) = 5 / number_infections_at_outbreak(year)
```

2. **Fit a GLM** with `r(year)` as the outcome and year:

```
g(r(year)) = β0 + β1 × year + ε
```

Where `g(·)` is an appropriate link function (e.g., logit, since r is bounded [0,1]).

3. **Project r into future years** using the fitted GLM for forward simulation under different surveillance scenarios.

This allows the observational model to reflect realistic temporal trends in surveillance quality, and lets intervention scenarios be evaluated against a plausible future reporting rate rather than a fixed historical one.

---

## Parameter Summary

| Parameter | Source | Notes |
|---|---|---|
| `r` | Calibrated via 5-rule + Rt > 1 | Can vary by year; projected forward via GLM |
| `positivity_rate(m)` | Empirical IgM data pooled by calendar month | Captures seasonal co-circulation of fever-rash illness |
| `suspected_reported(m)` | Empirical surveillance data by month | Already on surveillance scale; compatible with simulation by construction |
| `Se` | Diagnostic literature or local validation | Applied to true reported infections |
| `Sp` | Diagnostic literature or local validation | Applied to false reported suspects |

---


