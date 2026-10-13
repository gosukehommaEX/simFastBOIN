# Operating Characteristics of a BF-BOIN Design

Simulate the BOIN design with backfilling (BF-BOIN) of Zhao et al.
(2024) many times under one scenario of DLT and response probabilities,
and summarize how often each dose is selected as the MTD, how many
patients are treated and backfilled at each dose, how many DLTs and
responses occur, and how long the trials last.

## Usage

``` r
sim_bf_boin(
  target,
  p_true,
  p_resp,
  n_cohort,
  cohort_size,
  window,
  accrual_rate,
  n_cap = 12,
  backfill_dose = c("highest", "lowest"),
  conflict_dose = c("highest", "lowest"),
  no_slot = c("wait", "leave"),
  accrual = c("exponential", "uniform", "fixed"),
  dlt_time = c("weibull", "uniform"),
  late_fraction = 0.5,
  resp_window = window,
  resp_late_fraction = 0.5,
  resp_cor = 0,
  n_trials = 10000,
  start_dose = 1,
  n_earlystop = 18,
  p_saf = NULL,
  p_tox = NULL,
  cutoff_eli = 0.95,
  extrasafe = FALSE,
  offset = 0.05,
  stay_on_1_of_3 = FALSE,
  bound_mtd = FALSE,
  mtd_max_estimate = NULL,
  min_mtd_sample = 1,
  overdose_cutoff = NULL,
  n_earlystop_rule = c("with_stay", "simple"),
  keep_trials = FALSE,
  verbose = FALSE,
  seed = 123
)
```

## Arguments

- target:

  Numeric scalar. Target DLT probability, for example 0.25.

- p_true:

  Numeric vector. True DLT probability at each dose level, in increasing
  dose order. The probability of a DLT within the assessment window.

- p_resp:

  Numeric vector of the same length as `p_true`. True probability of a
  response within the response window at each dose. A response at or
  below a dose is what makes the dose eligible for backfilling.

- n_cohort:

  Integer scalar. Number of cohorts in the dose escalation.

- cohort_size:

  Integer scalar or vector. Number of patients per cohort of the dose
  escalation. A scalar is used for every cohort. A vector shorter than
  `n_cohort` is padded with its last element and a longer one is
  truncated.

- window:

  Numeric scalar. Length of the DLT assessment window.

- accrual_rate:

  Numeric scalar. Average number of patients arriving per unit of time,
  in the same unit as `window`.

- n_cap:

  Integer scalar. A dose is no longer backfilled once this many
  patients, from the dose escalation and from backfilling together, have
  been treated at it. The dose escalation can still treat patients at
  the dose. Defaults to 12.

- backfill_dose:

  Character scalar. Which dose receives a backfill patient when several
  are open: `"highest"` (the default, as in both articles) or
  `"lowest"`.

- conflict_dose:

  Character scalar. When the data of several backfilled doses conflict
  with those of the current dose, which of them anchors the pooled
  estimate: `"highest"` (the default here, following Zhao et al., 2024)
  or `"lowest"` (the default of
  [`be_boin_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/be_boin_simulate.md),
  following Chen et al., 2026). See the details.

- no_slot:

  Character scalar. What happens to a patient who arrives when neither
  the dose escalation nor any backfill dose can take a patient: `"wait"`
  (the default), the patient waits and is treated as soon as one can, as
  in
  [`tite_boin_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/tite_boin_simulate.md);
  or `"leave"`, the patient is turned away and the trial waits for the
  next arrival.

- accrual:

  Character scalar. Distribution of the times between arrivals:
  `"exponential"` (the default) for a Poisson process, `"uniform"` for a
  uniform distribution on `(0, 2 / accrual_rate)`, or `"fixed"` for
  arrivals exactly `1 / accrual_rate` apart.

- dlt_time:

  Character scalar. Distribution of the time to DLT within the window,
  and of the time to response within the response window: `"weibull"`
  (the default) or `"uniform"`.

- late_fraction:

  Numeric scalar strictly between 0 and 1. Share of the DLTs that occur
  in the second half of the assessment window under the Weibull
  distribution. Defaults to 0.5.

- resp_window:

  Numeric scalar. Length of the window within which a response is
  observed. Defaults to `window`.

- resp_late_fraction:

  Numeric scalar strictly between 0 and 1. Share of the responses that
  occur in the second half of the response window under the Weibull
  distribution. Defaults to 0.5.

- resp_cor:

  Numeric scalar strictly between -1 and 1. Correlation of the Gaussian
  copula that links the DLT and the response of a patient; a positive
  value makes a response more likely in a patient with a DLT. Defaults
  to 0, independence.

- n_trials:

  Integer scalar. Number of trials to simulate. Defaults to 10000.

- start_dose:

  Integer scalar. Dose level for the first cohort. Defaults to 1.

- n_earlystop:

  Integer scalar. The trial stops once this many patients, backfilled
  ones included, have been treated at the current dose and the design
  would stay there. Defaults to 18.

- p_saf:

  Numeric scalar. Highest DLT probability deemed subtherapeutic.
  Defaults to `0.6 * target`.

- p_tox:

  Numeric scalar. Lowest DLT probability deemed overly toxic. Defaults
  to `1.4 * target`.

- cutoff_eli:

  Numeric scalar. Posterior probability cutoff for dose elimination.
  Defaults to 0.95.

- extrasafe:

  Logical scalar. Apply the stricter safety stopping rule at the lowest
  dose. Defaults to `FALSE`.

- offset:

  Numeric scalar between 0 and 0.5. Amount by which `cutoff_eli` is
  relaxed for the safety stopping rule. Defaults to 0.05.

- stay_on_1_of_3:

  Logical scalar. When `TRUE`, one DLT out of three patients leads to
  staying rather than de-escalating, in every use of the decision table.
  Zhao et al. (2024) use this modification in their main simulation.
  Defaults to `FALSE`. See
  [`boin_boundary`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_boundary.md).

- bound_mtd:

  Logical scalar. Require the isotonic estimate at the selected dose to
  be at or below the de-escalation boundary. Defaults to `FALSE`.

- mtd_max_estimate:

  Numeric scalar or `NULL`. Largest isotonic estimate a dose may have
  and still be selected as the MTD. See
  [`sim_boin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin.md).

- min_mtd_sample:

  Integer scalar. Smallest number of patients a dose must have received
  to be eligible as the MTD. Defaults to 1.

- overdose_cutoff:

  Numeric scalar or `NULL`. Doses whose true DLT probability exceeds
  this value count as overdoses in the `overdose` component of the
  result. Defaults to `NULL`, which uses `target`.

- n_earlystop_rule:

  Character scalar, either `"with_stay"` or `"simple"`. See
  [`boin_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_simulate.md).

- keep_trials:

  Logical scalar. Keep the full trial by trial data in the result.
  Defaults to `FALSE`.

- verbose:

  Logical scalar. Report progress while the simulation runs. Defaults to
  `FALSE`.

- seed:

  Integer scalar or `NULL`. Random seed. When `NULL` the current state
  of the random number generator is used. The state of the calling
  session is restored on exit. Defaults to 123.

## Value

An object of class `c("bf_boin_oc", "backfill_oc", "boin_oc")`, which
has every component of the result of
[`sim_boin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin.md),
so that code written for
[`sim_boin()`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin.md)
also reads it, counting backfilled patients among the patients treated,
and in addition

- n_bf_dose:

  Average number of backfilled patients at each dose.

- total_n_bf:

  Average number of backfilled patients per trial.

- pct_trials_backfill:

  Percentage of trials with at least one backfilled patient.

- n_resp_dose:

  Average number of patients with a response at each dose.

- total_n_resp:

  Average number of patients with a response per trial.

- duration_mean:

  Average trial duration, from the first arrival until every enrolled
  patient has completed the DLT assessment, in the time unit of `window`
  and `accrual_rate`.

- duration_sd:

  Standard deviation of the trial duration.

- pct_trials_suspended:

  Percentage of trials in which the dose escalation had to wait at least
  once.

- avg_n_suspensions:

  Average number of decisions per trial at which the dose escalation had
  to wait.

- avg_time_suspended:

  Average total time per trial during which the dose escalation waited.

- avg_n_turned_away:

  Average number of patients turned away per trial, which is zero unless
  `no_slot = "leave"`.

- p_resp:

  The true response probabilities.

The `trials` component, when kept, holds the trial by trial data of
[`bf_boin_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/bf_boin_simulate.md)
together with the selected MTD.

## Details

The trials are simulated by
[`bf_boin_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/bf_boin_simulate.md),
whose details describe the conduct of a trial. Once every patient has
completed the assessment the MTD is selected from the data of all
patients, backfilled ones included, exactly as in
[`sim_boin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin.md),
and a trial stopped for safety selects no dose. The overdose cutoffs of
60 and 80 percent refer to the size of the dose escalation, `n_cohort`
cohorts.

When no dose is ever opened for backfilling, for example with every
`p_resp` equal to 0, the selection and the patient counts agree with
those of
[`sim_boin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin.md)
under the same seed.

To reproduce Table 4 of Zhao et al. (2024), use `target = 0.25`,
`n_cohort = 10`, `cohort_size = 3`, `window = 1`, `accrual_rate = 3`,
`n_earlystop = 9`, `stay_on_1_of_3 = TRUE`, `accrual = "uniform"` and
`no_slot = "leave"`.

## References

Zhao, Y., Yuan, Y., Korn, E. L. and Freidlin, B. (2024). Backfilling
Patients in Phase I Dose-Escalation Trials Using Bayesian Optimal
Interval Design (BOIN). Clinical Cancer Research, 30(4), 673-679.

## See also

[`bf_boin_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/bf_boin_simulate.md),
[`sim_bf_boin_multi`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_bf_boin_multi.md),
[`sim_be_boin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_be_boin.md),
[`sim_boin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin.md)

## Examples

``` r
# Scenario 3 of Zhao et al. (2024)
oc <- sim_bf_boin(
  target = 0.25,
  p_true = c(0.04, 0.12, 0.25, 0.43, 0.63),
  p_resp = c(0.10, 0.20, 0.30, 0.45, 0.58),
  n_cohort = 10,
  cohort_size = 3,
  window = 1,
  accrual_rate = 3,
  n_earlystop = 9,
  stay_on_1_of_3 = TRUE,
  accrual = "uniform",
  no_slot = "leave",
  n_trials = 500,
  seed = 123
)
oc
#> BF-BOIN operating characteristics (Zhao et al., 2024)
#>   target DLT rate : 0.25
#>   trials          : 500
#>   cohorts         : 10 of size 3
#>   escalation size : 30
#>   backfill cap    : 12 patients per dose
#>   lambda_e / _d   : 0.1968 / 0.2984
#>   window          : 1 (DLT), 1 (response)
#>   accrual         : 3 patients per unit of time (uniform)
#> 
#>                         DL1  DL2  DL3  DL4  DL5 Total / No MTD
#> True DLT rate (%)       4.0 12.0 25.0 43.0 63.0               
#> True response rate (%) 10.0 20.0 30.0 45.0 58.0               
#> MTD selected (%)        2.6 30.0 57.8  9.6  0.0            0.0
#> Patients treated        6.0 10.8  9.8  3.7  0.4           30.7
#>   of which backfilled   1.8  3.3  1.9  0.2  0.0            7.2
#> Patients with DLT       0.3  1.3  2.4  1.6  0.3            5.9
#> 
#> Doses above a true DLT rate of 0.25: DL4, DL5
#>   Patients treated there              : 13.3%
#>   Trials treating any patient there   : 49.2%
#>   Trials selecting an MTD there       : 9.6%
#> 
#> Backfilling and timing
#>   Trials with backfilled patients     : 86.4%
#>   Average trial duration              : 14.3
#>   Trials in which the escalation waited: 100%
#>   Average time the escalation waited  : 6.2
#>   Average patients turned away        : 11.1
```
