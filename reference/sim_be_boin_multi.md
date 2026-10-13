# Operating Characteristics of a BE-BOIN Design Across Several Scenarios

Run
[`sim_be_boin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_be_boin.md)
under a list of scenarios of DLT and response probabilities that share
the same design, and collect the results into one table.

## Usage

``` r
sim_be_boin_multi(
  target,
  scenarios,
  n_cohort,
  cohort_size,
  window,
  accrual_rate,
  n_cap = 12,
  backfill_dose = c("highest", "lowest"),
  conflict_dose = c("lowest", "highest"),
  no_slot = c("wait", "leave"),
  accrual = c("exponential", "uniform", "fixed"),
  dlt_time = c("weibull", "uniform"),
  late_fraction = 0.5,
  resp_window = window,
  resp_late_fraction = 0.5,
  resp_cor = 0,
  max_pending_ratio = 0.49,
  min_completed = 0,
  min_follow_up = 0.25,
  n_trials = 10000,
  start_dose = 1,
  n_earlystop = 18,
  p_saf = NULL,
  p_tox = NULL,
  cutoff_eli = 0.95,
  extrasafe = FALSE,
  offset = 0.05,
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

- scenarios:

  A list of scenarios. Each element is a list with a `p_true` vector of
  true DLT probabilities, a `p_resp` vector of true response
  probabilities of the same length and optionally a `name`; without a
  name the list names are used, or "Scenario 1", "Scenario 2" and so on.
  Every scenario must have the same number of doses.

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
  estimate: `"lowest"` (the default here, following Chen et al., 2026)
  or `"highest"` (the default of
  [`bf_boin_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/bf_boin_simulate.md),
  following Zhao et al., 2024).

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

- max_pending_ratio:

  Numeric scalar greater than 0 and at most 1. Rule 1 of the design: the
  escalation is suspended when the proportion of pending patients at the
  current dose exceeds this value, unless the dose is de-escalated
  whatever their outcomes. Defaults to 0.49, which suspends when fewer
  than 51 percent have completed the assessment, as in Chen et al.
  (2026). See
  [`tite_boin_decision_table`](https://gosukehommaex.github.io/simFastBOIN/reference/tite_boin_decision_table.md).

- min_completed:

  Integer scalar. While some patients are pending, escalation requires
  at least this many patients with a completed assessment. Defaults to
  0, which the design does not use. See
  [`tite_boin_decision_table`](https://gosukehommaex.github.io/simFastBOIN/reference/tite_boin_decision_table.md).

- min_follow_up:

  Numeric scalar between 0 and 1. Rule 2 of the design: escalation
  requires every pending patient at the current dose to have been
  followed for at least this fraction of the window. Defaults to 0.25,
  as in Chen et al. (2026). See
  [`tite_boin_decision_table`](https://gosukehommaex.github.io/simFastBOIN/reference/tite_boin_decision_table.md).

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

  Integer scalar or `NULL`. Random seed, applied to every scenario, so
  that a scenario simulated here matches the same scenario simulated on
  its own with
  [`sim_be_boin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_be_boin.md).
  Defaults to 123.

## Value

An object of class
`c("be_boin_oc_multi", "backfill_oc_multi", "boin_oc_multi")`, with the
components described in
[`sim_bf_boin_multi`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_bf_boin_multi.md);
`results` holds `be_boin_oc` objects.

## See also

[`sim_be_boin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_be_boin.md),
[`sim_bf_boin_multi`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_bf_boin_multi.md)

## Examples

``` r
scenarios <- list(
  list(name = "MTD at dose 3", p_true = c(0.04, 0.12, 0.25, 0.43, 0.63),
       p_resp = c(0.10, 0.20, 0.30, 0.45, 0.58)),
  list(name = "Plateau", p_true = c(0.04, 0.12, 0.25, 0.43, 0.63),
       p_resp = c(0.30, 0.32, 0.35, 0.36, 0.36))
)

oc <- sim_be_boin_multi(
  target = 0.25,
  scenarios = scenarios,
  n_cohort = 10,
  cohort_size = 3,
  window = 3,
  accrual_rate = 2,
  n_earlystop = 9,
  n_trials = 200,
  seed = 123
)
oc
#> BE-BOIN operating characteristics across 2 scenarios
#>   target DLT rate : 0.25
#>   trials each     : 200
#>   escalation size : 30
#>   backfill cap    : 12 patients per dose
#>   window          : 3 (DLT), 3 (response)
#>   accrual         : 2 patients per unit of time (exponential)
#> 
#>       Scenario                   Item  DL1  DL2  DL3  DL4  DL5 Total / No MTD
#>  MTD at dose 3      True DLT rate (%)  4.0 12.0 25.0 43.0 63.0               
#>                True response rate (%) 10.0 20.0 30.0 45.0 58.0               
#>                      MTD selected (%)  4.0 33.5 52.5 10.0  0.0            0.0
#>                      Patients treated  7.4 11.0  8.5  3.0  0.3           30.1
#>                   Backfilled patients  2.1  2.8  1.2  0.1  0.0            6.2
#>                     Patients with DLT  0.3  1.3  2.2  1.3  0.2            5.3
#>                 Trial duration (mean)                                    22.7
#>        Plateau      True DLT rate (%)  4.0 12.0 25.0 43.0 63.0               
#>                True response rate (%) 30.0 32.0 35.0 36.0 36.0               
#>                      MTD selected (%)  5.0 27.5 55.0 12.5  0.0            0.0
#>                      Patients treated  9.5 11.2  8.8  3.3  0.4           33.2
#>                   Backfilled patients  4.5  3.5  1.4  0.1  0.0            9.6
#>                     Patients with DLT  0.3  1.3  2.2  1.3  0.3            5.4
#>                 Trial duration (mean)                                    22.8
```
