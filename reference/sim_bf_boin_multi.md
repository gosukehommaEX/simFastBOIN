# Operating Characteristics of a BF-BOIN Design Across Several Scenarios

Run
[`sim_bf_boin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_bf_boin.md)
under a list of scenarios of DLT and response probabilities that share
the same design, and collect the results into one table.

## Usage

``` r
sim_bf_boin_multi(
  target,
  scenarios,
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

  Logical scalar. Keep the trial by trial data of every scenario.
  Defaults to `FALSE`.

- verbose:

  Logical scalar. Report progress while the simulations run. Defaults to
  `FALSE`.

- seed:

  Integer scalar or `NULL`. Random seed, applied to every scenario, so
  that a scenario simulated here matches the same scenario simulated on
  its own with
  [`sim_bf_boin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_bf_boin.md).
  Defaults to 123.

## Value

An object of class
`c("bf_boin_oc_multi", "backfill_oc_multi", "boin_oc_multi")`, which is
a list with components

- results:

  Named list of `bf_boin_oc` objects, one per scenario.

- summary_table:

  Data frame collecting all scenarios, with seven rows per scenario: the
  true DLT and response rates, the MTD selection, the patients treated,
  the backfilled patients, the patients with DLT and the average trial
  duration.

- scenario_names:

  Character vector of scenario names.

- n_doses:

  Number of dose levels.

- call:

  The matched call.

## See also

[`sim_bf_boin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_bf_boin.md),
[`sim_be_boin_multi`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_be_boin_multi.md)

## Examples

``` r
# Scenarios 2 and 6 of Zhao et al. (2024)
scenarios <- list(
  list(name = "Scenario 2", p_true = c(0.12, 0.25, 0.42, 0.49, 0.55),
       p_resp = c(0.20, 0.30, 0.40, 0.50, 0.60)),
  list(name = "Scenario 6", p_true = c(0.12, 0.25, 0.42, 0.49, 0.55),
       p_resp = c(0.30, 0.35, 0.36, 0.36, 0.36))
)

oc <- sim_bf_boin_multi(
  target = 0.25,
  scenarios = scenarios,
  n_cohort = 10,
  cohort_size = 3,
  window = 1,
  accrual_rate = 3,
  n_earlystop = 9,
  stay_on_1_of_3 = TRUE,
  n_trials = 200,
  seed = 123
)
oc
#> BF-BOIN operating characteristics across 2 scenarios
#>   target DLT rate : 0.25
#>   trials each     : 200
#>   escalation size : 30
#>   backfill cap    : 12 patients per dose
#>   window          : 1 (DLT), 1 (response)
#>   accrual         : 3 patients per unit of time (exponential)
#> 
#>    Scenario                   Item  DL1  DL2  DL3  DL4  DL5 Total / No MTD
#>  Scenario 2      True DLT rate (%) 12.0 25.0 42.0 49.0 55.0               
#>             True response rate (%) 20.0 30.0 40.0 50.0 60.0               
#>                   MTD selected (%) 30.0 59.5  9.0  0.5  0.0            1.0
#>                   Patients treated 10.8 10.9  3.9  0.6  0.0           26.2
#>                Backfilled patients  3.2  2.2  0.3  0.0  0.0            5.7
#>                  Patients with DLT  1.3  2.7  1.6  0.3  0.0            5.9
#>              Trial duration (mean)                                    11.9
#>  Scenario 6      True DLT rate (%) 12.0 25.0 42.0 49.0 55.0               
#>             True response rate (%) 30.0 35.0 36.0 36.0 36.0               
#>                   MTD selected (%) 29.5 59.0 11.0  0.0  0.0            0.5
#>                   Patients treated 11.6 10.8  3.7  0.4  0.0           26.5
#>                Backfilled patients  3.7  2.2  0.3  0.0  0.0            6.2
#>                  Patients with DLT  1.4  2.8  1.5  0.3  0.0            5.9
#>              Trial duration (mean)                                    11.9
```
