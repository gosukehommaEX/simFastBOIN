# Operating Characteristics of a BE-BOIN Design

Simulate the BOIN design with backfilling and late-onset toxicity
(BE-BOIN) of Chen et al. (2026) many times under one scenario of DLT and
response probabilities, and summarize how often each dose is selected as
the MTD, how many patients are treated and backfilled at each dose, how
many DLTs and responses occur, and how long the trials last.

## Usage

``` r
sim_be_boin(
  target,
  p_true,
  p_resp,
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

  Integer scalar or `NULL`. Random seed. When `NULL` the current state
  of the random number generator is used. The state of the calling
  session is restored on exit. Defaults to 123.

## Value

An object of class `c("be_boin_oc", "backfill_oc", "boin_oc")`, with the
components described in
[`sim_bf_boin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_bf_boin.md).

## Details

The trials are simulated by
[`be_boin_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/be_boin_simulate.md),
whose details describe the conduct of a trial. The MTD is selected as in
[`sim_bf_boin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_bf_boin.md).

When no dose is ever opened for backfilling, for example with every
`p_resp` equal to 0, and with `no_slot = "wait"`, every result that the
two functions share agrees with that of
[`sim_tite_boin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_tite_boin.md)
with `method = "imputation"`, the same suspension rules and the same
seed.

Only the MTD is selected. The second stage of Chen et al. (2026), which
randomizes patients between doses to choose the optimal biological dose,
is not part of the simulation.

## References

Chen, K., Zhao, Y., Takeda, K. and Yuan, Y. (2026). BE-BOIN: A Dose
Optimization Design Accommodating Backfill and Late-Onset Toxicity.
Therapeutic Innovation and Regulatory Science.
[doi:10.1007/s43441-026-00994-0](https://doi.org/10.1007/s43441-026-00994-0)

## See also

[`be_boin_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/be_boin_simulate.md),
[`sim_be_boin_multi`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_be_boin_multi.md),
[`sim_bf_boin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_bf_boin.md),
[`sim_tite_boin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_tite_boin.md)

## Examples

``` r
# A three month window with two patients a month
oc <- sim_be_boin(
  target = 0.25,
  p_true = c(0.04, 0.12, 0.25, 0.43, 0.63),
  p_resp = c(0.10, 0.20, 0.30, 0.45, 0.58),
  n_cohort = 10,
  cohort_size = 3,
  window = 3,
  accrual_rate = 2,
  n_earlystop = 9,
  n_trials = 500,
  seed = 123
)
oc
#> BE-BOIN operating characteristics (Chen et al., 2026)
#>   target DLT rate : 0.25
#>   trials          : 500
#>   cohorts         : 10 of size 3
#>   escalation size : 30
#>   backfill cap    : 12 patients per dose
#>   lambda_e / _d   : 0.1968 / 0.2984
#>   window          : 3 (DLT), 3 (response)
#>   accrual         : 2 patients per unit of time (exponential)
#> 
#>                         DL1  DL2  DL3  DL4  DL5 Total / No MTD
#> True DLT rate (%)       4.0 12.0 25.0 43.0 63.0               
#> True response rate (%) 10.0 20.0 30.0 45.0 58.0               
#> MTD selected (%)        3.6 32.4 52.2 11.6  0.2            0.0
#> Patients treated        7.2 10.7  8.4  3.2  0.4           29.9
#>   of which backfilled   2.0  2.8  1.2  0.1  0.0            6.1
#> Patients with DLT       0.3  1.4  2.1  1.3  0.3            5.3
#> 
#> Doses above a true DLT rate of 0.25: DL4, DL5
#>   Patients treated there              : 11.9%
#>   Trials treating any patient there   : 46.6%
#>   Trials selecting an MTD there       : 11.8%
#> 
#> Backfilling and timing
#>   Trials with backfilled patients     : 69.8%
#>   Average trial duration              : 22.4
#>   Trials in which the escalation waited: 100%
#>   Average time the escalation waited  : 8.1
```
