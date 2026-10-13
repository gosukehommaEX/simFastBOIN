# Simulate BE-BOIN Trials

Run the BOIN design with backfilling and late-onset toxicity (BE-BOIN)
of Chen et al. (2026) over many simulated trials. Dose escalation
follows TITE-BOIN, which decides while some patients are still pending;
patients who arrive while the escalation is suspended are backfilled to
a lower dose that is safe and has shown a response. Every estimate
imputes the pending patients, backfilled ones included. Return the
patient, DLT, backfill and response counts at every dose together with
the duration of each trial.

## Usage

``` r
be_boin_simulate(
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
  n_earlystop_rule = c("with_stay", "simple"),
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

- n_earlystop_rule:

  Character scalar, either `"with_stay"` or `"simple"`. See
  [`boin_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_simulate.md).

- seed:

  Integer scalar or `NULL`. Random seed. When `NULL` the current state
  of the random number generator is used. The state of the calling
  session is restored on exit. Defaults to 123.

## Value

An object of class `c("be_boin_trials", "backfill_trials")`, with the
components described in
[`bf_boin_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/bf_boin_simulate.md).

## Details

**Dose escalation.** The dose for a cohort is decided when its first
patient arrives, with the rules of
[`tite_boin_decision_table`](https://gosukehommaex.github.io/simFastBOIN/reference/tite_boin_decision_table.md)
for `method = "imputation"` and the suspension rules `max_pending_ratio`
and `min_follow_up`; every patient at the current dose counts,
backfilled ones included. A suspension lasts until a decision can be
taken again, as in
[`tite_boin_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/tite_boin_simulate.md).
Elimination counts every treated patient, the pending ones as without
DLT, and is checked at every dose.

**Estimates.** The DLT rate at a dose is the single mean imputation
estimate of TITE-BOIN, equation (1) of Chen et al. (2026), and the rate
pooled over several doses adds up the observed and imputed DLTs of each
dose, each imputed with its own posterior mean, and divides by the
number of treated patients, equation (2). Both enter the opening and
closing of backfill doses and the resolution of conflicts, which
otherwise follow
[`bf_boin_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/bf_boin_simulate.md).
The categories of a backfilled dose for the conflicts compare its
estimate with the escalation and de-escalation boundaries directly.

**End of the trial, arrivals and random numbers** are as in
[`bf_boin_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/bf_boin_simulate.md).
When no dose is ever opened for backfilling, for example with every
`p_resp` equal to 0, and with `no_slot = "wait"`, the trials are
identical to those of
[`tite_boin_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/tite_boin_simulate.md)
with `method = "imputation"`, the same suspension rules and the same
seed.

Chen et al. (2026) report their operating characteristics for scenarios
that are given only as figures, so they cannot be reproduced. The
decision table of their supplementary materials is reproduced by
[`tite_boin_decision_table`](https://gosukehommaex.github.io/simFastBOIN/reference/tite_boin_decision_table.md),
and their trial example is replayed patient by patient in the tests of
the package.

## References

Chen, K., Zhao, Y., Takeda, K. and Yuan, Y. (2026). BE-BOIN: A Dose
Optimization Design Accommodating Backfill and Late-Onset Toxicity.
Therapeutic Innovation and Regulatory Science.
[doi:10.1007/s43441-026-00994-0](https://doi.org/10.1007/s43441-026-00994-0)

Chen, K., Chen, T.-Y., Zhang, Y., Lin, R. and Yuan, Y. (2025). Practical
Considerations for Using the TITE-BOIN Design to Handle Late-Onset
Toxicity or Fast Accrual in Phase I Trials. Clinical Cancer Research,
31(13), 2573-2580.

Zhao, Y., Yuan, Y., Korn, E. L. and Freidlin, B. (2024). Backfilling
Patients in Phase I Dose-Escalation Trials Using Bayesian Optimal
Interval Design (BOIN). Clinical Cancer Research, 30(4), 673-679.

## See also

[`bf_boin_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/bf_boin_simulate.md),
[`tite_boin_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/tite_boin_simulate.md)

## Examples

``` r
# A three month window with two patients a month
trials <- be_boin_simulate(
  target = 0.25,
  p_true = c(0.04, 0.12, 0.25, 0.43, 0.63),
  p_resp = c(0.10, 0.20, 0.30, 0.45, 0.58),
  n_cohort = 10,
  cohort_size = 3,
  window = 3,
  accrual_rate = 2,
  n_earlystop = 9,
  n_trials = 200,
  seed = 123
)
trials
#> Simulated BE-BOIN (Chen et al., 2026) trials
#>   trials          : 200
#>   doses           : 5
#>   target DLT rate : 0.25
#>   cohorts         : 10 of size 3
#>   escalation size : 30
#>   backfill cap    : 12 patients per dose
#>   window          : 3 (DLT), 3 (response)
#>   accrual         : 2 patients per unit of time (exponential)
#> 
#> Average per trial:
#>                        DL1   DL2  DL3  DL4  DL5
#> Patients treated      7.38 10.96 8.47 2.98 0.30
#>   of which backfilled 2.08  2.82 1.20 0.07 0.00
#> DLTs observed         0.30  1.31 2.16 1.26 0.22
#> Responses             0.66  2.22 2.58 1.30 0.20
#> 
#> Patients per trial    : mean 30.09, of which backfilled 6.16
#> Trials with backfill  : 71.5%
#> Trial duration        : mean 22.72, median 23.19
#> Escalation waiting    : 100% of trials, mean time 7.93 per trial
#> 
#> Stopping reason (%):
#> 
#> max_sample_size     n_earlystop 
#>              35              65 

summary(trials$duration)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>   6.904  19.473  23.188  22.716  26.263  33.632 
```
