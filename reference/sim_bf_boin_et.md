# Operating Characteristics of the BF-BOIN-ET Design

Simulate trials of the backfill BOIN-ET (BF-BOIN-ET) design of Takeda,
Zhu and Hirakawa (2025) under one scenario of DLT and response
probabilities, select the optimal biological dose (OBD) of each trial,
and summarize how often each dose is selected, how many patients are
treated and backfilled at each dose, and how long the trials last.

## Usage

``` r
sim_bf_boin_et(
  target,
  target_resp,
  p_true,
  p_resp,
  n_cohort,
  cohort_size,
  window,
  accrual_rate,
  resp_window = window,
  accrual = c("uniform", "exponential", "fixed"),
  event_time = c("weibull", "uniform"),
  late_fraction = 0.5,
  resp_late_fraction = 0.5,
  resp_cor = 0.2,
  n_trials = 10000,
  start_dose = 1,
  n_earlystop = 12,
  n_elimination = 6,
  suspend = 0,
  p_saf = NULL,
  p_tox = NULL,
  p_saf_resp = NULL,
  cutoff_eli = 0.95,
  cutoff_eli_resp = 0.9,
  psi00 = 40,
  psi11 = 60,
  psi_resp_no_dlt = 100,
  psi_dlt_no_resp = 0,
  backfill_extension = 100,
  keep_trials = FALSE,
  verbose = FALSE,
  seed = 123
)
```

## Arguments

- target:

  Numeric scalar. Target DLT probability.

- target_resp:

  Numeric scalar. Lowest response probability of interest. A dose is
  excluded when its response rate is likely below it, and the backfill
  doses start from the lowest dose whose response rate reaches the
  boundary `eta1` derived from it.

- p_true:

  Numeric vector. True DLT probability at each dose level, in increasing
  dose order.

- p_resp:

  Numeric vector of the same length as `p_true`. True response
  probability at each dose level.

- n_cohort:

  Integer scalar. Number of cohorts in a trial.

- cohort_size:

  Integer scalar or vector. Number of patients per cohort. A scalar is
  used for every cohort. A vector shorter than `n_cohort` is padded with
  its last element and a longer one is truncated.

- window:

  Numeric scalar. Length of the DLT assessment window.

- accrual_rate:

  Numeric scalar. Average number of patients arriving per unit of time,
  in the same unit as `window`. The `accrual` argument of the 'boinet'
  package is the average time between arrivals, its inverse.

- resp_window:

  Numeric scalar. Length of the response assessment window. Defaults to
  `window`.

- accrual:

  Character scalar. Distribution of the times between arrivals:
  `"uniform"` (the default, as in 'boinet') on `(0, 2 / accrual_rate)`,
  `"exponential"` for a Poisson process, or `"fixed"` for arrivals
  exactly `1 / accrual_rate` apart.

- event_time:

  Character scalar. Distribution of the times to DLT and to response:
  `"weibull"` (the default), with the share of the events in the second
  half of the window set by `late_fraction` and `resp_late_fraction`, or
  `"uniform"` on `(0, window / p)` for an event of probability `p`. With
  `"weibull"` the probabilities must be below 1.

- late_fraction:

  Numeric scalar strictly between 0 and 1. Share of the DLTs that occur
  in the second half of the DLT window under the Weibull distribution.
  Defaults to 0.5.

- resp_late_fraction:

  Numeric scalar strictly between 0 and 1. The same for the responses.
  Defaults to 0.5.

- resp_cor:

  Numeric scalar strictly between -1 and 1. Correlation of the Gaussian
  copula that links the time to DLT and the time to response of a
  patient. Defaults to 0.2, the default of the 'boinet' package.

- n_trials:

  Integer scalar. Number of trials to simulate. Defaults to 10000.

- start_dose:

  Integer scalar. Dose level for the first cohort. Defaults to 1.

- n_earlystop:

  Integer scalar. Number of patients at a dose, dose escalation and
  backfill together, that closes it to backfilling, and that ends the
  trial when it is reached at the current dose and the decision would
  keep the trial there. Defaults to 12, the default of the 'bfboinet'
  package.

- n_elimination:

  Integer scalar. Number of patients a dose needs before it can be
  excluded. Defaults to 6.

- suspend:

  Numeric scalar from 0 to 1. Share of the patients at the current dose
  whose DLT assessment may still be running when the next cohort starts.
  Defaults to 0, which waits until the DLT outcomes of the cohort are
  known.

- p_saf:

  Numeric scalar or `NULL`. Highest DLT probability deemed
  subtherapeutic, used for the escalation boundary. Defaults to
  `0.6 * target`.

- p_tox:

  Numeric scalar or `NULL`. Lowest DLT probability deemed overly toxic,
  used for the de-escalation boundary. Defaults to `1.4 * target`.

- p_saf_resp:

  Numeric scalar or `NULL`. Response probability deemed subtherapeutic,
  used for the response boundary. Defaults to `0.6 * target_resp`.

- cutoff_eli:

  Numeric scalar. A dose with `n_elimination` patients or more is
  excluded, together with every higher dose, when the posterior
  probability that its DLT rate exceeds `target` is above this cutoff.
  Defaults to 0.95.

- cutoff_eli_resp:

  Numeric scalar. Such a dose is excluded alone when the posterior
  probability that its response rate is below `target_resp` is above
  this cutoff. Defaults to 0.90.

- psi00, psi11:

  Numeric scalars. Scores of a patient without DLT or response and of a
  patient with both, which define, with `psi_resp_no_dlt` and
  `psi_dlt_no_resp`, the utility used to choose the backfill dose and
  the OBD. Default to 40 and 60.

- psi_resp_no_dlt, psi_dlt_no_resp:

  Numeric scalars. Scores of a patient with a response and no DLT and of
  a patient with a DLT and no response. Default to 100 and 0. The
  'bfboinet' package calls them `ppsi10` and `ppsi01`.

- backfill_extension:

  Numeric scalar. Time, in the unit of `window`, during which
  backfilling continues after the last cohort has been enrolled.
  Defaults to 100, the value fixed in the 'bfboinet' package.

- keep_trials:

  Logical scalar. Keep the simulated trials and the selected OBDs in the
  result. Defaults to `FALSE`.

- verbose:

  Logical scalar. Report progress with messages. Defaults to `FALSE`.

- seed:

  Integer scalar or `NULL`. Random seed. When `NULL` the current state
  of the random number generator is used. The state of the calling
  session is restored on exit. Defaults to 123.

## Value

An object of class `bf_boin_et_oc`, which is a list with components

- design:

  `"BF-BOIN-ET"`.

- target, target_resp, p_true, p_resp:

  The targets and the scenario.

- n_doses, n_trials:

  Numbers of doses and of simulated trials.

- lambda1, lambda2, eta1:

  The escalation and de-escalation boundaries of the DLT rate and the
  boundary of the response rate.

- sel_percent:

  Percentage of trials selecting each dose as the OBD.

- percent_no_obd:

  Percentage of trials selecting no OBD.

- percent_early_stop:

  Percentage of trials stopped because every dose was excluded.

- n_pts_dose, n_bf_dose:

  Average numbers of patients treated and of patients backfilled at each
  dose.

- n_tox_dose, n_resp_dose:

  Average numbers of DLTs and responses at each dose observed at the
  last decision.

- n_tox_bf_dose, n_resp_bf_dose:

  The same among the backfilled patients.

- total_n_pts, total_n_bf, total_n_tox, total_n_resp:

  Their averages per trial.

- duration:

  Average time until the last patient finishes both assessments.

- stop_reason_percent:

  Percentage of trials ending for each reason, see Details.

- settings:

  List of the design parameters.

- call:

  The matched call.

- trials:

  When `keep_trials` is `TRUE`, a list with the matrices `n_pts`,
  `n_bf`, `n_tox`, `n_resp`, `n_tox_bf` and `n_resp_bf` (one row per
  trial), and the vectors `cohorts_used`, `stop_reason`, `early_stop`,
  `obd` (`NA` when none) and `duration`.

## Details

The rules are those of
[`bfboinet::get.oc.backboinet()`](https://rdrr.io/pkg/bfboinet/man/get.oc.backboinet.html)
(version 1.1.0). A cohort is treated at the current dose. The next
cohort starts with the first patient who arrives once the DLT outcomes
of the cohort are known (or, with `suspend`, once enough of them are),
and after the last cohort backfilling goes on for `backfill_extension`.
The patients who arrive in between are backfilled one by one, when the
current dose is not the lowest, to the backfill dose with the highest
utility (ties to the higher dose), the DLT rates being made monotone by
isotonic regression; doses with `n_earlystop` patients take no more.

The decision uses the DLT rate at the current dose, which counts a
patient still under assessment by the fraction of the window followed,
with the boundaries `lambda1` and `lambda2` of BOIN; after 1 DLT in 3
patients the trial stays. The backfill doses for the next cohort are the
admissible doses from the lowest dose whose response rate is at least
`eta1` up to the current dose (escalation), the dose below it (stay) or
two doses below it (de-escalation). Once excluded, a dose stays
excluded, and an escalation moves past a dose excluded by its response
data only. The OBD is the admissible dose up to the estimated MTD with
the highest utility, from the observed rates at the last decision; ties
go to the lower dose. The utilities are compared after rounding to ten
decimals, for the backfill dose and for the OBD, so that doses whose
utilities are equal as fractions tie whatever the rounding error of the
platform; 'bfboinet' and versions of this package up to 2.2.0 compare
them as computed.

Possible values of `stop_reason` are `"no_admissible_dose"` (the only
early stop, which selects no OBD when fewer patients than the maximum
sample size have been treated), `"lowest_dose_excluded"`,
`"top_dose_excluded"`, `"no_lower_admissible_dose"`, `"n_earlystop"` and
`"max_cohorts"`.

The 'bfboinet' package reports the position of the OBD among the doses
tried, which is the dose itself when the trial starts at the lowest
dose; the dose is reported here in every case.

## References

Takeda, K., Zhu, J. and Hirakawa, A. (2025). BF-BOIN-ET: A backfill
Bayesian optimal interval design using efficacy and toxicity outcomes
for dose optimization. Pharmaceutical Statistics, 24, e2470.
[doi:10.1002/pst.2470](https://doi.org/10.1002/pst.2470)

## See also

[`sim_bf_boin_et_multi`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_bf_boin_et_multi.md),
[`sim_tite_boin_et`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_tite_boin_et.md)

## Examples

``` r
# The example of bfboinet::get.oc.backboinet(): windows of one month and
# three patients a month
oc <- sim_bf_boin_et(
  target = 0.30, target_resp = 0.25,
  p_true = c(0.03, 0.05, 0.20, 0.22, 0.45),
  p_resp = c(0.05, 0.10, 0.50, 0.68, 0.70),
  n_cohort = 10, cohort_size = 3, window = 1, accrual_rate = 3,
  n_trials = 300, seed = 123
)
oc
#> BF-BOIN-ET operating characteristics
#>   target DLT rate      : 0.3
#>   target response rate : 0.25
#>   trials               : 300
#>   cohorts              : 10 of size 3
#>   patients per dose    : 12 at most with backfill
#>   lambda1 / lambda2    : 0.2365 / 0.3585
#>   eta1                 : 0.1968
#> 
#>                        DL1  DL2  DL3  DL4  DL5 Total / No OBD
#> True DLT rate (%)      3.0  5.0 20.0 22.0 45.0               
#> True response rate (%) 5.0 10.0 50.0 68.0 70.0               
#> OBD selected (%)       3.0  2.7 41.3 48.3  4.7            0.0
#> Patients treated       4.1  6.0 11.4  9.5  5.0           36.1
#> Patients backfilled    0.7  1.5  4.5  2.5  0.1            9.4
#> Patients with DLT      0.1  0.3  2.3  2.0  2.4            7.1
#> Patients with response 0.2  0.6  5.7  6.5  3.6           16.6
#> 
#> Trials stopped early: 0%
#> Average trial duration: 16.6
```
