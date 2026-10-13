# Operating Characteristics of a TITE-BOIN Design

Simulate a time-to-event BOIN (TITE-BOIN) trial many times under one
dose-toxicity scenario and summarize how often each dose is selected as
the MTD, how many patients are treated and how many DLTs are observed at
each dose, how long the trials last and how often accrual is suspended.

## Usage

``` r
sim_tite_boin(
  target,
  p_true,
  n_cohort,
  cohort_size,
  window,
  accrual_rate,
  method = c("imputation", "ess"),
  accrual = c("exponential", "uniform", "fixed"),
  dlt_time = c("weibull", "uniform"),
  late_fraction = 0.5,
  prior_weights = c(1, 1, 1)/3,
  max_pending_ratio = NULL,
  min_completed = NULL,
  min_follow_up = 0,
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

  Numeric scalar. Target DLT probability, for example 0.30.

- p_true:

  Numeric vector. True DLT probability at each dose level, in increasing
  dose order. The probability of a DLT within the assessment window.

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
  in the same unit as `window`.

- method:

  Character scalar, `"imputation"` (the default) or `"ess"`. How the
  pending patients enter the estimate of the DLT probability. See
  [`tite_boin_decision_table`](https://gosukehommaex.github.io/simFastBOIN/reference/tite_boin_decision_table.md).

- accrual:

  Character scalar. Distribution of the times between arrivals:
  `"exponential"` (the default) for a Poisson process, `"uniform"` for a
  uniform distribution on `(0, 2 / accrual_rate)`, as in the code of Lin
  and Yuan (2020), or `"fixed"` for arrivals exactly `1 / accrual_rate`
  apart.

- dlt_time:

  Character scalar. Distribution of the time to DLT within the window:
  `"weibull"` (the default), with the share of DLTs in the second half
  of the window set by `late_fraction`, or `"uniform"`.

- late_fraction:

  Numeric scalar strictly between 0 and 1. Share of the DLTs that occur
  in the second half of the assessment window under the Weibull
  distribution. Defaults to 0.5, the setting of Yuan et al. (2018) and
  Lin and Yuan (2020).

- prior_weights:

  Numeric vector of three non-negative values. Prior probabilities that
  a DLT occurs in the first, second and last third of the assessment
  window, used to weight the follow-up of the pending patients in STFT
  (Yuan et al., 2018, Supplementary Appendix D) and ESS (Lin and Yuan,
  2020, Supplementary S1). They are rescaled to sum to one. Defaults to
  equal weights, which give the unweighted STFT. The decision table does
  not change with the weights; only the value compared with its
  boundaries does. This describes the design and is separate from
  `dlt_time` and `late_fraction`, which describe how the simulated DLTs
  actually occur.

- max_pending_ratio:

  Numeric scalar or `NULL`. Suspension rule on the proportion of pending
  patients. See
  [`tite_boin_decision_table`](https://gosukehommaex.github.io/simFastBOIN/reference/tite_boin_decision_table.md).

- min_completed:

  Integer scalar or `NULL`. Number of completed assessments required for
  escalation. See
  [`tite_boin_decision_table`](https://gosukehommaex.github.io/simFastBOIN/reference/tite_boin_decision_table.md).

- min_follow_up:

  Numeric scalar between 0 and 1. Shortest follow-up, as a fraction of
  `window`, that every pending patient at the current dose must have
  reached before the dose is escalated; accrual is suspended until then.
  Defaults to 0, which never suspends on this ground. See
  [`tite_boin_decision_table`](https://gosukehommaex.github.io/simFastBOIN/reference/tite_boin_decision_table.md).

- n_trials:

  Integer scalar. Number of trials to simulate. Defaults to 10000.

- start_dose:

  Integer scalar. Dose level for the first cohort. Defaults to 1.

- n_earlystop:

  Integer scalar. The trial stops once this many patients have been
  treated at the current dose and the design would stay there. Defaults
  to 18.

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

An object of class `c("tite_boin_oc", "boin_oc")`, which has every
component of the result of
[`sim_boin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin.md),
so that code written for
[`sim_boin()`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin.md)
also reads it, and in addition

- duration_mean:

  Average trial duration, from the first arrival until every enrolled
  patient has completed the assessment, in the time unit of `window` and
  `accrual_rate`.

- duration_sd:

  Standard deviation of the trial duration.

- pct_trials_suspended:

  Percentage of trials in which accrual was suspended at least once.

- avg_n_suspensions:

  Average number of suspensions per trial.

- avg_time_suspended:

  Average total time with accrual suspended per trial.

The `trials` component, when kept, also holds `duration`,
`n_suspensions` and `time_suspended` for every trial.

## Details

The trials are simulated by
[`tite_boin_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/tite_boin_simulate.md),
whose details describe the conduct of a trial. Once every patient has
completed the assessment the MTD is selected from the final data exactly
as in
[`sim_boin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin.md),
with
[`boin_select_mtd`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_select_mtd.md),
and a trial stopped for safety selects no dose.

When no patient is ever pending at a decision, for example with
`accrual = "fixed"` and `window` shorter than `1 / accrual_rate`, every
result agrees with that of
[`sim_boin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin.md)
under the same seed.

## References

Yuan, Y., Lin, R., Li, D., Nie, L. and Warren, K. E. (2018).
Time-to-Event Bayesian Optimal Interval Design to Accelerate Phase I
Trials. Clinical Cancer Research, 24(20), 4921-4930.

Lin, R. and Yuan, Y. (2020). Time-to-Event Model-Assisted Designs for
Dose-Finding Trials with Delayed Toxicity. Biostatistics, 21(4),
807-824.

## See also

[`tite_boin_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/tite_boin_simulate.md),
[`sim_tite_boin_multi`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_tite_boin_multi.md),
[`sim_boin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin.md)

## Examples

``` r
# A three month window with two patients a month
oc <- sim_tite_boin(
  target = 0.30,
  p_true = c(0.05, 0.15, 0.30, 0.45, 0.60),
  n_cohort = 10,
  cohort_size = 3,
  window = 3,
  accrual_rate = 2,
  n_trials = 500,
  seed = 123
)
oc
#> TITE-BOIN operating characteristics
#>   method          : single mean imputation (Yuan et al., 2018)
#>   target DLT rate : 0.3
#>   trials          : 500
#>   cohorts         : 10 of size 3
#>   max sample size : 30
#>   lambda_e / _d   : 0.2365 / 0.3585
#>   window          : 3
#>   accrual         : 2 patients per unit of time (exponential)
#> 
#>                   DL1  DL2  DL3  DL4  DL5 Total / No MTD
#> True DLT rate (%) 5.0 15.0 30.0 45.0 60.0               
#> MTD selected (%)  0.8 28.8 55.0 15.2  0.2            0.0
#> Patients treated  5.6 10.3  9.6  3.4  0.4           29.3
#> Patients with DLT 0.3  1.5  2.9  1.5  0.3            6.6
#> 
#> Doses above a true DLT rate of 0.3: DL4, DL5
#>   Patients treated there              : 13.1%
#>   Trials treating any patient there   : 41.6%
#>   Trials treating over 60% there      : 2.4%
#>   Trials treating over 80% there      : 0%
#>   Trials selecting an MTD there       : 15.4%
#> 
#> Trial duration and suspension of accrual
#>   Average trial duration              : 23.5
#>   Trials with accrual suspended       : 100%
#>   Average number of suspensions       : 3.9
#>   Average time suspended              : 6.5

# \donttest{
# The effective sample size of Lin and Yuan (2020)
oc_ess <- sim_tite_boin(
  target = 0.30,
  p_true = c(0.05, 0.15, 0.30, 0.45, 0.60),
  n_cohort = 12,
  cohort_size = 3,
  window = 3,
  accrual_rate = 2,
  method = "ess",
  n_trials = 10000,
  seed = 123
)
oc_ess
#> TITE-BOIN operating characteristics
#>   method          : effective sample size (Lin and Yuan, 2020)
#>   target DLT rate : 0.3
#>   trials          : 10000
#>   cohorts         : 12 of size 3
#>   max sample size : 36
#>   lambda_e / _d   : 0.2365 / 0.3585
#>   window          : 3
#>   accrual         : 2 patients per unit of time (exponential)
#> 
#>                   DL1  DL2  DL3  DL4  DL5 Total / No MTD
#> True DLT rate (%) 5.0 15.0 30.0 45.0 60.0               
#> MTD selected (%)  1.1 23.5 57.7 16.6  1.1            0.0
#> Patients treated  5.6 10.9 11.7  4.8  0.9           33.9
#> Patients with DLT 0.3  1.6  3.5  2.2  0.5            8.1
#> 
#> Doses above a true DLT rate of 0.3: DL4, DL5
#>   Patients treated there              : 16.7%
#>   Trials treating any patient there   : 55.9%
#>   Trials treating over 60% there      : 3%
#>   Trials treating over 80% there      : 0%
#>   Trials selecting an MTD there       : 17.7%
#> 
#> Trial duration and suspension of accrual
#>   Average trial duration              : 24.1
#>   Trials with accrual suspended       : 99.5%
#>   Average number of suspensions       : 2.9
#>   Average time suspended              : 4.8
# }
```
