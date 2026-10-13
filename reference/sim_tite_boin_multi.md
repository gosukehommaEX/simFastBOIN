# Operating Characteristics of a TITE-BOIN Design Across Several Scenarios

Run
[`sim_tite_boin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_tite_boin.md)
under a list of dose-toxicity scenarios that share the same design, and
collect the results into one table.

## Usage

``` r
sim_tite_boin_multi(
  target,
  scenarios,
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

- scenarios:

  A list of scenarios. Each element is a list with a `name` and a
  `p_true` vector, or a plain numeric vector of true DLT probabilities
  in which case the list names are used as scenario names. Every
  scenario must have the same number of doses.

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

  Logical scalar. Keep the trial by trial data of every scenario.
  Defaults to `FALSE`.

- verbose:

  Logical scalar. Report progress while the simulations run. Defaults to
  `FALSE`.

- seed:

  Integer scalar or `NULL`. Random seed, applied to every scenario, so
  that a scenario simulated here matches the same scenario simulated on
  its own with
  [`sim_tite_boin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_tite_boin.md).
  Defaults to 123.

## Value

An object of class `c("tite_boin_oc_multi", "boin_oc_multi")`, which is
a list with components

- results:

  Named list of `tite_boin_oc` objects, one per scenario.

- summary_table:

  Data frame collecting all scenarios, with the four rows of
  [`sim_boin_multi`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin_multi.md)
  and two more per scenario for the average trial duration and the
  percentage of trials with accrual suspended, both in the last column.

- scenario_names:

  Character vector of scenario names.

- n_doses:

  Number of dose levels.

- call:

  The matched call.

## See also

[`sim_tite_boin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_tite_boin.md),
[`sim_boin_multi`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin_multi.md)

## Examples

``` r
scenarios <- list(
  list(name = "MTD at dose 3", p_true = c(0.05, 0.15, 0.30, 0.45, 0.60)),
  list(name = "All toxic",     p_true = c(0.35, 0.45, 0.55, 0.65, 0.75))
)

oc <- sim_tite_boin_multi(
  target = 0.30,
  scenarios = scenarios,
  n_cohort = 10,
  cohort_size = 3,
  window = 3,
  accrual_rate = 2,
  n_trials = 200,
  seed = 123
)
oc
#> TITE-BOIN operating characteristics across 2 scenarios
#>   method          : single mean imputation (Yuan et al., 2018)
#>   target DLT rate : 0.3
#>   trials each     : 200
#>   max sample size : 30
#>   window          : 3
#>   accrual         : 2 patients per unit of time (exponential)
#> 
#>       Scenario                         Item  DL1  DL2  DL3  DL4  DL5
#>  MTD at dose 3            True DLT rate (%)  5.0 15.0 30.0 45.0 60.0
#>                            MTD selected (%)  2.0 30.5 56.0 11.5  0.0
#>                            Patients treated  5.8 10.6  9.1  3.2  0.4
#>                           Patients with DLT  0.3  1.6  2.8  1.5  0.3
#>                       Trial duration (mean)                         
#>                Accrual suspended (% trials)                         
#>      All toxic            True DLT rate (%) 35.0 45.0 55.0 65.0 75.0
#>                            MTD selected (%) 75.0  8.0  1.0  0.0  0.0
#>                            Patients treated 16.4  4.1  0.4  0.0  0.0
#>                           Patients with DLT  5.7  1.9  0.2  0.0  0.0
#>                       Trial duration (mean)                         
#>                Accrual suspended (% trials)                         
#>  Total / No MTD
#>                
#>             0.0
#>            29.2
#>             6.5
#>            23.5
#>           100.0
#>                
#>            16.0
#>            21.0
#>             7.8
#>            15.6
#>            98.0
```
