# Operating Characteristics of the BF-BOIN-ET Design Across Scenarios

Run
[`sim_bf_boin_et`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_bf_boin_et.md)
for several scenarios of DLT and response probabilities with the same
design and collect the results in one table.

## Usage

``` r
sim_bf_boin_et_multi(
  target,
  target_resp,
  scenarios,
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
  psi01 = 0,
  psi10 = 100,
  psi11 = 60,
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

- scenarios:

  List of scenarios. Each element is a list with the components `p_true`
  and `p_resp`, and optionally `name`. Names of the list are used when
  `name` is absent. All scenarios must have the same number of doses.

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

- psi00, psi01, psi10, psi11:

  Numeric scalars. Scores of a patient without DLT or response, with a
  DLT and no response, with a response and no DLT, and with both, which
  define the utility used to choose the backfill dose and the OBD.
  Default to 40, 0, 100 and 60.

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

An object of class `bf_boin_et_oc_multi`, which is a list with
components

- results:

  Named list of `bf_boin_et_oc` objects, one per scenario.

- summary_table:

  Data frame with seven rows per scenario, as shown by
  [`print.bf_boin_et_oc_multi`](https://gosukehommaex.github.io/simFastBOIN/reference/print.bf_boin_et_oc_multi.md).

- scenario_names:

  Character vector of scenario names.

- n_doses:

  Number of doses.

- call:

  The matched call.

## Details

Every scenario is simulated with the same `seed`.

## See also

[`sim_bf_boin_et`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_bf_boin_et.md)

## Examples

``` r
scenarios <- list(
  "Example" = list(p_true = c(0.03, 0.05, 0.20, 0.22, 0.45),
                   p_resp = c(0.05, 0.10, 0.50, 0.68, 0.70)),
  "Plateau" = list(p_true = c(0.05, 0.10, 0.15, 0.30, 0.45),
                   p_resp = c(0.10, 0.30, 0.45, 0.50, 0.52))
)
res <- sim_bf_boin_et_multi(
  target = 0.30, target_resp = 0.25, scenarios = scenarios,
  n_cohort = 10, cohort_size = 3, window = 1, accrual_rate = 3,
  n_trials = 200, seed = 123
)
res
#> BF-BOIN-ET operating characteristics across 2 scenarios
#>   target DLT rate      : 0.3
#>   target response rate : 0.25
#>   trials each          : 200
#>   cohorts              : 10 of size 3
#> 
#>  Scenario                   Item  DL1  DL2  DL3  DL4  DL5 Total / No OBD
#>   Example      True DLT rate (%)  3.0  5.0 20.0 22.0 45.0               
#>           True response rate (%)  5.0 10.0 50.0 68.0 70.0               
#>                 OBD selected (%)  4.0  2.5 42.0 47.0  4.5            0.0
#>                 Patients treated  4.2  6.0 11.4  9.3  4.8           35.7
#>              Patients backfilled  0.8  1.4  4.4  2.5  0.2            9.2
#>                Patients with DLT  0.1  0.3  2.3  2.1  2.3            7.0
#>           Patients with response  0.2  0.5  5.8  6.4  3.5           16.4
#>   Plateau      True DLT rate (%)  5.0 10.0 15.0 30.0 45.0               
#>           True response rate (%) 10.0 30.0 45.0 50.0 52.0               
#>                 OBD selected (%)  6.0 17.5 51.0 22.0  3.0            0.5
#>                 Patients treated  5.4  9.0 11.8  9.3  3.5           38.9
#>              Patients backfilled  1.8  4.4  4.5  1.5  0.2           12.3
#>                Patients with DLT  0.3  0.9  1.9  2.9  1.6            7.5
#>           Patients with response  0.6  2.6  5.5  4.7  2.0           15.4
```
