# Operating Characteristics of the U-BOIN Design

Simulate a trial run with the utility-based BOIN (U-BOIN) phase I/II
design of Zhou, Lee and Yuan (2019) many times under one scenario of
true DLT and response probabilities, and summarize how often each dose
is selected as the optimal biological dose (OBD), how many patients are
treated at each dose and how many DLTs and responses are observed.

## Usage

``` r
sim_uboin(
  tox_upper,
  resp_lower,
  p_true = NULL,
  p_resp = NULL,
  n_cohort,
  cohort_size,
  n_trials = 10000,
  target = NULL,
  resp_assoc = 0,
  p_outcome = NULL,
  utility = c(0, 30, 50, 100),
  n_stage1 = 12,
  n_stop = NULL,
  start_dose = 1,
  p_saf = NULL,
  p_tox = NULL,
  cutoff_eli = 0.95,
  cutoff_eli_resp = 0.9,
  prior_tox = c(1, 1),
  prior_resp = c(1, 1),
  keep_trials = FALSE,
  verbose = FALSE,
  seed = 123
)
```

## Arguments

- tox_upper:

  Numeric scalar. Highest acceptable DLT probability, for example 0.30.

- resp_lower:

  Numeric scalar. Lowest acceptable response probability, for example
  0.20.

- p_true:

  Numeric vector or `NULL`. True DLT probability of each dose.

- p_resp:

  Numeric vector or `NULL`. True response probability of each dose.

- n_cohort:

  Integer scalar. Number of cohorts in a trial; the maximum sample size
  is `n_cohort * cohort_size`.

- cohort_size:

  Integer scalar. Number of patients per cohort.

- n_trials:

  Integer scalar. Number of trials to simulate. Defaults to 10000.

- target:

  Numeric scalar or `NULL`. Target DLT probability of the BOIN design in
  stage I. Defaults to `NULL`, which uses `tox_upper - 0.05`.

- resp_assoc:

  Numeric scalar. Association parameter of the Gumbel model that joins
  the DLT and the response of a patient; 0 makes them independent.
  Defaults to 0.

- p_outcome:

  Numeric matrix or `NULL`. Probabilities of the four outcomes of a
  patient at each dose, one row per dose, in the columns (no response,
  DLT), (no response, no DLT), (response, DLT) and (response, no DLT).
  When given, `p_true`, `p_resp` and `resp_assoc` must be `NULL` or are
  not used. Defaults to `NULL`.

- utility:

  Numeric vector of length four. Utility of the four outcomes, in the
  order of `p_outcome`. Defaults to `c(0, 30, 50, 100)`, example 1 of
  the article.

- n_stage1:

  Integer scalar. Stage II starts once a dose has this many patients (s1
  of the article). Defaults to 12.

- n_stop:

  Integer scalar or `NULL`. The trial ends once a dose has this many
  patients (s2 of the article). Defaults to `NULL`, the maximum sample
  size.

- start_dose:

  Integer scalar. Dose level of the first cohort. Defaults to 1.

- p_saf:

  Numeric scalar. Highest DLT probability deemed subtherapeutic, for the
  boundaries of stage I. Defaults to `0.6 * target`.

- p_tox:

  Numeric scalar. Lowest DLT probability deemed overly toxic, for the
  boundaries of stage I. Defaults to `1.4 * target`.

- cutoff_eli:

  Numeric scalar. A dose is toxic when the posterior probability that
  its DLT probability exceeds `tox_upper` is above this value. Defaults
  to 0.95.

- cutoff_eli_resp:

  Numeric scalar. A dose is futile when the posterior probability that
  its response probability is below `resp_lower` is above this value.
  Defaults to 0.90.

- prior_tox:

  Numeric vector of length two. Beta prior of the DLT probability in the
  toxicity rule. Defaults to `c(1, 1)`.

- prior_resp:

  Numeric vector of length two. Beta prior of the response probability
  in the futility rule. Defaults to `c(1, 1)`.

- keep_trials:

  Logical scalar. Keep the trial by trial data in the result. Defaults
  to `FALSE`.

- verbose:

  Logical scalar. Report progress while the simulation runs. Defaults to
  `FALSE`.

- seed:

  Integer scalar or `NULL`. Random seed. The state of the calling
  session is restored on exit. Defaults to 123.

## Value

An object of class `uboin_oc`, which is a list with components

- sel_percent:

  Percentage of trials selecting each dose as the OBD.

- percent_no_obd:

  Percentage of trials without an OBD, which the article calls the
  percentage of early stopping.

- n_pts_dose, n_tox_dose, n_resp_dose:

  Average numbers of patients, DLTs and responses at each dose.

- total_n_pts, total_n_tox, total_n_resp:

  Their totals.

- true_utility:

  True mean utility of each dose, `p_outcome %*% utility`.

- percent_stage2:

  Percentage of trials that reached stage II.

- stop_reason_percent:

  Percentage of trials by reason for ending.

- trials, obd:

  The trials and their OBDs when `keep_trials` is `TRUE`.

together with the design parameters and the call.

## Details

The trials are simulated with
[`uboin_simulate`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/uboin_simulate.md),
whose help page describes the rules, and the OBD is selected with
[`uboin_select_obd`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/uboin_select_obd.md):
the admissible dose with the largest posterior mean utility at the end
of the trial. A trial stopped because the lowest dose was eliminated or
no dose was admissible selects none.

With `resp_assoc = 0.2`, cohorts of 3, `n_cohort = 18`,
`tox_upper = 0.30`, `resp_lower = 0.20` and the other defaults, the
results agree with Table 3 of the article and Table S9 of its supporting
information within simulation error
(`inst/reproduce/reproduce_Zhou2019.R`).

## References

Zhou, Y., Lee, J. J. and Yuan, Y. (2019). A utility-based Bayesian
optimal interval (U-BOIN) phase I/II design to identify the optimal
biological dose for targeted and immune therapies. Statistics in
Medicine, 38(28), S5299-S5316.
[doi:10.1002/sim.8361](https://doi.org/10.1002/sim.8361)

## See also

[`sim_uboin_multi`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_uboin_multi.md),
[`uboin_simulate`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/uboin_simulate.md),
[`uboin_select_obd`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/uboin_select_obd.md),
[`uboin_utility_table`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/uboin_utility_table.md)

## Examples

``` r
# Scenario 1 of Simulation A of the article
oc <- sim_uboin(
  tox_upper = 0.30,
  resp_lower = 0.20,
  p_true = c(0.02, 0.15, 0.30, 0.45, 0.60),
  p_resp = c(0.20, 0.65, 0.65, 0.65, 0.65),
  n_cohort = 18,
  cohort_size = 3,
  resp_assoc = 0.2,
  n_trials = 500,
  seed = 123
)
oc
#> U-BOIN operating characteristics
#>   DLT upper limit      : 0.3
#>   response lower limit : 0.2
#>   stage I target       : 0.25
#>   trials               : 500
#>   cohorts              : 18 of size 3
#>   max sample size      : 54
#>   stage II from / stop : 12 / 54 patients at a dose
#>   lambda_e / lambda_d  : 0.1968 / 0.2984
#>   utilities            : 0, 30, 50, 100
#>   cutoffs (DLT / resp) : 0.95 / 0.9
#> 
#>                         DL1  DL2  DL3  DL4  DL5 Total / No OBD
#> True DLT rate (%)       2.0 15.0 30.0 45.0 60.0               
#> True response rate (%) 20.0 65.0 65.0 65.0 65.0               
#> True utility           43.3 69.0 62.5 56.0 49.6               
#> OBD selected (%)        2.4 77.2 19.2  0.8  0.4            0.0
#> Patients treated        6.3 31.2 12.8  3.3  0.5           54.0
#> Patients with DLT       0.1  4.7  3.9  1.5  0.3           10.5
#> Patients with response  1.2 20.4  8.3  2.1  0.3           32.3
#> 
#> Trials reaching stage II: 100%
#> Reason for ending (%):
#>        max_sample_size                 n_stop lowest_dose_eliminated 
#>                    100                      0                      0 
#>     no_admissible_dose 
#>                      0 
```
