# Operating Characteristics of the TITE-BOIN-ET Design Across Scenarios

Run
[`sim_tite_boin_et`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_tite_boin_et.md)
for several scenarios of DLT and response probabilities with the same
design and collect the results in one table.

## Usage

``` r
sim_tite_boin_et_multi(
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
  min_completed_share = 0.5,
  n_trials = 10000,
  start_dose = 1,
  n_earlystop = NULL,
  p_saf = NULL,
  p_tox = NULL,
  p_saf_resp = NULL,
  cutoff_eli = 0.95,
  cutoff_eli_resp = 0.99,
  p_eli = NULL,
  p_eli_resp = NULL,
  boinet_compatible = FALSE,
  prior = rep(1/6, 6),
  est_method = c("observed", "multi_iso", "fp_logistic"),
  obd_method = c("max_resp", "utility_weighted", "utility_truncated_linear",
    "utility_scoring"),
  w1 = 0.33,
  w2 = 1.09,
  t_low = NULL,
  t_upp = NULL,
  e_low = NULL,
  e_upp = NULL,
  psi00 = 40,
  psi11 = 60,
  keep_trials = FALSE,
  verbose = FALSE,
  seed = 123
)
```

## Arguments

- target:

  Numeric scalar. Target DLT probability, for example 0.30.

- target_resp:

  Numeric scalar. Target response probability, for example 0.60.

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

- min_completed_share:

  Numeric scalar from 0 to 1. Share of the patients at the current dose
  who must have finished both assessments before the next decision; the
  decision waits until then. Defaults to 0.5, the rule of the 'boinet'
  package. Takeda et al. (2020) did not wait in their simulations, which
  is `min_completed_share = 0`, and recommended such a rule when the
  accrual is fast.

- n_trials:

  Integer scalar. Number of trials to simulate. Defaults to 10000.

- start_dose:

  Integer scalar. Dose level for the first cohort. Defaults to 1.

- n_earlystop:

  Integer scalar or `NULL`. The trial stops once this many patients have
  been treated at any dose. Defaults to `NULL`, the maximum sample size,
  which is the default of the 'boinet' package.

- p_saf:

  Numeric scalar. Highest DLT probability deemed subtherapeutic.
  Defaults to `0.1 * target`, the value recommended by Takeda et al.
  (2018), which is far below the `0.6 * target` of the BOIN design.

- p_tox:

  Numeric scalar. Lowest DLT probability deemed overly toxic. Defaults
  to `1.4 * target`.

- p_saf_resp:

  Numeric scalar. Highest response probability deemed subtherapeutic.
  Defaults to `0.6 * target_resp`.

- cutoff_eli:

  Numeric scalar. A dose is no longer admissible when the posterior
  probability that its DLT rate is below `p_eli` falls below
  `1 - cutoff_eli`. Defaults to 0.95.

- cutoff_eli_resp:

  Numeric scalar. A dose is no longer admissible when the posterior
  probability that its response rate is above `p_eli_resp` falls below
  `1 - cutoff_eli_resp`. Defaults to 0.99.

- p_eli:

  Numeric scalar or `NULL`. DLT probability of the rule above. Defaults
  to `NULL`, which uses `target` as the 'boinet' package does. Takeda et
  al. (2018) used 0.40 with a cutoff of 0.90.

- p_eli_resp:

  Numeric scalar or `NULL`. Response probability of the rule above.
  Defaults to `NULL`, which uses `p_saf_resp` as the 'boinet' package
  does. Takeda et al. (2018) used 0.20 with a cutoff of 0.90.

- boinet_compatible:

  Logical scalar. Apply the rules of the 'boinet' package (version
  1.6.0) exactly. Defaults to `FALSE`, which keeps to Takeda et al.
  (2018), who do not allow dose skipping during escalation; see Details.

- prior:

  Numeric vector of length 6. Prior probabilities of the six hypotheses.
  Defaults to `rep(1 / 6, 6)`, the noninformative choice of the article.
  The third hypothesis, under which every decision is correct, does not
  enter the criterion.

- est_method:

  Character scalar. Estimator of the response rates: `"observed"` (the
  default) for the observed rates, `"multi_iso"` for the model-averaged
  unimodal isotonic regression of Lin and Yin (2017), or `"fp_logistic"`
  for the fractional polynomial logistic regression of Takeda et al.
  (2018), which needs the 'mfp' package.

- obd_method:

  Character scalar. Selection rule: `"max_resp"` (the default) for the
  admissible dose with the highest estimated response rate up to the
  estimated MTD, or the dose that maximizes the utility
  `"utility_weighted"`, `"utility_truncated_linear"` or
  `"utility_scoring"` among the admissible doses.

- w1, w2:

  Numeric scalars. Weights of the weighted utility
  `q - w1 * p - w2 * p * (p > p_tox)`, where `p` and `q` are the
  estimated DLT and response rates. Default to 0.33 and 1.09.

- t_low, t_upp:

  Numeric scalars or `NULL`. DLT rates between which the truncated
  linear utility falls from one to zero. Default to `p_saf` and `p_tox`.

- e_low, e_upp:

  Numeric scalars or `NULL`. Response rates between which the truncated
  linear utility rises from zero to one. Default to `p_saf_resp / 2` and
  `target_resp`.

- psi00, psi11:

  Numeric scalars. Scores of the scoring utility for a patient with
  neither DLT nor response and for a patient with both, the scores of a
  response without DLT and of a DLT without response being 100 and 0.
  Default to 40 and 60.

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

An object of class `boin_et_oc_multi`, as returned by
[`sim_boin_et_multi`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin_et_multi.md),
whose `results` are the `boin_et_oc` objects of
[`sim_tite_boin_et`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_tite_boin_et.md).

## Details

Every scenario is simulated with the same `seed`.

## See also

[`sim_tite_boin_et`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_tite_boin_et.md)

## Examples

``` r
scenarios <- list(
  "Scenario 1" = list(p_true = c(0.05, 0.15, 0.25, 0.30, 0.35, 0.40),
                      p_resp = c(0.05, 0.35, 0.55, 0.65, 0.75, 0.70)),
  "Scenario 7" = list(p_true = c(0.05, 0.15, 0.25, 0.35, 0.45, 0.55),
                      p_resp = c(0.30, 0.70, 0.60, 0.50, 0.40, 0.30))
)
res <- sim_tite_boin_et_multi(
  target = 0.30, target_resp = 0.60, scenarios = scenarios,
  n_cohort = 12, cohort_size = 3, window = 30, resp_window = 45,
  accrual_rate = 1 / 10, n_trials = 300, seed = 123
)
res
#> TITE-BOIN-ET operating characteristics across 2 scenarios
#>   target DLT rate      : 0.3
#>   target response rate : 0.6
#>   trials each          : 300
#>   max sample size      : 36
#>   OBD selection        : max_resp with observed response rates
#> 
#>    Scenario                   Item  DL1  DL2  DL3  DL4  DL5  DL6 Total / No OBD
#>  Scenario 1      True DLT rate (%)  5.0 15.0 25.0 30.0 35.0 40.0               
#>             True response rate (%)  5.0 35.0 55.0 65.0 75.0 70.0               
#>                   OBD selected (%)  0.7 18.7 43.0 26.7  7.3  0.3            3.3
#>                   Patients treated  4.2  9.5 13.3  6.8  1.7  0.1           35.6
#>                  Patients with DLT  0.2  1.4  3.4  2.1  0.6  0.0            7.7
#>             Patients with response  0.2  3.4  7.4  4.4  1.3  0.1           16.8
#>  Scenario 7      True DLT rate (%)  5.0 15.0 25.0 35.0 45.0 55.0               
#>             True response rate (%) 30.0 70.0 60.0 50.0 40.0 30.0               
#>                   OBD selected (%)  4.0 78.3 14.3  2.7  0.7  0.0            0.0
#>                   Patients treated  6.6 21.9  5.4  1.5  0.5  0.1           36.0
#>                  Patients with DLT  0.3  3.4  1.3  0.6  0.2  0.1            5.9
#>             Patients with response  2.0 15.5  3.3  0.8  0.2  0.0           21.9
```
