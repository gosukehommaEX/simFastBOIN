# Operating Characteristics of the BOIN-ET Design Across Scenarios

Run
[`sim_boin_et`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin_et.md)
for several scenarios of DLT and response probabilities with the same
design and collect the results in one table.

## Usage

``` r
sim_boin_et_multi(
  target,
  target_resp,
  scenarios,
  n_cohort,
  cohort_size,
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
  resp_cor = 0.2,
  resp_cor_type = c("copula", "binary"),
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

- resp_cor:

  Numeric scalar strictly between -1 and 1. Correlation between the DLT
  and the response of a patient. Defaults to 0.2, the default of the
  'boinet' package.

- resp_cor_type:

  Character scalar. How `resp_cor` links the two outcomes: `"copula"`
  (the default) uses it as the correlation of a Gaussian copula, which
  is how the 'boinet' package generates its data; `"binary"` uses it as
  the correlation of the two binary outcomes, as in the simulations of
  Takeda et al. (2018).

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

An object of class `boin_et_oc_multi`, which is a list with components

- results:

  Named list of `boin_et_oc` objects, one per scenario.

- summary_table:

  Data frame with six rows per scenario, as shown by
  [`print.boin_et_oc_multi`](https://gosukehommaex.github.io/simFastBOIN/reference/print.boin_et_oc_multi.md).

- scenario_names:

  Character vector of scenario names.

- n_doses:

  Number of doses.

- call:

  The matched call.

## Details

Every scenario is simulated with the same `seed`.

## See also

[`sim_boin_et`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin_et.md)

## Examples

``` r
scenarios <- list(
  list(name = "Monotone",
       p_true = c(0.05, 0.10, 0.15, 0.20, 0.25, 0.40),
       p_resp = c(0.05, 0.10, 0.25, 0.40, 0.65, 0.70)),
  list(name = "Peaked",
       p_true = c(0.05, 0.15, 0.25, 0.35, 0.45, 0.55),
       p_resp = c(0.30, 0.70, 0.60, 0.50, 0.40, 0.30))
)
sim_boin_et_multi(
  target = 0.30, target_resp = 0.60, scenarios = scenarios,
  n_cohort = 12, cohort_size = 3, n_trials = 500, seed = 123
)
#> BOIN-ET operating characteristics across 2 scenarios
#>   target DLT rate      : 0.3
#>   target response rate : 0.6
#>   trials each          : 500
#>   max sample size      : 36
#>   OBD selection        : max_resp with observed response rates
#> 
#>  Scenario                   Item  DL1  DL2  DL3  DL4  DL5  DL6 Total / No OBD
#>  Monotone      True DLT rate (%)  5.0 10.0 15.0 20.0 25.0 40.0               
#>           True response rate (%)  5.0 10.0 25.0 40.0 65.0 70.0               
#>                 OBD selected (%)  2.4  2.0  7.8 20.6 59.4  7.6            0.2
#>                 Patients treated  3.1  3.9  5.3  9.2 12.5  2.0           36.0
#>                Patients with DLT  0.2  0.4  0.8  1.8  3.1  0.8            7.1
#>           Patients with response  0.1  0.3  1.3  3.8  8.3  1.4           15.3
#>    Peaked      True DLT rate (%)  5.0 15.0 25.0 35.0 45.0 55.0               
#>           True response rate (%) 30.0 70.0 60.0 50.0 40.0 30.0               
#>                 OBD selected (%)  2.6 81.6 12.4  3.2  0.2  0.0            0.0
#>                 Patients treated  5.5 25.5  3.6  1.0  0.3  0.1           36.0
#>                Patients with DLT  0.3  3.9  0.9  0.4  0.1  0.1            5.6
#>           Patients with response  1.6 18.0  2.2  0.5  0.1  0.0           22.4
```
