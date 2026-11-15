# Operating Characteristics of the BOIN-ET Design

Simulate trials of the BOIN-ET design of Takeda, Taguri and Morita
(2018) under one scenario of DLT and response probabilities, select the
optimal biological dose (OBD) of each trial, and summarize how often
each dose is selected and how many patients, DLTs and responses occur at
each dose.

## Usage

``` r
sim_boin_et(
  target,
  target_resp,
  p_true,
  p_resp,
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

An object of class `boin_et_oc`, which is a list with components

- target, target_resp, p_true, p_resp:

  The targets and the scenario.

- n_doses, n_trials:

  Numbers of doses and of simulated trials.

- lambda1, lambda2, eta1:

  The boundaries, see
  [`boin_et_boundary`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_et_boundary.md).

- sel_percent:

  Percentage of trials selecting each dose as the OBD.

- percent_no_obd:

  Percentage of trials selecting no OBD.

- percent_early_stop:

  Percentage of trials stopped early.

- n_pts_dose, n_tox_dose, n_resp_dose:

  Average numbers of patients, DLTs and responses at each dose.

- total_n_pts, total_n_tox, total_n_resp:

  Their averages per trial.

- stop_reasons:

  Percentage of trials ending for each reason, see
  [`boin_et_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_et_simulate.md).

- est_method, obd_method:

  The OBD selection used.

- settings:

  List of the design parameters.

- call:

  The matched call.

- trials, obd:

  The simulated trials and the selected OBDs, when `keep_trials` is
  `TRUE`.

## Details

The trials are simulated by
[`boin_et_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_et_simulate.md)
and the OBD is selected by
[`boin_et_select_obd`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_et_select_obd.md);
see there for the rules. The defaults are those of
[`boinet::boinet()`](https://rdrr.io/pkg/boinet/man/boinet.html):
observed response rates, and the admissible dose with the highest
response rate up to the estimated MTD.

## References

Takeda, K., Taguri, M. and Morita, S. (2018). BOIN-ET: Bayesian optimal
interval design for dose finding based on both efficacy and toxicity
outcomes. Pharmaceutical Statistics, 17(4), 383-395.
[doi:10.1002/pst.1864](https://doi.org/10.1002/pst.1864)

Yamaguchi, Y., Takeda, K. and Maruo, K. (2025). boinet: An R Package for
Bayesian Optimal Interval Design for Dose Finding Based on Both Efficacy
and Toxicity Outcomes. The R Journal, 17(2), 15-34.
[doi:10.32614/RJ-2025-012](https://doi.org/10.32614/RJ-2025-012)

## See also

[`sim_boin_et_multi`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin_et_multi.md),
[`boin_et_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_et_simulate.md),
[`boin_et_select_obd`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_et_select_obd.md)

## Examples

``` r
# Scenario 1 of Takeda et al. (2018)
oc <- sim_boin_et(
  target = 0.30, target_resp = 0.60,
  p_true = c(0.05, 0.15, 0.25, 0.30, 0.35, 0.40),
  p_resp = c(0.05, 0.35, 0.55, 0.65, 0.75, 0.70),
  n_cohort = 12, cohort_size = 3, n_trials = 1000, seed = 123
)
oc
#> BOIN-ET operating characteristics
#>   target DLT rate      : 0.3
#>   target response rate : 0.6
#>   trials               : 1000
#>   cohorts              : 12 of size 3
#>   max sample size      : 36
#>   lambda1 / lambda2    : 0.14 / 0.35
#>   eta1                 : 0.48
#>   OBD selection        : max_resp with observed response rates
#> 
#>                        DL1  DL2  DL3  DL4  DL5  DL6 Total / No OBD
#> True DLT rate (%)      5.0 15.0 25.0 30.0 35.0 40.0               
#> True response rate (%) 5.0 35.0 55.0 65.0 75.0 70.0               
#> OBD selected (%)       1.4 15.1 48.0 26.5  7.3  0.9            0.8
#> Patients treated       3.5  8.0 14.4  7.6  2.0  0.3           35.9
#> Patients with DLT      0.2  1.3  3.7  2.4  0.7  0.1            8.4
#> Patients with response 0.2  2.8  7.9  5.0  1.5  0.2           17.7
#> 
#> Trials stopped early: 0.8%

# The same trials under the scoring utility
sim_boin_et(
  target = 0.30, target_resp = 0.60,
  p_true = c(0.05, 0.15, 0.25, 0.30, 0.35, 0.40),
  p_resp = c(0.05, 0.35, 0.55, 0.65, 0.75, 0.70),
  n_cohort = 12, cohort_size = 3, n_trials = 1000,
  obd_method = "utility_scoring", seed = 123
)$sel_percent
#>  DL1  DL2  DL3  DL4  DL5  DL6 
#>  2.1 10.3 39.2 31.8 13.1  2.7 
```
