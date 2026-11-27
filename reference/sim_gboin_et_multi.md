# Operating Characteristics of the gBOIN-ET Design in Several Scenarios

Run
[`sim_gboin_et`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_gboin_et.md)
for several scenarios of probabilities of toxicity and efficacy
categories with the same design, and collect the results in one table.

## Usage

``` r
sim_gboin_et_multi(
  target,
  target_resp,
  scenarios,
  score,
  score_resp,
  n_cohort,
  cohort_size,
  n_trials = 10000,
  start_dose = 1,
  n_earlystop = NULL,
  p_saf = NULL,
  p_tox = NULL,
  p_saf_resp = NULL,
  lambda1 = NULL,
  lambda2 = NULL,
  eta1 = NULL,
  cutoff_eli = 0.95,
  cutoff_eli_resp = 0.99,
  p_eli = NULL,
  p_eli_resp = NULL,
  resp_cor = 0.2,
  boinet_compatible = FALSE,
  prior = rep(1/6, 6),
  est_method = c("observed", "fp_linear", "fp_logistic"),
  obd_method = c("max_resp", "utility_weighted", "utility_truncated_linear",
    "utility_scoring"),
  min_obd_sample = 1,
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

  Numeric scalar. Target value of the mean normalized toxicity score,
  for example `0.47 / 1.5` for a target equivalent toxicity score of
  0.47 with a largest score of 1.5.

- target_resp:

  Numeric scalar. Target value of the mean normalized efficacy score,
  for example `0.95 / 1.5` for a target equivalent efficacy score of
  0.95 with a largest score of 1.5.

- scenarios:

  A list of scenarios, each a list with the components `p_grade` and
  `p_resp_grade`, matrices as the arguments of the same names of
  [`sim_gboin_et`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_gboin_et.md),
  and optionally `name`. The names of the list are used as scenario
  names. Every scenario must have the same numbers of categories and of
  doses.

- score:

  Numeric vector, one value per row of `p_grade`. Score of each toxicity
  category, for example `c(0, 0.5, 1, 1.5)` for grades 0 or 1, 2, 3
  and 4. The scores are divided by their maximum.

- score_resp:

  Numeric vector, one value per row of `p_resp_grade`. Score of each
  efficacy category, for example `c(0, 1, 1.5)` for stable or
  progressive disease, partial response and complete response. The
  scores are divided by their maximum.

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

- lambda1, lambda2, eta1:

  Numeric scalars or `NULL`. Boundaries of the mean toxicity score and
  of the mean efficacy score. Defaults to `NULL`, the boundaries of
  [`boin_et_boundary`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_et_boundary.md)
  for the targets. Give all three to use other values, such as those
  printed in an article.

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

  Numeric scalar strictly between -1 and 1. Correlation of the Gaussian
  copula that links the toxicity event and the efficacy event of a
  patient. Defaults to 0.2, the default of the 'boinet' package.

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

  Character scalar. Estimator of the mean efficacy scores: `"observed"`
  (the default) for the observed means, `"fp_linear"` for the linear
  fractional polynomial of Takeda et al. (2022), or `"fp_logistic"` for
  the fractional polynomial logistic regression of the 'boinet' package,
  which needs the 'mfp' package.

- obd_method:

  Character scalar. Selection rule: `"max_resp"` (the default) for the
  admissible dose with the highest estimated response rate up to the
  estimated MTD, or the dose that maximizes the utility
  `"utility_weighted"`, `"utility_truncated_linear"` or
  `"utility_scoring"` among the admissible doses.

- min_obd_sample:

  Integer scalar. Only doses given to at least this many patients can be
  selected. Defaults to 1. The simulations of Takeda et al. (2022) use
  6.

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

  The result of
  [`sim_gboin_et`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_gboin_et.md)
  for every scenario.

- summary_table:

  Data frame with the true mean scores, the selection percentages and
  the average numbers of patients and total scores of every scenario.

- scenario_names:

  The names of the scenarios.

- n_doses:

  The number of doses.

- call:

  The call.

## Details

Every scenario is simulated with the same seed, as in
[`sim_boin_et_multi`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin_et_multi.md).

## References

Takeda, K., Morita, S. and Taguri, M. (2022). gBOIN-ET: The generalized
Bayesian optimal interval design for optimal dose-finding accounting for
ordinal graded efficacy and toxicity in early clinical trials.
Biometrical Journal, 64(7), 1178-1191.
[doi:10.1002/bimj.202100263](https://doi.org/10.1002/bimj.202100263)

## See also

[`sim_gboin_et`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_gboin_et.md),
[`sim_boin_et_multi`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin_et_multi.md)

## Examples

``` r
# The same toxicity grades with an increasing and with a plateauing
# efficacy
p_grade <- rbind(
  c(0.83, 0.62, 0.51, 0.34),
  c(0.12, 0.18, 0.19, 0.16),
  c(0.04, 0.11, 0.14, 0.15),
  c(0.01, 0.09, 0.16, 0.35)
)
scenarios <- list(
  "Increasing" = list(
    p_grade = p_grade,
    p_resp_grade = rbind(
      c(0.90, 0.50, 0.30, 0.25),
      c(0.05, 0.35, 0.20, 0.20),
      c(0.05, 0.15, 0.50, 0.55)
    )
  ),
  "Plateau" = list(
    p_grade = p_grade,
    p_resp_grade = rbind(
      c(0.60, 0.30, 0.30, 0.30),
      c(0.30, 0.20, 0.20, 0.20),
      c(0.10, 0.50, 0.50, 0.50)
    )
  )
)
res <- sim_gboin_et_multi(
  target = 0.47 / 1.5, target_resp = 0.95 / 1.5, scenarios = scenarios,
  score = c(0, 0.5, 1, 1.5), score_resp = c(0, 1, 1.5),
  n_cohort = 10, cohort_size = 3, n_trials = 300, seed = 123
)
res
#> gBOIN-ET operating characteristics across 2 scenarios
#>   target toxicity score: 0.3133333
#>   target efficacy score: 0.6333333
#>   trials each          : 300
#>   max sample size      : 30
#>   OBD selection        : max_resp with observed efficacy scores
#> 
#>    Scenario                    Item  DL1  DL2  DL3  DL4 Total / No OBD
#>  Increasing True toxicity score (%)  7.7 22.3 31.7 50.3               
#>             True efficacy score (%)  8.3 38.3 63.3 68.3               
#>                    OBD selected (%)  2.7 32.0 57.3  4.3            3.7
#>                    Patients treated  4.4 10.6 13.3  1.3           29.6
#>                      Toxicity score  0.3  2.4  4.3  0.7            7.6
#>                      Efficacy score  0.4  4.0  8.5  0.9           13.8
#>     Plateau True toxicity score (%)  7.7 22.3 31.7 50.3               
#>             True efficacy score (%) 30.0 63.3 63.3 63.3               
#>                    OBD selected (%)  2.0 82.3 15.0  0.7            0.0
#>                    Patients treated  5.8 19.9  3.8  0.5           30.0
#>                      Toxicity score  0.4  4.5  1.3  0.2            6.4
#>                      Efficacy score  1.7 13.0  2.4  0.3           17.3
```
