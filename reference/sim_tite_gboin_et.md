# Operating Characteristics of the TITE-gBOIN-ET Design

Simulate trials of the TITE-gBOIN-ET design of Takeda, Yamaguchi, Taguri
and Morita (2023) under one scenario of probabilities of toxicity and
efficacy categories, select the optimal biological dose (OBD) of each
trial, and summarize how often each dose is selected, how many patients
are treated at each dose, with their total scores, and how long the
trials last.

## Usage

``` r
sim_tite_gboin_et(
  target,
  target_resp,
  p_grade,
  p_resp_grade,
  score,
  score_resp,
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
  lambda1 = NULL,
  lambda2 = NULL,
  eta1 = NULL,
  cutoff_eli = 0.95,
  cutoff_eli_resp = 0.99,
  p_eli = NULL,
  p_eli_resp = NULL,
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

- p_grade:

  Numeric matrix. True probability of each toxicity category (rows) at
  each dose (columns, in increasing dose order). Every column sums to 1.
  The first row is the category without a toxicity event.

- p_resp_grade:

  Numeric matrix with one column per dose. True probability of each
  efficacy category (rows). Every column sums to 1. The first row is the
  category without an efficacy event.

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

- window:

  Numeric scalar. Length of the toxicity assessment window.

- accrual_rate:

  Numeric scalar. Average number of patients arriving per unit of time,
  in the same unit as `window`. The `accrual` argument of the 'boinet'
  package is the average time between arrivals, its inverse.

- resp_window:

  Numeric scalar. Length of the efficacy assessment window. Defaults to
  `window`.

- accrual:

  Character scalar. Distribution of the times between arrivals:
  `"uniform"` (the default, as in 'boinet') on `(0, 2 / accrual_rate)`,
  `"exponential"` for a Poisson process, or `"fixed"` for arrivals
  exactly `1 / accrual_rate` apart.

- event_time:

  Character scalar. Distribution of the times to the toxicity event and
  to the efficacy event: `"weibull"` (the default), with the share of
  the events in the second half of the window set by `late_fraction` and
  `resp_late_fraction`, or `"uniform"` on `(0, window / p)` for an event
  of probability `p`. With `"weibull"` the first category of both
  outcomes, the one without an event, must have a positive probability
  at every dose.

- late_fraction:

  Numeric scalar strictly between 0 and 1. Share of the toxicity events
  that occur in the second half of the toxicity window under the Weibull
  distribution. Defaults to 0.5.

- resp_late_fraction:

  Numeric scalar strictly between 0 and 1. The same for the efficacy
  events. Defaults to 0.5.

- resp_cor:

  Numeric scalar strictly between -1 and 1. Correlation of the Gaussian
  copula that links the time to the toxicity event and the time to the
  efficacy event of a patient. Defaults to 0.2, the default of the
  'boinet' package.

- min_completed_share:

  Numeric scalar from 0 to 1. Share of the patients at the current dose
  who must have finished both assessments before the next decision; the
  decision waits until then. Defaults to 0.5, the rule of Takeda et
  al. (2023) and of the 'boinet' package.

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

An object of class `boin_et_oc`, as returned by
[`sim_gboin_et`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_gboin_et.md),
with `design` set to `"TITE-gBOIN-ET"` and one more component

- duration:

  Average trial duration, in the unit of `window`.

## Details

The trials are simulated by
[`tite_gboin_et_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/tite_gboin_et_simulate.md)
and the OBD is selected by
[`gboin_et_select_obd`](https://gosukehommaex.github.io/simFastBOIN/reference/gboin_et_select_obd.md)
from the data at the last decision, as in
[`boinet::tite.gboinet()`](https://rdrr.io/pkg/boinet/man/tite.gboinet.html);
see there for the rules. The result prints like that of
[`sim_gboin_et`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_gboin_et.md)
and enters
[`oc_results`](https://gosukehommaex.github.io/simFastBOIN/reference/oc_results.md)
and
[`oc_compare`](https://gosukehommaex.github.io/simFastBOIN/reference/oc_compare.md)
under the name TITE-gBOIN-ET.

Takeda et al. (2023) select, among the doses up to the MTD, the dose
with the highest observed mean efficacy score, which is the default
(`est_method = "observed"`, `obd_method = "max_resp"`), and exclude a
dose for its efficacy with a cutoff of 0.95 (`cutoff_eli_resp = 0.95`);
the default of 0.99 is that of the other functions of the BOIN-ET family
and of 'boinet'. `inst/reproduce/reproduce_Takeda2023.R` compares the
results with the tables of the article and lists the two points in which
the program of the authors differs from the design.

## References

Takeda, K., Yamaguchi, Y., Taguri, M. and Morita, S. (2023).
TITE-gBOIN-ET: Time-to-event generalized Bayesian optimal interval
design to accelerate dose-finding accounting for ordinal graded efficacy
and toxicity outcomes. Biometrical Journal, 65(7), 2200265.
[doi:10.1002/bimj.202200265](https://doi.org/10.1002/bimj.202200265)

## See also

[`sim_tite_gboin_et_multi`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_tite_gboin_et_multi.md),
[`tite_gboin_et_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/tite_gboin_et_simulate.md),
[`gboin_et_select_obd`](https://gosukehommaex.github.io/simFastBOIN/reference/gboin_et_select_obd.md),
[`sim_gboin_et`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_gboin_et.md),
[`sim_tite_boin_et`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_tite_boin_et.md)

## Examples

``` r
# Toxicity grades 0 or 1, 2, 3 and 4, and progressive disease, stable
# disease, partial response and complete response, at four doses, with a
# 30-day toxicity window, a 60-day efficacy window and one patient every
# 10 days on average
p_grade <- rbind(
  c(0.94, 0.79, 0.62, 0.50),
  c(0.05, 0.15, 0.20, 0.20),
  c(0.01, 0.05, 0.15, 0.25),
  c(0.00, 0.01, 0.03, 0.05)
)
p_resp_grade <- rbind(
  c(0.64, 0.45, 0.20, 0.05),
  c(0.30, 0.40, 0.40, 0.15),
  c(0.05, 0.10, 0.20, 0.35),
  c(0.01, 0.05, 0.20, 0.45)
)
oc <- sim_tite_gboin_et(
  target = 0.47 / 1.5, target_resp = 1.75 / 3,
  p_grade = p_grade, p_resp_grade = p_resp_grade,
  score = c(0, 0.5, 1, 1.5), score_resp = c(0, 0.25, 1, 3),
  n_cohort = 12, cohort_size = 3, window = 30, resp_window = 60,
  accrual_rate = 1 / 10, n_trials = 500, seed = 123
)
oc
#> TITE-gBOIN-ET operating characteristics
#>   target toxicity score: 0.3133333
#>   target efficacy score: 0.5833333
#>   trials               : 500
#>   cohorts              : 12 of size 3
#>   max sample size      : 36
#>   lambda1 / lambda2    : 0.1413 / 0.3733
#>   eta1                 : 0.46
#>   OBD selection        : max_resp with observed efficacy scores
#> 
#>                         DL1  DL2  DL3  DL4 Total / No OBD
#> True toxicity score (%) 2.3  9.3 19.7 28.3               
#> True efficacy score (%) 5.2 11.7 30.0 57.9               
#> OBD selected (%)        0.6  1.0 14.2 84.2            0.0
#> Patients treated        3.1  4.2  9.5 19.2           36.0
#> Toxicity score          0.1  0.4  1.8  5.5            7.8
#> Efficacy score          0.2  0.5  2.9 11.1           14.7
#> 
#> Trials stopped early: 0%
#> Average trial duration: 503.6
```
