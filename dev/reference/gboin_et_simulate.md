# Simulate gBOIN-ET Trials

Run the gBOIN-ET dose-finding algorithm of Takeda, Morita and Taguri
(2022) over many simulated trials, for toxicity and efficacy outcomes
with several categories that each carry a score, and return the
patients, their categories and the total scores at every dose. This is
the simulation engine behind
[`sim_gboin_et`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_gboin_et.md),
exposed for users who want the raw trial data.

## Usage

``` r
gboin_et_simulate(
  target,
  target_resp,
  p_grade,
  p_resp_grade,
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
  [`boin_et_boundary`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/boin_et_boundary.md)
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

- seed:

  Integer scalar or `NULL`. Random seed. When `NULL` the current state
  of the random number generator is used. The state of the calling
  session is restored on exit. Defaults to 123.

## Value

An object of class `boin_et_trials`, which is a list with components

- n_pts:

  Integer matrix, `n_trials` by number of doses, of patients treated.

- n_tox:

  Numeric matrix of the same shape: the total of the normalized toxicity
  scores of the patients treated at each dose.

- n_resp:

  Numeric matrix of the same shape: the total of the normalized efficacy
  scores.

- n_grade:

  Integer array, `n_trials` by doses by toxicity categories, of the
  patients in each category.

- n_resp_grade:

  Integer array, `n_trials` by doses by efficacy categories. The
  categories take the row names of `p_grade` and `p_resp_grade` when
  every row has one, and `G1`, `G2`, ... and `R1`, `R2`, ... otherwise.

- cohorts_used:

  Integer vector of cohorts enrolled in each trial.

- stop_reason:

  Character vector giving why each trial ended.

- early_stop:

  Logical vector, `TRUE` for the trials stopped early, which select no
  OBD.

- design:

  `"gBOIN-ET"`.

- boundary:

  The boundaries used, an object of class `boin_et_boundary`.

- settings:

  List of the design parameters, with the true mean normalized scores of
  each dose in `p_true` and `p_resp`.

## Details

The gBOIN-ET design is the BOIN-ET design for ordinal outcomes. Every
toxicity category and every efficacy category carries a score, such as
the equivalent toxicity score and the equivalent efficacy score of the
article. The scores divided by their maximum lie between 0 and 1 and are
treated like fractional DLTs and responses (quasi-Bernoulli outcomes):
the mean normalized scores of the patients at the current dose take the
place of the observed DLT and response rates in the rules of
[`boin_et_simulate`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/boin_et_simulate.md).
With `pt` and `pe` the mean toxicity and efficacy scores at the current
dose, the design

- escalates when `pt <= lambda1` and `pe <= eta1`,

- stays when `pt < lambda2` and `pe > eta1`,

- de-escalates when `pt >= lambda2`,

- and otherwise moves to the dose with the highest mean efficacy score
  among the current dose and its neighbors, an untried higher dose being
  preferred and ties being broken at random.

A dose is excluded when the posterior probability that its mean toxicity
score is below `p_eli` falls below `1 - cutoff_eli`, or when the
posterior probability that its mean efficacy score is above `p_eli_resp`
falls below `1 - cutoff_eli_resp`; the posterior distribution of a mean
score with total `x` among `n` patients is Beta(`1 + x`, `1 + n - x`).
How the next dose is found among the admissible doses, when the trial
ends, and how the rules differ from those of 'boinet' unless
`boinet_compatible` is `TRUE` is described in
[`boin_et_simulate`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/boin_et_simulate.md);
with `boinet_compatible = TRUE` the rules are those of
[`boinet::gboinet()`](https://rdrr.io/pkg/boinet/man/gboinet.html).

`target`, `target_resp` and the other design parameters are on the
normalized scale. The boundaries of
[`boin_et_boundary`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/boin_et_boundary.md)
come from a grid in steps of 0.01 that starts at `p_saf`, `target` and
`p_saf_resp`, so that they depend on the decimals of the targets: the
targets `0.47 / 1.5` and `0.95 / 1.5` give 0.1413, 0.3733 and 0.50, the
values used by Takeda et al. (2022), who print 0.14, 0.37 and 0.50. The
mean scores are not multiples of `1 / n`, so boundaries rounded to two
decimals can give other decisions; `lambda1`, `lambda2` and `eta1` take
the boundaries of an article directly.

A patient has a toxicity event when the toxicity category is not the
first one, and an efficacy event likewise. The two events are linked by
a Gaussian copula with correlation `resp_cor`, and the category of a
patient with an event is drawn from the probabilities of the categories
above the first, independently of everything else, which is the joint
distribution of
[`boinet::gboinet()`](https://rdrr.io/pkg/boinet/man/gboinet.html). Each
patient draws one uniform variate for the toxicity event, one for the
efficacy event, and one more for each outcome with more than two
categories; a tie under the rule for the neighboring doses draws one
more. With two categories of scores 0 and 1 for both outcomes the trials
are, for the same seed, those of
[`boin_et_simulate`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/boin_et_simulate.md).
The 'boinet' package draws its random numbers differently, so the trials
agree with it in distribution, not one by one.

Possible values of `stop_reason` are those of
[`boin_et_simulate`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/boin_et_simulate.md).

## References

Takeda, K., Morita, S. and Taguri, M. (2022). gBOIN-ET: The generalized
Bayesian optimal interval design for optimal dose-finding accounting for
ordinal graded efficacy and toxicity in early clinical trials.
Biometrical Journal, 64(7), 1178-1191.
[doi:10.1002/bimj.202100263](https://doi.org/10.1002/bimj.202100263)

Yamaguchi, Y., Takeda, K. and Maruo, K. (2025). boinet: An R Package for
Bayesian Optimal Interval Design for Dose Finding Based on Both Efficacy
and Toxicity Outcomes. The R Journal, 17(2), 15-34.
[doi:10.32614/RJ-2025-012](https://doi.org/10.32614/RJ-2025-012)

## See also

[`sim_gboin_et`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_gboin_et.md),
[`gboin_et_select_obd`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/gboin_et_select_obd.md),
[`boin_et_simulate`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/boin_et_simulate.md)

## Examples

``` r
# Toxicity grades 0 or 1, 2, 3 and 4, and stable or progressive disease,
# partial response and complete response, at four doses
p_grade <- rbind(
  c(0.83, 0.62, 0.51, 0.34),
  c(0.12, 0.18, 0.19, 0.16),
  c(0.04, 0.11, 0.14, 0.15),
  c(0.01, 0.09, 0.16, 0.35)
)
p_resp_grade <- rbind(
  c(0.90, 0.50, 0.30, 0.25),
  c(0.05, 0.35, 0.20, 0.20),
  c(0.05, 0.15, 0.50, 0.55)
)
trials <- gboin_et_simulate(
  target = 0.47 / 1.5, target_resp = 0.95 / 1.5,
  p_grade = p_grade, p_resp_grade = p_resp_grade,
  score = c(0, 0.5, 1, 1.5), score_resp = c(0, 1, 1.5),
  n_cohort = 10, cohort_size = 3, n_trials = 200, seed = 123
)
trials
#> Simulated gBOIN-ET trials
#>   trials               : 200
#>   doses                : 4
#>   target toxicity score: 0.3133333
#>   target efficacy score: 0.6333333
#>   cohorts              : 10 of size 3
#>   max sample size      : 30
#> 
#> Average per trial:
#>                   DL1   DL2   DL3  DL4
#> Patients treated 4.41 10.66 13.24 1.35
#> Toxicity score   0.37  2.33  4.22 0.74
#> Efficacy score   0.43  3.99  8.38 0.94
#> 
#> Stopping reason (%):
#> 
#>              max_cohorts       no_admissible_dose no_lower_admissible_dose 
#>                       97                        1                        2 

# Patients of the first trial by dose and toxicity category
trials$n_grade[1, , ]
#>     G1 G2 G3 G4
#> DL1  3  0  0  0
#> DL2  2  0  1  0
#> DL3  7  6  0  2
#> DL4  5  0  0  4
```
