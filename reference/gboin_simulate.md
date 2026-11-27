# Simulate gBOIN Trials

Run the gBOIN dose-finding algorithm of Mu et al. (2019) over many
simulated trials, for a toxicity outcome with several grades that each
carry a score, and return the patients, their grades and the total score
at every dose. This is the simulation engine behind
[`sim_gboin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_gboin.md),
exposed for users who want the raw trial data.

## Usage

``` r
gboin_simulate(
  target,
  p_grade,
  score,
  n_cohort,
  cohort_size,
  n_trials = 10000,
  start_dose = 1,
  n_earlystop = 18,
  p_saf = NULL,
  p_tox = NULL,
  cutoff_eli = 0.95,
  extrasafe = FALSE,
  offset = 0.05,
  elim_rule = c("posterior", "count"),
  n_earlystop_rule = c("with_stay", "simple"),
  seed = 123
)
```

## Arguments

- target:

  Numeric scalar. Target value of the mean normalized toxicity score,
  for example 0.313 for a target equivalent toxicity score of 0.47 with
  a largest score of 1.5. See the details.

- p_grade:

  Numeric matrix. True probability of each toxicity category (rows, from
  the lowest to the highest) at each dose (columns, in increasing dose
  order). Every column sums to 1.

- score:

  Numeric vector, one value per row of `p_grade`. Score of each toxicity
  category, for example `c(0, 0.5, 1, 1.5)` for grades 0 or 1, 2, 3
  and 4. The scores are divided by their maximum.

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

  Integer scalar. The trial stops once this many patients have been
  treated at the current dose and the design would stay there. Defaults
  to 18. A value above the maximum sample size switches the rule off.

- p_saf:

  Numeric scalar. Highest mean score deemed subtherapeutic. Defaults to
  `0.6 * target`.

- p_tox:

  Numeric scalar. Lowest mean score deemed overly toxic. Defaults to
  `1.4 * target`.

- cutoff_eli:

  Numeric scalar. Posterior probability cutoff for dose elimination.
  Defaults to 0.95.

- extrasafe:

  Logical scalar. Apply the stricter safety stopping rule at the lowest
  dose. Defaults to `FALSE`.

- offset:

  Numeric scalar between 0 and 0.5. Amount by which `cutoff_eli` is
  relaxed for the safety stopping rule. Defaults to 0.05.

- elim_rule:

  Character scalar, either `"posterior"` or `"count"`. How a dose is
  eliminated. Under `"posterior"`, the default and the rule of the
  article, the posterior probability is evaluated at the total score.
  Under `"count"` the dose is eliminated when the total score reaches
  the elimination boundary of the number of DLTs of
  [`boin_boundary`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_boundary.md),
  a whole number, as in the 'UnifiedDoseFinding' package. See the
  details.

- n_earlystop_rule:

  Character scalar, either `"with_stay"` or `"simple"`. Under
  `"with_stay"` the trial stops at `n_earlystop` only when the design
  would remain at the current dose. Under `"simple"` it stops as soon as
  the threshold is reached.

- seed:

  Integer scalar or `NULL`. Random seed. When `NULL` the current state
  of the random number generator is used. The state of the calling
  session is restored on exit. Defaults to 123.

## Value

An object of class `c("gboin_trials", "boin_trials")`, which is a list
with components

- n_pts:

  Integer matrix, `n_trials` by number of doses, of patients treated.

- n_tox:

  Numeric matrix of the same shape: the total of the normalized scores
  of the patients treated at each dose (the total of the scores as
  given, divided by the largest score), which takes the place of the
  number of DLTs.

- n_grade:

  Integer array, `n_trials` by doses by toxicity categories, of the
  patients in each category. The categories take the row names of
  `p_grade` when every row has one, and `G1`, `G2`, ... otherwise.

- eliminated:

  Logical matrix marking doses eliminated during the trial.

- cohorts_used:

  Integer vector of cohorts completed in each trial.

- stop_reason:

  Character vector giving why each trial ended.

- lambda_e, lambda_d:

  The escalation and de-escalation boundaries of the mean score.

- design:

  The name of the design, `"gBOIN"`.

- settings:

  List of the design parameters, with the true mean normalized score of
  each dose in `p_true`.

## Details

The gBOIN design generalizes the BOIN design to toxicity outcomes from
the exponential family. This function covers the quasi-Bernoulli
outcome: every toxicity category carries a score, such as the equivalent
toxicity score (ETS) of Yuan et al. (2007), and the scores divided by
their maximum, which lie between 0 and 1, are treated like fractional
DLTs. With `x` the total normalized score of the `n` patients treated at
the current dose, the design

- escalates when `x / n <= lambda_e`,

- de-escalates when `x / n >= lambda_d`,

- stays otherwise,

where `lambda_e` and `lambda_d` are the BOIN boundaries of
[`boin_lambda`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_lambda.md)
for `target`, `p_saf` and `p_tox`; they are the optimal boundaries for
the quasi-Bernoulli outcome as well (Theorem 1 of the article). A dose
with at least three patients is eliminated, together with every higher
dose, when the posterior probability that its mean score exceeds
`target`, from the distribution Beta(`x + 1`, `n - x + 1`), is above
`cutoff_eli`; the next cohort is then treated at the dose below, and the
trial stops when the lowest dose is eliminated.

`target`, `p_saf` and `p_tox` are on the normalized scale. The article
computes its boundaries for the ETS example (0.373 and 0.566) from the
unnormalized target of 0.47; on the normalized scale the target is 0.47
/ 1.5 = 0.313 and the boundaries, multiplied by 1.5, are 0.371 and
0.562.

The reference implementation in the 'UnifiedDoseFinding' package applies
the elimination boundary of the number of DLTs, a whole number, to the
total score, which eliminates a dose a little later than the posterior
probability evaluated at the total score itself, the rule the article
describes. `elim_rule = "count"` gives the rule of the package; it is
meant for comparisons with it. For whole numbers of DLTs the two rules
are the same.

One uniform variate is drawn per patient, in enrollment order: the
patient falls in the highest category for which the probability of that
category or a higher one exceeds the variate. A cohort cut short by the
maximum sample size still draws a full cohort's worth. With two
categories of scores 0 and 1 the trials are therefore, for the same
seed, those of
[`boin_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_simulate.md)
without titration.

The total normalized score of a dose is the total of the scores as
given, divided once by the largest score, as in
[`gboin_et_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/gboin_et_simulate.md).
With scores that are binary fractions, such as 0, 0.5, 1 and 1.5, the
total of the scores as given is exact, so that the same patients give
the same total, and the same decision, on every platform.

Possible values of `stop_reason` are `"lowest_dose_eliminated"`,
`"lowest_dose_too_toxic"`, `"n_earlystop"`, `"max_sample_size"` and
`"max_cohorts"`.

## References

Mu, R., Yuan, Y., Xu, J., Mandrekar, S. J. and Yin, J. (2019). gBOIN: a
unified model-assisted phase I trial design accounting for toxicity
grades, and binary or continuous end points. Journal of the Royal
Statistical Society: Series C, 68, 289-308.
[doi:10.1111/rssc.12263](https://doi.org/10.1111/rssc.12263)

Yuan, Z., Chappell, R. and Bailey, H. (2007). The continual reassessment
method for multiple toxicity grades: a Bayesian quasi-likelihood
approach. Biometrics, 63, 173-179.

## See also

[`sim_gboin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_gboin.md),
[`gboin_select_mtd`](https://gosukehommaex.github.io/simFastBOIN/reference/gboin_select_mtd.md),
[`boin_lambda`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_lambda.md)

## Examples

``` r
# Scenario 1 of Mu et al. (2019): grades 0 or 1, 2, 3 and 4 at six doses
p_grade <- rbind(
  c(0.83, 0.75, 0.62, 0.51, 0.34, 0.19),
  c(0.12, 0.15, 0.18, 0.19, 0.16, 0.11),
  c(0.04, 0.07, 0.11, 0.14, 0.15, 0.11),
  c(0.01, 0.03, 0.09, 0.16, 0.35, 0.59)
)
trials <- gboin_simulate(
  target = 0.47 / 1.5,
  p_grade = p_grade,
  score = c(0, 0.5, 1, 1.5),
  n_cohort = 10,
  cohort_size = 3,
  n_trials = 200,
  seed = 123
)
trials
#> Simulated gBOIN trials
#>   trials          : 200
#>   doses           : 6
#>   target score    : 0.3133333
#>   cohorts         : 10 of size 3
#>   max sample size : 30
#> 
#> Average per trial:
#>                   DL1  DL2  DL3  DL4  DL5  DL6
#> Patients treated 3.42 5.20 8.85 8.87 2.79 0.33
#> Toxicity score   0.24 0.68 1.95 2.83 1.39 0.26
#> 
#> Stopping reason (%):
#> 
#> max_sample_size     n_earlystop 
#>              87              13 

# Average number of patients of each grade at each dose
apply(trials$n_grade, c(3, 2), mean)
#>      DL1   DL2   DL3   DL4   DL5   DL6
#> G1 2.845 3.910 5.545 4.425 0.990 0.035
#> G2 0.430 0.740 1.565 1.830 0.395 0.035
#> G3 0.135 0.375 0.950 1.165 0.435 0.035
#> G4 0.010 0.180 0.790 1.445 0.970 0.225
```
