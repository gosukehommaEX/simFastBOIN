# Operating Characteristics of the gBOIN Design

Simulate trials of the gBOIN design of Mu et al. (2019) for a toxicity
outcome with several grades that each carry a score, and summarize how
often each dose is selected as the MTD and how the patients are
distributed over the doses.

## Usage

``` r
sim_gboin(
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
  [`boin_boundary`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/boin_boundary.md),
  a whole number, as in the 'UnifiedDoseFinding' package. See the
  details.

- min_mtd_sample:

  Integer scalar. Smallest number of patients a dose must have been
  given to be selected as the MTD. Defaults to 1.

- overdose_cutoff:

  Numeric scalar or `NULL`. Doses with a true mean normalized score
  above this value count as overdoses in the summary. Defaults to
  `NULL`, the target.

- n_earlystop_rule:

  Character scalar, either `"with_stay"` or `"simple"`. Under
  `"with_stay"` the trial stops at `n_earlystop` only when the design
  would remain at the current dose. Under `"simple"` it stops as soon as
  the threshold is reached.

- keep_trials:

  Logical scalar. Keep the trial level data in the result. Defaults to
  `FALSE`.

- verbose:

  Logical scalar. Report the progress in messages. Defaults to `FALSE`.

- seed:

  Integer scalar or `NULL`. Random seed. When `NULL` the current state
  of the random number generator is used. The state of the calling
  session is restored on exit. Defaults to 123.

## Value

An object of class `c("gboin_oc", "boin_oc")`, a list with the
components of the result of
[`sim_boin`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_boin.md),
in which

- p_true:

  is the true mean normalized score of each dose,

- n_tox_dose, total_n_tox:

  are averages of the total normalized score, in the place of the number
  of DLTs,

and in addition

- design:

  The name of the design, `"gBOIN"`.

- p_grade, score:

  The true grade probabilities and the scores.

- n_grade_dose:

  Matrix, toxicity categories by doses, of the average number of
  patients in each category.

With `keep_trials = TRUE` the component `trials` also holds the array
`n_grade` of
[`gboin_simulate`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/gboin_simulate.md).

## Details

The trials are simulated by
[`gboin_simulate`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/gboin_simulate.md),
which describes the design, and the MTD is selected by
[`gboin_select_mtd`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/gboin_select_mtd.md).
A trial stopped because the lowest dose is eliminated selects no dose.

With two toxicity categories of scores 0 and 1 the design is the BOIN
design, and the result is, for the same seed, that of
[`sim_boin`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_boin.md)
without titration.

The result can be passed to
[`oc_results`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/oc_results.md)
and
[`oc_compare`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/oc_compare.md),
where it appears under the name gBOIN.

## References

Mu, R., Yuan, Y., Xu, J., Mandrekar, S. J. and Yin, J. (2019). gBOIN: a
unified model-assisted phase I trial design accounting for toxicity
grades, and binary or continuous end points. Journal of the Royal
Statistical Society: Series C, 68, 289-308.
[doi:10.1111/rssc.12263](https://doi.org/10.1111/rssc.12263)

## See also

[`gboin_simulate`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/gboin_simulate.md),
[`gboin_select_mtd`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/gboin_select_mtd.md),
[`sim_gboin_multi`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_gboin_multi.md),
[`sim_boin`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_boin.md)

## Examples

``` r
# Scenario 1 of Mu et al. (2019): grades 0 or 1, 2, 3 and 4 at six doses,
# equivalent toxicity scores 0, 0.5, 1 and 1.5, target score 0.47
p_grade <- rbind(
  c(0.83, 0.75, 0.62, 0.51, 0.34, 0.19),
  c(0.12, 0.15, 0.18, 0.19, 0.16, 0.11),
  c(0.04, 0.07, 0.11, 0.14, 0.15, 0.11),
  c(0.01, 0.03, 0.09, 0.16, 0.35, 0.59)
)
oc <- sim_gboin(
  target = 0.47 / 1.5,
  p_grade = p_grade,
  score = c(0, 0.5, 1, 1.5),
  n_cohort = 10,
  cohort_size = 3,
  n_earlystop = 100,
  n_trials = 500,
  seed = 123
)
oc
#> gBOIN operating characteristics
#>   target score    : 0.3133333
#>   trials          : 500
#>   cohorts         : 10 of size 3
#>   max sample size : 30
#>   lambda_e / _d   : 0.2471 / 0.3746
#> 
#>                  DL1  DL2  DL3  DL4  DL5  DL6 Total / No MTD
#> True score (%)   7.7 12.7 22.3 31.7 50.3 70.0               
#> MTD selected (%) 0.2  3.4 32.0 53.8 10.4  0.2            0.0
#> Patients treated 3.6  5.2  9.2  9.0  2.7  0.3           30.0
#> Toxicity score   0.3  0.6  2.0  2.9  1.4  0.2            7.5
#> 
#> Doses above a true score of 0.3133333: DL4, DL5, DL6
#>   Patients treated there              : 40.1%
#>   Trials treating any patient there   : 86.2%
#>   Trials treating over 60% there      : 23.4%
#>   Trials treating over 80% there      : 0%
#>   Trials selecting an MTD there       : 64.4%
```
