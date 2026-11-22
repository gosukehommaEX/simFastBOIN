# Operating Characteristics of the gBOIN Design in Several Scenarios

Run
[`sim_gboin`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_gboin.md)
for several sets of true grade probabilities with the same design, and
collect the results in one table.

## Usage

``` r
sim_gboin_multi(
  target,
  scenarios,
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

- scenarios:

  A list of scenarios, each a matrix of true grade probabilities as the
  argument `p_grade` of
  [`sim_gboin`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_gboin.md),
  or a list with the components `name` and `p_grade`. The names of the
  list are used as scenario names. Every scenario must have the same
  numbers of toxicity categories and of doses.

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

An object of class `boin_oc_multi`, which is a list with components

- results:

  The result of
  [`sim_gboin`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_gboin.md)
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
[`sim_boin_multi`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_boin_multi.md).

## References

Mu, R., Yuan, Y., Xu, J., Mandrekar, S. J. and Yin, J. (2019). gBOIN: a
unified model-assisted phase I trial design accounting for toxicity
grades, and binary or continuous end points. Journal of the Royal
Statistical Society: Series C, 68, 289-308.
[doi:10.1111/rssc.12263](https://doi.org/10.1111/rssc.12263)

## See also

[`sim_gboin`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_gboin.md),
[`sim_boin_multi`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_boin_multi.md)

## Examples

``` r
# Scenarios 8 and 10 of Mu et al. (2019) have the same mean scores but
# different grades
scenarios <- list(
  "Grade 3 only" = rbind(
    c(0.92, 0.76, 0.68, 0.57, 0.45, 0.25),
    c(0.00, 0.00, 0.00, 0.00, 0.00, 0.00),
    c(0.08, 0.24, 0.32, 0.43, 0.55, 0.75),
    c(0.00, 0.00, 0.00, 0.00, 0.00, 0.00)
  ),
  "Grade 2, then grade 3" = rbind(
    c(0.84, 0.52, 0.36, 0.14, 0.45, 0.25),
    c(0.16, 0.48, 0.64, 0.86, 0.00, 0.00),
    c(0.00, 0.00, 0.00, 0.00, 0.55, 0.75),
    c(0.00, 0.00, 0.00, 0.00, 0.00, 0.00)
  )
)
res <- sim_gboin_multi(
  target = 0.47 / 1.5,
  scenarios = scenarios,
  score = c(0, 0.5, 1, 1.5),
  n_cohort = 10,
  cohort_size = 3,
  n_earlystop = 100,
  n_trials = 300,
  seed = 123
)
res
#> gBOIN operating characteristics across 2 scenarios
#>   target score    : 0.3133333
#>   trials each     : 300
#>   max sample size : 30
#> 
#>               Scenario             Item DL1  DL2  DL3  DL4  DL5  DL6
#>           Grade 3 only   True score (%) 5.3 16.0 21.3 28.7 36.7 50.0
#>                        MTD selected (%) 0.3  5.3 19.0 47.0 26.7  1.7
#>                        Patients treated 3.7  5.5  7.5  7.7  4.3  1.3
#>                          Toxicity score 0.2  0.9  1.6  2.2  1.6  0.7
#>  Grade 2, then grade 3   True score (%) 5.3 16.0 21.3 28.7 36.7 50.0
#>                        MTD selected (%) 0.0  0.3  6.0 73.7 19.3  0.7
#>                        Patients treated 3.0  3.8  5.5 13.9  2.8  0.9
#>                          Toxicity score 0.2  0.6  1.2  4.0  1.0  0.5
#>  Total / No MTD
#>                
#>             0.0
#>            30.0
#>             7.1
#>                
#>             0.0
#>            30.0
#>             7.5
```
