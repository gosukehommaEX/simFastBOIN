# Operating Characteristics of the iBOIN Design in Several Scenarios

Run
[`sim_iboin`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_iboin.md)
for several true dose-toxicity curves with the same design and the same
prior information, and collect the results in one table.

## Usage

``` r
sim_iboin_multi(
  target,
  scenarios,
  skeleton,
  prior_ess,
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
  robust = FALSE,
  select_with_prior = TRUE,
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

  Numeric scalar. Target DLT probability, for example 0.30.

- scenarios:

  A list of scenarios, each a numeric vector of true DLT probabilities
  or a list with the components `name` and `p_true`, as for
  [`sim_boin_multi`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_boin_multi.md).
  Every scenario must have the same number of doses.

- skeleton:

  Numeric vector, one value per dose. Prior DLT probability at each
  dose, strictly between 0 and 1.

- prior_ess:

  Whole number or vector of whole numbers, one per dose. Prior effective
  sample size. See
  [`iboin_boundary`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/iboin_boundary.md).

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
  to 18. The simulations of the article do not stop early; a value above
  the maximum sample size switches the rule off.

- p_saf:

  Numeric scalar. Highest DLT probability deemed subtherapeutic.
  Defaults to `0.6 * target`.

- p_tox:

  Numeric scalar. Lowest DLT probability deemed overly toxic. Defaults
  to `1.4 * target`.

- cutoff_eli:

  Numeric scalar. Posterior probability cutoff for dose elimination.
  Defaults to 0.95.

- extrasafe:

  Logical scalar. Apply the stricter safety stopping rule at the lowest
  dose. Defaults to `FALSE`.

- offset:

  Numeric scalar between 0 and 0.5. Amount by which `cutoff_eli` is
  relaxed for the safety stopping rule. Defaults to 0.05.

- robust:

  Logical scalar. Use the robust prior, which drops the prior
  information above the prior MTD when the prior MTD lies in the upper
  half of the doses. Defaults to `FALSE`. See
  [`iboin_boundary`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/iboin_boundary.md).

- select_with_prior:

  Logical scalar. Add the prior information to the estimates of the DLT
  probabilities from which the MTD is selected, as
  [`iboin_select_mtd`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/iboin_select_mtd.md)
  does. Defaults to `TRUE`, which reproduces the operating
  characteristics of the article. With `FALSE` the MTD is selected from
  the data alone, by
  [`boin_select_mtd`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/boin_select_mtd.md).

- min_mtd_sample:

  Integer scalar. Smallest number of patients a dose must have been
  given to be selected as the MTD. Defaults to 1.

- overdose_cutoff:

  Numeric scalar or `NULL`. Doses with a true DLT probability above this
  value count as overdoses in the summary. Defaults to `NULL`, the
  target.

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
  [`sim_iboin`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_iboin.md)
  for every scenario.

- summary_table:

  Data frame with the true DLT rates, the selection percentages and the
  average numbers of patients and DLTs of every scenario.

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

Zhou, Y., Lee, J. J., Wang, S., Bailey, S. and Yuan, Y. (2021).
Incorporating historical information to improve phase I clinical trials.
Pharmaceutical Statistics, 20, 1017-1034.
[doi:10.1002/pst.2121](https://doi.org/10.1002/pst.2121)

## See also

[`sim_iboin`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_iboin.md),
[`sim_boin_multi`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_boin_multi.md)

## Examples

``` r
# The prior MTD is dose 3; it is the true MTD in the first scenario only
res <- sim_iboin_multi(
  target = 0.30,
  scenarios = list(
    "MTD at dose 3" = c(0.08, 0.15, 0.31, 0.45, 0.55),
    "MTD at dose 4" = c(0.09, 0.12, 0.15, 0.30, 0.45)
  ),
  skeleton = c(0.10, 0.19, 0.30, 0.42, 0.54),
  prior_ess = 3,
  n_cohort = 10,
  cohort_size = 3,
  n_earlystop = 100,
  n_trials = 300,
  seed = 123
)
res
#> iBOIN operating characteristics across 2 scenarios
#>   target DLT rate : 0.3
#>   trials each     : 300
#>   max sample size : 30
#> 
#>       Scenario              Item DL1  DL2  DL3  DL4  DL5 Total / No MTD
#>  MTD at dose 3 True DLT rate (%) 8.0 15.0 31.0 45.0 55.0               
#>                 MTD selected (%) 2.7 22.3 62.0 12.3  0.3            0.3
#>                 Patients treated 4.0  9.2 12.9  3.5  0.3           29.9
#>                Patients with DLT 0.4  1.4  4.0  1.6  0.2            7.7
#>  MTD at dose 4 True DLT rate (%) 9.0 12.0 15.0 30.0 45.0               
#>                 MTD selected (%) 1.0  3.0 30.3 59.7  5.7            0.3
#>                 Patients treated 3.7  5.2 10.6  8.2  2.3           29.9
#>                Patients with DLT 0.4  0.6  1.6  2.5  1.1            6.2
```
