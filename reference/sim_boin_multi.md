# Operating Characteristics Across Several Scenarios

Run
[`sim_boin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin.md)
under a list of dose-toxicity scenarios that share the same design, and
collect the results into one table.

## Usage

``` r
sim_boin_multi(
  target,
  scenarios,
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
  titration = FALSE,
  stay_on_1_of_3 = FALSE,
  bound_mtd = FALSE,
  mtd_max_estimate = NULL,
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

  Numeric scalar. Target DLT probability.

- scenarios:

  A list of scenarios. Each element is a list with a `name` and a
  `p_true` vector, or a plain numeric vector of true DLT probabilities
  in which case the list names are used as scenario names. Every
  scenario must have the same number of doses.

- n_cohort:

  Integer scalar. Number of cohorts in a trial.

- cohort_size:

  Integer scalar or vector. Number of patients per cohort.

- n_trials:

  Integer scalar. Number of trials per scenario. Defaults to 10000.

- start_dose:

  Integer scalar. Dose level for the first cohort. Defaults to 1. It is
  ignored when `titration` is `TRUE`.

- n_earlystop:

  Integer scalar. Early stopping sample size at the current dose.
  Defaults to 18.

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
  dose.

- offset:

  Numeric scalar between 0 and 0.5. Relaxation of `cutoff_eli` for the
  safety stopping rule. Defaults to 0.05.

- titration:

  Logical scalar. Start with single patient cohorts until the first DLT.

- stay_on_1_of_3:

  Logical scalar. Make one DLT out of three patients a stay rather than
  a de-escalation. See
  [`boin_boundary`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_boundary.md).

- bound_mtd:

  Logical scalar. Bound the isotonic estimate at the selected MTD by the
  de-escalation boundary.

- mtd_max_estimate:

  Numeric scalar or `NULL`. Largest isotonic estimate a dose may have
  and still be selected as the MTD. See
  [`sim_boin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin.md).

- min_mtd_sample:

  Integer scalar. Smallest number of patients for a dose to be eligible.

- overdose_cutoff:

  Numeric scalar or `NULL`. Doses whose true DLT probability exceeds
  this value count as overdoses. Defaults to `NULL`, which uses
  `target`.

- n_earlystop_rule:

  Character scalar, either `"with_stay"` or `"simple"`.

- keep_trials:

  Logical scalar. Keep the trial by trial data of every scenario.

- verbose:

  Logical scalar. Report progress while the simulations run.

- seed:

  Integer scalar or `NULL`. Random seed, applied to every scenario, so
  that a scenario simulated here matches the same scenario simulated on
  its own with
  [`sim_boin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin.md).

## Value

An object of class `boin_oc_multi`, which is a list with components

- results:

  Named list of `boin_oc` objects, one per scenario.

- summary_table:

  Data frame collecting all scenarios.

- scenario_names:

  Character vector of scenario names.

- n_doses:

  Number of dose levels.

- call:

  The matched call.

## See also

[`sim_boin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin.md)

## Examples

``` r
scenarios <- list(
  list(name = "MTD at dose 3", p_true = c(0.05, 0.15, 0.30, 0.45, 0.60)),
  list(name = "All toxic",     p_true = c(0.35, 0.45, 0.55, 0.65, 0.75))
)

oc <- sim_boin_multi(
  target = 0.30,
  scenarios = scenarios,
  n_cohort = 10,
  cohort_size = 3,
  n_trials = 200,
  seed = 123
)
oc
#> BOIN operating characteristics across 2 scenarios
#>   target DLT rate : 0.3
#>   trials each     : 200
#>   max sample size : 30
#> 
#>       Scenario              Item  DL1  DL2  DL3  DL4  DL5 Total / No MTD
#>  MTD at dose 3 True DLT rate (%)  5.0 15.0 30.0 45.0 60.0               
#>                 MTD selected (%)  1.0 26.0 52.0 21.0  0.0            0.0
#>                 Patients treated  4.2  9.1 10.6  4.5  0.6           29.1
#>                Patients with DLT  0.2  1.4  3.2  2.1  0.4            7.3
#>      All toxic True DLT rate (%) 35.0 45.0 55.0 65.0 75.0               
#>                 MTD selected (%) 67.0  9.0  2.0  0.0  0.0           22.0
#>                 Patients treated 14.5  5.4  0.8  0.1  0.0           20.8
#>                Patients with DLT  5.0  2.5  0.4  0.1  0.0            8.0
```
