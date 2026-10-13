# Operating Characteristics of the BOIN Drug-Combination Design Across Scenarios

Run
[`sim_comb_boin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_comb_boin.md)
under a list of matrices of true DLT probabilities that share the same
design, and collect the main results into one table.

## Usage

``` r
sim_comb_boin_multi(
  target,
  scenarios,
  n_cohort,
  cohort_size,
  n_trials = 10000,
  start_dose = c(1, 1),
  n_earlystop = 18,
  p_saf = NULL,
  p_tox = NULL,
  cutoff_eli = 0.95,
  extrasafe = FALSE,
  offset = 0.05,
  titration = FALSE,
  bound_mtd = FALSE,
  overdose_cutoff = NULL,
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
  `p_true` matrix, or a plain matrix of true DLT probabilities in which
  case the list names are used as scenario names. Every matrix must have
  the same numbers of rows (levels of drug A) and columns (levels of
  drug B).

- n_cohort:

  Integer scalar. Number of cohorts in a trial.

- cohort_size:

  Integer scalar or vector. Number of patients per cohort.

- n_trials:

  Integer scalar. Number of trials per scenario. Defaults to 10000.

- start_dose:

  Integer vector of length two. Dose levels of drug A and drug B for the
  first cohort. Defaults to `c(1, 1)`.

- n_earlystop:

  Integer scalar. Early stopping sample size at the current combination.
  Defaults to 18.

- p_saf:

  Numeric scalar. Highest DLT probability deemed subtherapeutic.
  Defaults to `0.6 * target`.

- p_tox:

  Numeric scalar. Lowest DLT probability deemed overly toxic. Defaults
  to `1.4 * target`.

- cutoff_eli:

  Numeric scalar. Posterior probability cutoff for eliminating a
  combination. Defaults to 0.95.

- extrasafe:

  Logical scalar. Apply the stricter safety stopping rule at the lowest
  combination.

- offset:

  Numeric scalar between 0 and 0.5. Relaxation of `cutoff_eli` for the
  safety stopping rule. Defaults to 0.05.

- titration:

  Logical scalar. Treat one patient per combination until the first DLT.

- bound_mtd:

  Logical scalar. Require the isotonic estimate at the selected
  combination to be below the de-escalation boundary.

- overdose_cutoff:

  Numeric scalar or `NULL`. Combinations whose true DLT probability
  exceeds this value count as overdoses. Defaults to `NULL`, which uses
  `target`.

- keep_trials:

  Logical scalar. Keep the trial by trial data of every scenario.

- verbose:

  Logical scalar. Report progress while the simulations run.

- seed:

  Integer scalar or `NULL`. Random seed, applied to every scenario, so
  that a scenario simulated here matches the same scenario simulated on
  its own with
  [`sim_comb_boin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_comb_boin.md).

## Value

An object of class `comb_boin_oc_multi`, which is a list with components

- results:

  Named list of `comb_boin_oc` objects, one per scenario.

- summary_table:

  Data frame with one row per scenario: the true MTD combinations (those
  closest to `target`), the percentage of trials selecting one of them,
  the percentage of patients treated at them, the percentage of trials
  without an MTD, the percentage of trials selecting a combination above
  `overdose_cutoff`, the percentage of patients treated above it, and
  the average numbers of patients and DLTs.

- scenario_names:

  Character vector of scenario names.

- n_rows:

  Number of dose levels of drug A.

- n_cols:

  Number of dose levels of drug B.

- call:

  The matched call.

## See also

[`sim_comb_boin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_comb_boin.md)

## Examples

``` r
scenarios <- list(
  "Two MTD combinations" = matrix(c(0.05, 0.10, 0.20,
                                    0.10, 0.20, 0.30,
                                    0.20, 0.30, 0.45), 3, byrow = TRUE),
  "All toxic" = matrix(c(0.40, 0.50, 0.60,
                         0.50, 0.60, 0.70,
                         0.60, 0.70, 0.80), 3, byrow = TRUE)
)

oc <- sim_comb_boin_multi(
  target = 0.30,
  scenarios = scenarios,
  n_cohort = 12,
  cohort_size = 3,
  n_trials = 200,
  seed = 123
)
oc
#> BOIN drug-combination operating characteristics across 2 scenarios
#>   target DLT rate : 0.3
#>   dose levels     : 3 of drug A, 3 of drug B
#>   trials each     : 200
#>   max sample size : 36
#>   overdose cutoff : 0.3
#> 
#>              Scenario   True MTD MTD selected (%) Patients at MTD (%)
#>  Two MTD combinations A3B2, A2B3             57.0                33.1
#>             All toxic       A1B1             60.5                76.9
#>  No MTD (%) Above cutoff selected (%) Patients above cutoff (%) Patients DLTs
#>         0.0                      11.5                      11.9     34.7  8.1
#>        34.5                      65.5                     100.0     19.0  8.0

# The full results of one scenario
oc$results[["All toxic"]]
#> BOIN drug-combination operating characteristics
#>   target DLT rate : 0.3
#>   dose levels     : 3 of drug A, 3 of drug B
#>   trials          : 200
#>   cohorts         : 12 of size 3
#>   max sample size : 36
#>   lambda_e / _d   : 0.2365 / 0.3585
#> 
#>               Item Drug A   B1   B2   B3
#>  True DLT rate (%)     A1 40.0 50.0 60.0
#>                        A2 50.0 60.0 70.0
#>                        A3 60.0 70.0 80.0
#>   MTD selected (%)     A1 60.5  2.0  0.0
#>                        A2  3.0  0.0  0.0
#>                        A3  0.0  0.0  0.0
#>   Patients treated     A1 14.6  2.2  0.1
#>                        A2  1.8  0.2  0.0
#>                        A3  0.1  0.0  0.0
#>  Patients with DLT     A1  5.7  1.1  0.0
#>                        A2  0.9  0.1  0.0
#>                        A3  0.1  0.0  0.0
#> 
#> No MTD selected                      : 34.5%
#> Average number of patients           : 19
#> Average number of DLTs               : 8
#> 
#> True MTD (closest to the target)     : A1B1
#>   Selected as the MTD                : 60.5%
#>   Patients treated there             : 76.9%
#> 
#> Combinations above a true DLT rate of 0.3: 9 of 9
#>   Patients treated there             : 100%
#>   Trials treating any patient there  : 100%
#>   Trials treating over 60% there     : 31.5%
#>   Trials treating over 80% there     : 17.5%
#>   Trials selecting an MTD there      : 65.5%
```
