# Operating Characteristics of the BOIN-ETC Design Across Scenarios

Run
[`sim_boin_etc`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin_etc.md)
for several scenarios of DLT and response probabilities with the same
design and collect the results in one table.

## Usage

``` r
sim_boin_etc_multi(
  target,
  target_resp,
  scenarios,
  n_cohort,
  cohort_size,
  n_trials = 10000,
  subtrials = c("etc1", "etc2", "etc3"),
  n_earlystop = 6,
  etc1_start = c("diagonal", "below"),
  rule4 = c("paper", "authors_code"),
  p_saf = NULL,
  p_tox = NULL,
  p_saf_resp = NULL,
  prior = rep(1/6, 6),
  cutoff_eli = 0.95,
  p_eli = NULL,
  resp_cor = 0,
  resp_cor_type = c("copula", "binary"),
  psi00 = 40,
  psi11 = 60,
  utility_prior = c(1, 1),
  odc_tie = c("lowest", "all"),
  keep_trials = FALSE,
  verbose = FALSE,
  seed = 123
)
```

## Arguments

- target:

  Numeric scalar. Target DLT probability, phi.

- target_resp:

  Numeric scalar. Target response probability, delta, as in
  [`sim_boin_et`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin_et.md).

- scenarios:

  List of scenarios. Each element is a list with the matrices `p_true`
  and `p_resp` (one row per level of drug A, one column per level of
  drug B), and optionally `name`. Names of the list are used when `name`
  is absent. All scenarios must have the same numbers of levels.

- n_cohort:

  Integer scalar. Number of cohorts in a trial; with `cohort_size` it
  gives the maximum sample size.

- cohort_size:

  Integer scalar or vector. Number of patients per cohort. A scalar is
  used for every cohort; a shorter vector is padded with its last
  element.

- n_trials:

  Integer scalar. Number of trials to simulate. Defaults to 10000.

- subtrials:

  Character scalar. The division of the dose matrix into subtrials and
  their starting combinations: `"etc1"` (the default), `"etc2"` or
  `"etc3"`, BOIN-ETC1, 2 and 3 of the figures and tables of Kakizume et
  al. (2024). See Details.

- n_earlystop:

  Integer scalar. A subtrial ends when the combination chosen for the
  next cohort has this many patients (n_stop of Kakizume et al., 2024).
  Defaults to 6.

- etc1_start:

  Character scalar. Starting combination of a row subtrial of `"etc1"`:
  `"diagonal"` (the default) one level of drug B below the candidate of
  the previous subtrial, as the text and the figures of Kakizume et
  al. (2024) state; `"below"` the same level of drug B, the rule of the
  code of the authors that produced their published results. See
  Details.

- rule4:

  Character scalar. Choice among the neighboring combinations when the
  DLT rate is between `lambda1` and `lambda2` and the response rate at
  most `eta1`: `"paper"` (the default) the rule of the text;
  `"authors_code"` the rule of the code of the authors, needed to
  reproduce their results for BOIN-ETC1. See Details.

- p_saf, p_tox, p_saf_resp, prior:

  Design parameters of the BOIN-ET boundaries, as in
  [`boin_et_boundary`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_et_boundary.md).
  The defaults give `(lambda1, lambda2, eta1) = (0.14, 0.35, 0.48)` for
  a target DLT probability of 0.30 and a target response probability of
  0.60.

- cutoff_eli:

  Numeric scalar. A combination is excluded, with the higher
  combinations of its subtrial, when the posterior probability that its
  DLT rate exceeds `p_eli` is above this cutoff. Defaults to 0.95.

- p_eli:

  Numeric scalar or `NULL`. DLT probability of the exclusion rule.
  Defaults to `NULL`, `target`.

- resp_cor:

  Numeric scalar between -1 and 1. Correlation between the DLT and the
  response of a patient, see `resp_cor_type`. Defaults to 0, as in the
  simulations of Kakizume et al. (2024); the BOIN-ET functions default
  to 0.2.

- resp_cor_type:

  Character scalar. `"copula"` (the default) or `"binary"`, as in
  [`sim_boin_et`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin_et.md).

- psi00, psi11:

  Numeric scalars from 0 to 100. Utilities of a patient with neither a
  DLT nor a response and with both; a response without DLT scores 100
  and a DLT without response 0. Default to 40 and 60, the utilities
  (100, 40, 60, 0) of Kakizume et al. (2024).

- utility_prior:

  Numeric vector of two positive numbers. Beta prior of the utility of a
  combination divided by 100. Defaults to `c(1, 1)`.

- odc_tie:

  Character scalar. `"lowest"` (the default) counts one ODC per trial,
  the tied candidate with the lowest sum of the levels; `"all"` counts
  every tied candidate as selected, the definition of the published
  percentages of Kakizume et al. (2024), so that `sel_percent` can add
  up to more than 100.

- keep_trials:

  Logical scalar. Keep the simulated trials (`trials`) and the selection
  (`obd`, `tied`) in the result. Defaults to `FALSE`.

- verbose:

  Logical scalar. Print progress messages. Defaults to `FALSE`.

- seed:

  Integer scalar or `NULL`. Random seed. The state of the calling
  session is restored on exit. Defaults to 123.

## Value

An object of class `boin_etc_oc_multi`, which is a list with components

- results:

  Named list of `boin_etc_oc` objects, one per scenario.

- summary_table:

  Data frame with one row per scenario, as shown by
  [`print.boin_etc_oc_multi`](https://gosukehommaex.github.io/simFastBOIN/reference/print.boin_etc_oc_multi.md).

- scenario_names:

  Character vector of scenario names.

- n_rows, n_cols:

  Numbers of levels of drug A and drug B.

- call:

  The matched call.

## Details

Every scenario is simulated with the same `seed`, so that a scenario
simulated here matches the same scenario simulated on its own with
[`sim_boin_etc`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin_etc.md).

## See also

[`sim_boin_etc`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin_etc.md),
[`print.boin_etc_oc_multi`](https://gosukehommaex.github.io/simFastBOIN/reference/print.boin_etc_oc_multi.md)

## Examples

``` r
# Scenarios 1 and 2 of Kakizume et al. (2024)
scenarios <- list(
  "Scenario 1" = list(
    p_true = matrix(c(0.05, 0.10, 0.10, 0.30, 0.30, 0.50), nrow = 2),
    p_resp = matrix(c(0.10, 0.25, 0.25, 0.60, 0.60, 0.65), nrow = 2)),
  "Scenario 2" = list(
    p_true = matrix(c(0.05, 0.15, 0.15, 0.35, 0.30, 0.55), nrow = 2),
    p_resp = matrix(c(0.25, 0.60, 0.60, 0.65, 0.65, 0.70), nrow = 2))
)
res <- sim_boin_etc_multi(target = 0.30, target_resp = 0.60,
                          scenarios = scenarios, n_cohort = 8,
                          cohort_size = 3, n_trials = 200)
res
#> BOIN-ETC1 operating characteristics across 2 scenarios
#>   target DLT rate      : 0.3
#>   target response rate : 0.6
#>   dose levels          : 2 of drug A, 3 of drug B
#>   trials each          : 200
#>   max sample size      : 24
#>   subtrial end (n_stop): 6
#>   start of a row       : diagonal
#>   rule 4               : paper
#>   utility scores       : 100, 40, 60, 0; ties lowest
#> 
#>    Scenario   True OBD OBD selected (%) Patients at OBD (%) No OBD (%)
#>  Scenario 1 A2B2, A1B3               72                42.6          0
#>  Scenario 2 A2B1, A1B2               81                60.9          0
#>  Above limit selected (%) Patients above limit (%) Patients DLTs Responses
#>                       3.0                      5.9     22.5  4.5       9.0
#>                       7.5                      9.0     17.5  2.9       9.3
```
