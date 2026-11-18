# Operating Characteristics of the Comb-BOIN12 Design Across Scenarios

Run
[`sim_comb_boin12`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_comb_boin12.md)
for several scenarios of DLT and response probabilities with the same
design and collect the results in one table.

## Usage

``` r
sim_comb_boin12_multi(
  target,
  target_resp,
  scenarios,
  n_cohort,
  cohort_size,
  n_trials = 10000,
  start_dose = c(1, 1),
  n_earlystop = NULL,
  n_earlystop_rule = c("with_stay", "simple"),
  p_saf = NULL,
  p_tox = NULL,
  cutoff_eli = 0.95,
  cutoff_eli_resp = 0.9,
  extrasafe = FALSE,
  offset = 0.05,
  psi00 = 40,
  psi11 = 60,
  benchmark_weight = 0.5,
  utility_prior = c(1, 1),
  n_star = 6,
  tie_break = c("higher", "random"),
  no_candidate = c("nearest", "stop"),
  resp_cor = 0,
  resp_cor_type = c("copula", "binary"),
  est_method = c("logistic", "observed"),
  obd_bound = c("rectangle", "contour"),
  dose_a = NULL,
  dose_b = NULL,
  logistic_scale = c(2.5, 2.5),
  n_draws = 4000,
  keep_trials = FALSE,
  verbose = FALSE,
  seed = 123
)
```

## Arguments

- target:

  Numeric scalar. Highest acceptable DLT probability, phi_T. It is the
  target of the BOIN boundaries `lambda_e` and `lambda_d`, of the
  admissibility of a combination and of the MTD combination that bounds
  the OBD combination.

- target_resp:

  Numeric scalar. Lowest acceptable response probability, phi_E, as in
  [`sim_boin12`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin12.md).

- scenarios:

  List of scenarios. Each element is a list with the matrices `p_true`
  and `p_resp` (one row per level of drug A, one column per level of
  drug B), and optionally `name`. Names of the list are used when `name`
  is absent. All scenarios must have the same numbers of levels.

- n_cohort:

  Integer scalar. Number of cohorts in a trial.

- cohort_size:

  Integer scalar or vector. Number of patients per cohort. A scalar is
  used for every cohort; a shorter vector is padded with its last
  element.

- n_trials:

  Integer scalar. Number of trials to simulate. Defaults to 10000.

- start_dose:

  Integer vector of length 2. Levels of drug A and drug B of the first
  cohort. Defaults to `c(1, 1)`.

- n_earlystop:

  Integer scalar or `NULL`. Number of patients at the current
  combination that ends the trial when the next cohort would stay there
  (see `n_earlystop_rule`). Defaults to `NULL`, which runs every trial
  to the maximum sample size, as in Lu et al. (2025).

- n_earlystop_rule:

  Character scalar. `"with_stay"` (the default) ends the trial at
  `n_earlystop` only when the next cohort would stay at the current
  combination; `"simple"` ends it in any case.

- p_saf, p_tox:

  Numeric scalars or `NULL`. Highest DLT probability deemed
  subtherapeutic and lowest deemed overly toxic, which give the BOIN
  boundaries. Default to `0.6 * target` and `1.4 * target`.

- cutoff_eli:

  Numeric scalar. A combination is inadmissible, together with every
  combination with at least the same levels of both drugs, when the
  posterior probability that its DLT rate exceeds `target` is at least
  this cutoff. Defaults to 0.95.

- cutoff_eli_resp:

  Numeric scalar. A combination is inadmissible when the posterior
  probability that its response rate is below `target_resp` is at least
  this cutoff. Defaults to 0.90.

- extrasafe:

  Logical scalar. Apply the cutoff `cutoff_eli - offset` to the lowest
  combination, as Lu et al. (2025) do with 0.85. Defaults to `FALSE`.

- offset:

  Numeric scalar. Amount by which `cutoff_eli` is lowered at the lowest
  combination when `extrasafe` is `TRUE`. Defaults to 0.05, as in
  [`sim_comb_boin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_comb_boin.md);
  Lu et al. (2025) use 0.10.

- psi00, psi11:

  Numeric scalars from 0 to 100. Utility of a patient without DLT and
  without response, and with both, as in
  [`sim_boin12`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin12.md).
  Default to 40 and 60. The utility `pE - wT * pT` of Lu et al. (2025)
  is `psi00 = 100 * wT / (1 + wT)` and `psi11 = 100 / (1 + wT)`; their
  `wT = 0.11` is `psi00 = 9.90991` and `psi11 = 90.09009`.

- benchmark_weight:

  Numeric scalar above 0 and at most 1. Weight of the lowest acceptable
  utility in the utility benchmark, as in
  [`sim_boin12`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin12.md);
  it is the weight `w_u` of Lu et al. (2025), who use 0.6. Defaults to
  0.5.

- utility_prior:

  Numeric vector of two positive numbers. Beta prior of the utility of a
  combination divided by 100. Defaults to `c(1, 1)`.

- n_star:

  Integer scalar. With at least this many patients at the current
  combination and a DLT rate between `lambda_e` and `lambda_d`, the next
  combination is chosen from the current one and the two below it.
  Defaults to 6.

- tie_break:

  Character scalar. Choice among the combinations of equal desirability:
  `"higher"` (the default) keeps those with the highest sum of the two
  levels and draws one of them at random when several remain; `"random"`
  draws one of all the tied combinations. See Details.

- no_candidate:

  Character scalar. What happens when none of the combinations allowed
  for the next cohort is admissible: `"nearest"` (the default) goes to
  the admissible combination nearest to the current one; `"stop"` ends
  the trial without an OBD combination. See Details.

- resp_cor:

  Numeric scalar between -1 and 1. Correlation between the DLT and the
  response of a patient, see `resp_cor_type`. Defaults to 0.

- resp_cor_type:

  Character scalar. `"copula"` (the default) or `"binary"`, as in
  [`sim_boin12`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin12.md).

- est_method:

  Character scalar. Estimate of the response probability: `"logistic"`
  (the default) the posterior mean under the logistic model of Lu et al.
  (2025), see Details; `"observed"` the observed response rate.

- obd_bound:

  Character scalar. Combinations from which the OBDC is chosen:
  `"rectangle"` (the default) those with both levels at most those of
  the MTDC, as Lu et al. (2025) write; `"contour"` those, in every
  column (level of drug B), at most the combination of that column whose
  isotonic DLT rate is closest to `target`, the rule that reproduces the
  results published by Lu et al. (2025). See Details.

- dose_a, dose_b:

  Numeric vectors or `NULL`. Dose values of the levels of drug A and of
  drug B in the logistic model. Default to `NULL`, the level numbers.
  Each is standardized to mean 0 and standard deviation 0.5.

- logistic_scale:

  Numeric vector of two positive numbers. Scales of the Cauchy priors of
  the intercept and of the other coefficients of the logistic model.
  Defaults to `c(2.5, 2.5)`, as in Lu et al. (2025).

- n_draws:

  Integer scalar. Number of points of the importance sampling of the
  logistic model. Defaults to 4000.

- keep_trials:

  Logical scalar. Keep the simulated trials (`trials`) and the selected
  combinations (`obd`) in the result. Defaults to `FALSE`.

- verbose:

  Logical scalar. Print progress messages. Defaults to `FALSE`.

- seed:

  Integer scalar or `NULL`. Random seed. The state of the calling
  session is restored on exit. Defaults to 123.

## Value

An object of class `comb_boin12_oc_multi`, which is a list with
components

- results:

  Named list of `comb_boin12_oc` objects, one per scenario.

- summary_table:

  Data frame with one row per scenario, as shown by
  [`print.comb_boin12_oc_multi`](https://gosukehommaex.github.io/simFastBOIN/reference/print.comb_boin12_oc_multi.md).

- scenario_names:

  Character vector of scenario names.

- n_rows, n_cols:

  Numbers of levels of drug A and drug B.

- call:

  The matched call.

## Details

Every scenario is simulated with the same `seed`, so that a scenario
simulated here matches the same scenario simulated on its own with
[`sim_comb_boin12`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_comb_boin12.md).

## See also

[`sim_comb_boin12`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_comb_boin12.md)

## Examples

``` r
# Scenarios 16 and 17 of Lu et al. (2025)
scenarios <- list(
  "Scenario 16" = list(
    p_true = matrix(c(0.05, 0.10, 0.15, 0.15, 0.20, 0.45), nrow = 2),
    p_resp = matrix(c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60), nrow = 2)),
  "Scenario 17" = list(
    p_true = matrix(c(0.05, 0.10, 0.15, 0.20, 0.40, 0.50), nrow = 2),
    p_resp = matrix(c(0.10, 0.20, 0.30, 0.40, 0.50, 0.55), nrow = 2))
)
oc <- sim_comb_boin12_multi(
  target = 0.33, target_resp = 0.2, scenarios = scenarios,
  n_cohort = 14, cohort_size = 3, n_trials = 100,
  psi00 = 100 * 0.11 / 1.11, psi11 = 100 / 1.11, benchmark_weight = 0.6,
  extrasafe = TRUE, offset = 0.10, seed = 1
)
oc
#> Comb-BOIN12 operating characteristics across 2 scenarios
#>   DLT limit            : 0.33
#>   response limit       : 0.2
#>   dose levels          : 2 of drug A, 3 of drug B
#>   trials each          : 100
#>   max sample size      : 42
#>   ties / no candidate  : higher / nearest
#>   OBD selection        : logistic response, rectangle bound
#> 
#>     Scenario True OBD OBD selected (%) Patients at OBD (%) No OBD (%)
#>  Scenario 16     A1B3               40                31.3          2
#>  Scenario 17     A2B2               28                20.9          1
#>  Above limit selected (%) Patients above limit (%) Patients DLTs Responses
#>                        16                     22.8     41.2  9.4      17.8
#>                        40                     38.0     41.6 11.0      15.7
```
