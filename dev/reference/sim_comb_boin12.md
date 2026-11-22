# Operating Characteristics of the Comb-BOIN12 Design

Simulate trials of the Comb-BOIN12 design of Lu et al. (2025) for two
drugs given together with
[`comb_boin12_simulate`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/comb_boin12_simulate.md),
select the optimal biological dose combination (OBDC) of every trial
with
[`comb_boin12_select_obd`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/comb_boin12_select_obd.md),
and summarize the selection percentages and the numbers of patients,
DLTs and responses at every combination.

## Usage

``` r
sim_comb_boin12(
  target,
  target_resp,
  p_true,
  p_resp,
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
  [`sim_boin12`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_boin12.md).

- p_true:

  Numeric matrix. True DLT probability of every combination, with one
  row per level of drug A and one column per level of drug B.

- p_resp:

  Numeric matrix of the same size. True response probability of every
  combination.

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
  [`sim_comb_boin`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_comb_boin.md);
  Lu et al. (2025) use 0.10.

- psi00, psi11:

  Numeric scalars from 0 to 100. Utility of a patient without DLT and
  without response, and with both, as in
  [`sim_boin12`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_boin12.md).
  Default to 40 and 60. The utility `pE - wT * pT` of Lu et al. (2025)
  is `psi00 = 100 * wT / (1 + wT)` and `psi11 = 100 / (1 + wT)`; their
  `wT = 0.11` is `psi00 = 9.90991` and `psi11 = 90.09009`.

- benchmark_weight:

  Numeric scalar above 0 and at most 1. Weight of the lowest acceptable
  utility in the utility benchmark, as in
  [`sim_boin12`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_boin12.md);
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
  [`sim_boin12`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_boin12.md).

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

An object of class `comb_boin12_oc`, which is a list with components

- design:

  `"Comb-BOIN12"`.

- target, target_resp, p_true, p_resp:

  The limits and the true probabilities.

- n_rows, n_cols, n_trials:

  Numbers of levels of drug A and drug B, and of trials.

- lambda_e, lambda_d, u_benchmark:

  The BOIN boundaries and the utility benchmark (scale 0 to 100).

- sel_percent:

  Matrix of the percentage of trials selecting each combination as the
  OBDC.

- percent_no_obd, percent_early_stop:

  Percentages of trials selecting no OBDC and stopped early.

- n_pts_dose, n_tox_dose, n_resp_dose:

  Matrices of the average numbers of patients, DLTs and responses at
  each combination.

- total_n_pts, total_n_tox, total_n_resp:

  Their totals.

- true_obd:

  Two-column matrix of the true OBDCs: the combinations with the highest
  true utility among those with a DLT probability at most `target`.

- pct_true_obd_selected, pct_patients_true_obd:

  Percentage of trials selecting a true OBDC and percentage of the
  patients treated there.

- overdose:

  List with the `cutoff` (`target`), the `combinations` above it, and
  the percentages of trials selecting one of them (`pct_selected`) and
  of patients treated there (`pct_patients`).

- stop_reason_percent:

  Percentage of trials by stopping reason.

- est_method, obd_bound, settings, call:

  The settings used.

- trials, obd:

  With `keep_trials = TRUE`, the `comb_boin12_trials` object and a
  two-column matrix of the selected combinations.

## Details

The percentages of patients are those of the average numbers of
patients. Lu et al. (2025) report the average over trials of the number
of patients as a percentage of the maximum sample size; with
`keep_trials = TRUE` it can be computed from `trials`, as in
`inst/reproduce/reproduce_Lu2025.R`.

With the logistic estimate, each trial takes a few milliseconds to
select; `est_method = "observed"` is much faster.

## References

Lu, M., Zhang, J., Yuan, Y. and Lin, R. (2025). Comb-BOIN12: A
utility-based Bayesian optimal interval design for dose optimization in
cancer drug-combination trials. Statistics in Biopharmaceutical
Research, 17(2), 266-276.
[doi:10.1080/19466315.2024.2370403](https://doi.org/10.1080/19466315.2024.2370403)

## See also

[`sim_comb_boin12_multi`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_comb_boin12_multi.md),
[`comb_boin12_simulate`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/comb_boin12_simulate.md),
[`comb_boin12_select_obd`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/comb_boin12_select_obd.md),
[`print.comb_boin12_oc`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/print.comb_boin12_oc.md)

## Examples

``` r
# Scenario 16 of Lu et al. (2025), with their design settings
p_true <- matrix(c(0.05, 0.10, 0.15, 0.15, 0.20, 0.45), nrow = 2)
p_resp <- matrix(c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60), nrow = 2)
oc <- sim_comb_boin12(
  target = 0.33, target_resp = 0.2, p_true = p_true, p_resp = p_resp,
  n_cohort = 14, cohort_size = 3, n_trials = 200,
  psi00 = 100 * 0.11 / 1.11, psi11 = 100 / 1.11, benchmark_weight = 0.6,
  extrasafe = TRUE, offset = 0.10, dose_a = c(100, 200),
  dose_b = c(50, 100, 200), seed = 123
)
oc
#> Comb-BOIN12 operating characteristics
#>   DLT limit            : 0.33
#>   response limit       : 0.2
#>   dose levels          : 2 of drug A, 3 of drug B
#>   trials               : 200
#>   cohorts              : 14 of size 3
#>   max sample size      : 42
#>   lambda_e / lambda_d  : 0.2604 / 0.3947
#>   utility scores       : 100, 9.9099, 90.0901, 0; benchmark 54.79
#>   ties / no candidate  : higher / nearest
#>   OBD selection        : logistic response, rectangle bound
#> 
#>                    Item Drug A   B1   B2   B3
#>       True DLT rate (%)     A1  5.0 15.0 20.0
#>                             A2 10.0 15.0 45.0
#>  True response rate (%)     A1 10.0 30.0 50.0
#>                             A2 20.0 40.0 60.0
#>        OBD selected (%)     A1  5.0 12.0 41.0
#>                             A2  5.5 18.5 17.0
#>        Patients treated     A1  4.2  5.6 12.1
#>                             A2  3.5  8.2  8.0
#>       Patients with DLT     A1  0.2  0.9  2.5
#>                             A2  0.3  1.3  3.6
#>  Patients with response     A1  0.4  1.7  6.2
#>                             A2  0.7  3.2  4.8
#> 
#> No OBD selected                      : 1%
#> Trials stopped early                 : 1%
#> Average number of patients           : 41.6
#> 
#> True OBD (highest true utility)      : A1B3
#>   Selected as the OBD                : 41%
#>   Patients treated there             : 29.1%
#> Combinations above the DLT limit     : 1 of 6
#>   Selected as the OBD                : 17%
#>   Patients treated there             : 19.2%
```
