# Operating Characteristics of the BOIN-ETC Design

Simulate trials of the BOIN-ETC design of Kakizume et al. (2024) for two
drugs given together with
[`boin_etc_simulate`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/boin_etc_simulate.md),
select the optimal dose combination (ODC) of every trial with
[`boin_etc_select_obd`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/boin_etc_select_obd.md),
and summarize the selection percentages and the numbers of patients,
DLTs and responses at every combination.

## Usage

``` r
sim_boin_etc(
  target,
  target_resp,
  p_true,
  p_resp,
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
  [`sim_boin_et`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_boin_et.md).

- p_true:

  Numeric matrix. True DLT probability of every combination, with one
  row per level of drug A and one column per level of drug B. At least
  two levels of each drug.

- p_resp:

  Numeric matrix of the same size. True response probability of every
  combination.

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
  [`boin_et_boundary`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/boin_et_boundary.md).
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
  [`sim_boin_et`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_boin_et.md).

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

An object of class `boin_etc_oc`, which is a list with components

- design:

  `"BOIN-ETC1"`, `"BOIN-ETC2"` or `"BOIN-ETC3"`.

- target, target_resp, p_true, p_resp:

  The targets and the true probabilities.

- n_rows, n_cols, n_trials:

  Numbers of levels of drug A and drug B, and of trials.

- lambda1, lambda2, eta1:

  The BOIN-ET boundaries.

- sel_percent:

  Matrix of the percentage of trials selecting each combination as the
  ODC.

- percent_no_obd, percent_early_stop:

  Percentages of trials selecting no ODC and stopped early.

- n_pts_dose, n_tox_dose, n_resp_dose:

  Matrices of the average numbers of patients, DLTs and responses at
  each combination.

- total_n_pts, total_n_tox, total_n_resp:

  Their totals.

- mean_subtrials:

  Average number of subtrials run.

- true_obd:

  Two-column matrix of the true optimal combinations: the combinations
  with the highest true utility among those with a DLT probability at
  most `target`.

- pct_true_obd_selected, pct_patients_true_obd:

  Percentage of trials selecting a true optimal combination (with
  `odc_tie = "all"`, the sum of `sel_percent` over them) and percentage
  of the patients treated there.

- overdose:

  List with the `cutoff` (`target`), the `combinations` above it, and
  the percentages of trials selecting one of them (`pct_selected`) and
  of patients treated there (`pct_patients`).

- stop_reason_percent:

  Percentage of trials by stopping reason.

- odc_tie, settings, call:

  The settings used.

- trials, obd, tied:

  With `keep_trials = TRUE`, the `boin_etc_trials` object, a two-column
  matrix of the selected combinations and the logical array of the tied
  candidates.

## Details

The percentages of patients are those of the average numbers of
patients. Kakizume et al. (2024) report the average over trials of the
percentage of the patients of a trial treated at the optimal
combinations ( from `trials`, as in
`inst/reproduce/reproduce_Kakizume2024.R`. Their optimal combinations
(bold in their Table 2) are not always those with the highest true
utility;
[`oc_compare`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/oc_compare.md)
accepts them through its argument `truth`.

## References

Kakizume, T., Takeda, K., Taguri, M. and Morita, S. (2024). BOIN-ETC: A
Bayesian optimal interval design considering efficacy and toxicity to
identify the optimal dose combinations. Statistical Methods in Medical
Research, 33(4), 716-727.
[doi:10.1177/09622802241236936](https://doi.org/10.1177/09622802241236936)

## See also

[`sim_boin_etc_multi`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_boin_etc_multi.md),
[`boin_etc_simulate`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/boin_etc_simulate.md),
[`boin_etc_select_obd`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/boin_etc_select_obd.md),
[`print.boin_etc_oc`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/print.boin_etc_oc.md)

## Examples

``` r
# Scenario 1 of Kakizume et al. (2024), BOIN-ETC1 with n_stop = 6
p_true <- matrix(c(0.05, 0.10, 0.10, 0.30, 0.30, 0.50), nrow = 2)
p_resp <- matrix(c(0.10, 0.25, 0.25, 0.60, 0.60, 0.65), nrow = 2)
oc <- sim_boin_etc(
  target = 0.30, target_resp = 0.60, p_true = p_true, p_resp = p_resp,
  n_cohort = 8, cohort_size = 3, n_trials = 200, seed = 123
)
oc
#> BOIN-ETC1 operating characteristics
#>   target DLT rate      : 0.3
#>   target response rate : 0.6
#>   dose levels          : 2 of drug A, 3 of drug B
#>   trials               : 200
#>   cohorts              : 8 of size 3
#>   max sample size      : 24
#>   subtrial end (n_stop): 6
#>   lambda1/lambda2/eta1 : 0.14 / 0.35 / 0.48
#>   start of a row       : diagonal
#>   rule 4               : paper
#>   utility scores       : 100, 40, 60, 0; ties lowest
#> 
#>                    Item Drug A   B1   B2   B3
#>       True DLT rate (%)     A1  5.0 10.0 30.0
#>                             A2 10.0 30.0 50.0
#>  True response rate (%)     A1 10.0 25.0 60.0
#>                             A2 25.0 60.0 65.0
#>        OBD selected (%)     A1  1.5 14.0 31.0
#>                             A2  9.5 41.0  3.0
#>        Patients treated     A1  3.2  4.2  4.4
#>                             A2  4.2  5.2  1.3
#>       Patients with DLT     A1  0.2  0.4  1.3
#>                             A2  0.4  1.5  0.7
#>  Patients with response     A1  0.3  1.0  2.7
#>                             A2  1.0  3.2  0.8
#> 
#> No OBD selected                      : 0%
#> Trials stopped early                 : 0%
#> Average number of patients           : 22.5
#> Average number of subtrials          : 2
#> 
#> True OBD (highest true utility)      : A2B2, A1B3
#>   Selected as the OBD                : 72%
#>   Patients treated there             : 42.6%
#> Combinations above the target DLT    : 1 of 6
#>   Selected as the OBD                : 3%
#>   Patients treated there             : 5.9%
```
