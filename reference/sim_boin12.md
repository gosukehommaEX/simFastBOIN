# Operating Characteristics of the BOIN12 Design

Simulate trials of the BOIN12 design of Lin et al. (2020) under one
scenario of DLT and response probabilities, select the optimal
biological dose (OBD) of each trial, and summarize how often each dose
is selected and how many patients, DLTs and responses occur at each
dose.

## Usage

``` r
sim_boin12(
  target,
  target_resp,
  p_true,
  p_resp,
  n_cohort,
  cohort_size,
  n_trials = 10000,
  start_dose = 1,
  n_earlystop = NULL,
  n_earlystop_rule = c("with_stay", "simple"),
  p_saf = NULL,
  p_tox = NULL,
  cutoff_eli = 0.95,
  cutoff_eli_resp = 0.9,
  psi00 = 40,
  psi11 = 60,
  benchmark_weight = 0.5,
  utility_prior = c(1, 1),
  n_star = 6,
  n_explore = 9,
  resp_cor = 0,
  resp_cor_type = c("copula", "binary"),
  est_method = c("posterior_mean", "observed", "tail_probability"),
  keep_trials = FALSE,
  verbose = FALSE,
  seed = 123
)
```

## Arguments

- target:

  Numeric scalar. Highest acceptable DLT probability, called phi_T by
  Lin et al. (2020). It is the target of the BOIN boundaries `lambda_e`
  and `lambda_d`, of the admissibility of a dose and of the MTD that
  bounds the OBD. Lin et al. recommend a value slightly above the target
  of a design that only considers toxicity, for example 0.35 instead of
  0.30.

- target_resp:

  Numeric scalar. Lowest acceptable response probability, phi_E. Note
  that this is a lower limit, as in
  [`sim_bf_boin_et`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_bf_boin_et.md),
  whereas the `target_resp` of
  [`sim_boin_et`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin_et.md)
  is a target.

- p_true:

  Numeric vector. True DLT probability at each dose.

- p_resp:

  Numeric vector. True response probability at each dose.

- n_cohort:

  Integer scalar. Number of cohorts in a trial.

- cohort_size:

  Integer scalar or vector. Number of patients per cohort. A scalar is
  used for every cohort; a shorter vector is padded with its last
  element.

- n_trials:

  Integer scalar. Number of trials to simulate. Defaults to 10000.

- start_dose:

  Integer scalar. Dose level of the first cohort. Defaults to 1.

- n_earlystop:

  Integer scalar or `NULL`. Number of patients at the current dose that
  ends the trial when the next cohort would stay there (see
  `n_earlystop_rule`). Defaults to `NULL`, which runs every trial to the
  maximum sample size, as in Lin et al. (2020).

- n_earlystop_rule:

  Character scalar. `"with_stay"` (the default) ends the trial at
  `n_earlystop` only when the next cohort would stay at the current
  dose, as in the example of Zhou et al. (2022); `"simple"` ends it in
  any case.

- p_saf, p_tox:

  Numeric scalars or `NULL`. Highest DLT probability deemed
  subtherapeutic and lowest deemed overly toxic, which give the BOIN
  boundaries. Default to `0.6 * target` and `1.4 * target`.

- cutoff_eli:

  Numeric scalar. A dose is inadmissible, together with every higher
  dose, when the posterior probability that its DLT rate exceeds
  `target` is at least this cutoff. Defaults to 0.95.

- cutoff_eli_resp:

  Numeric scalar. A dose is inadmissible when the posterior probability
  that its response rate is below `target_resp` is at least this cutoff.
  Defaults to 0.90.

- psi00, psi11:

  Numeric scalars from 0 to 100. Utility of a patient without DLT and
  without response, and with both. The utility of a response without DLT
  is 100 and that of a DLT without response is 0. Default to 40 and 60,
  the example of Lin et al. (2020).

- benchmark_weight:

  Numeric scalar above 0 and at most 1. The utility benchmark is
  `benchmark_weight * u_low + (1 - benchmark_weight) * 100`, where
  `u_low` is the expected utility of a dose whose DLT and response
  probabilities are `target` and `target_resp`, the two outcomes being
  independent. Defaults to 0.5, the benchmark of Lin et al. (2020); Lu
  et al. (2025) use 0.6.

- utility_prior:

  Numeric vector of two positive numbers. Beta prior of the utility of a
  dose divided by 100. Defaults to `c(1, 1)`.

- n_star:

  Integer scalar. With at least this many patients at the current dose
  and a DLT rate between `lambda_e` and `lambda_d`, the next dose is
  chosen from the current dose and the dose below. Defaults to 6.

- n_explore:

  Integer scalar. With at least this many patients at the current dose,
  a DLT rate below `lambda_d` and an untried admissible dose above, the
  next cohort goes to that dose. Defaults to 9.

- resp_cor:

  Numeric scalar between -1 and 1. Correlation between the DLT and the
  response of a patient, see `resp_cor_type`. Defaults to 0.

- resp_cor_type:

  Character scalar. `"copula"` (the default) links the two outcomes by a
  Gaussian copula with correlation `resp_cor`, the bivariate normal
  model of Lin et al. (2020); `"binary"` sets the correlation of the two
  binary outcomes, as in
  [`boin_et_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_et_simulate.md).

- est_method:

  Character scalar. Estimate of the utility of a dose:
  `"posterior_mean"` (the default), the posterior mean of the
  quasi-binomial model given in Section 4 of the Supplement of Lin et
  al. (2020); `"observed"`, the observed mean utility, which the example
  of Section 6 of that Supplement reports; or `"tail_probability"`, the
  desirability used during the trial, which the 'escalation' package
  uses.

- keep_trials:

  Logical scalar. Keep the simulated trials and the selected OBDs in the
  result. Defaults to `FALSE`.

- verbose:

  Logical scalar. Report progress with messages. Defaults to `FALSE`.

- seed:

  Integer scalar or `NULL`. Random seed. The state of the calling
  session is restored on exit. Defaults to 123.

## Value

An object of class `boin12_oc`, which is a list with components

- design:

  `"BOIN12"`.

- target, target_resp, p_true, p_resp:

  The limits and the scenario.

- n_doses, n_trials:

  Numbers of doses and of simulated trials.

- lambda_e, lambda_d, u_benchmark:

  The BOIN boundaries and the utility benchmark.

- sel_percent:

  Percentage of trials selecting each dose as the OBD.

- percent_no_obd:

  Percentage of trials selecting no OBD.

- percent_early_stop:

  Percentage of trials stopped because no dose was admissible.

- n_pts_dose, n_tox_dose, n_resp_dose:

  Average numbers of patients, DLTs and responses at each dose.

- total_n_pts, total_n_tox, total_n_resp:

  Their averages per trial.

- stop_reason_percent:

  Percentage of trials ending for each reason, see
  [`boin12_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/boin12_simulate.md).

- est_method:

  The estimate of the utility used to select the OBD.

- settings:

  List of the design parameters.

- call:

  The matched call.

- trials, obd:

  The simulated trials and the selected OBDs, when `keep_trials` is
  `TRUE`.

## Details

The trials are simulated by
[`boin12_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/boin12_simulate.md)
and the OBD is selected by
[`boin12_select_obd`](https://gosukehommaex.github.io/simFastBOIN/reference/boin12_select_obd.md);
see there for the rules.

## References

Lin, R., Zhou, Y., Yan, F., Li, D. and Yuan, Y. (2020). BOIN12: Bayesian
optimal interval phase I/II trial design for utility-based dose finding
in immunotherapy and targeted therapies. JCO Precision Oncology, 4,
1393-1402.
[doi:10.1200/PO.20.00257](https://doi.org/10.1200/PO.20.00257)

## See also

[`sim_boin12_multi`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin12_multi.md),
[`boin12_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/boin12_simulate.md),
[`boin12_select_obd`](https://gosukehommaex.github.io/simFastBOIN/reference/boin12_select_obd.md)

## Examples

``` r
# Scenario 1 of Lin et al. (2020), 36 patients in cohorts of 3
oc <- sim_boin12(
  target = 0.35, target_resp = 0.25,
  p_true = c(0.05, 0.12, 0.27, 0.35, 0.50),
  p_resp = c(0.20, 0.35, 0.36, 0.37, 0.38),
  n_cohort = 12, cohort_size = 3, n_trials = 1000, seed = 123
)
oc
#> BOIN12 operating characteristics
#>   DLT limit            : 0.35
#>   response limit       : 0.25
#>   trials               : 1000
#>   cohorts              : 12 of size 3
#>   max sample size      : 36
#>   lambda_e / lambda_d  : 0.2763 / 0.4189
#>   utility scores       : 100, 40, 60, 0; benchmark 70.5
#>   OBD selection        : utility by posterior_mean
#> 
#>                         DL1  DL2  DL3  DL4  DL5 Total / No OBD
#> True DLT rate (%)       5.0 12.0 27.0 35.0 50.0               
#> True response rate (%) 20.0 35.0 36.0 37.0 38.0               
#> OBD selected (%)       22.2 47.7 20.4  8.7  1.0            0.0
#> Patients treated        7.8 12.5  8.4  5.3  2.1           36.0
#> Patients with DLT       0.4  1.5  2.3  1.9  1.0            7.1
#> Patients with response  1.6  4.4  3.1  1.9  0.8           11.7
#> 
#> Trials stopped early: 0%
```
