# Operating Characteristics of the BOIN12 Design Across Scenarios

Run
[`sim_boin12`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin12.md)
for several scenarios of DLT and response probabilities with the same
design and collect the results in one table.

## Usage

``` r
sim_boin12_multi(
  target,
  target_resp,
  scenarios,
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
  tie_break = c("higher", "random"),
  no_candidate = c("nearest", "stop"),
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

- scenarios:

  List of scenarios. Each element is a list with the components `p_true`
  and `p_resp`, and optionally `name`. Names of the list are used when
  `name` is absent. All scenarios must have the same number of doses.

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

- tie_break:

  Character scalar. Choice among the doses of equal desirability:
  `"higher"` (the default) takes the highest of them, the rule of the
  BOIN12 software according to the 'escalation' package; `"random"`
  takes one of them at random, as stated in Table S2 of the Supplement
  of Lin et al. (2020). See Details.

- no_candidate:

  Character scalar. What happens when none of the doses allowed for the
  next cohort is admissible: `"nearest"` (the default) goes to the
  admissible dose nearest to the current one; `"stop"` ends the trial
  without an OBD. See Details.

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

An object of class `boin12_oc_multi`, which is a list with components

- results:

  Named list of `boin12_oc` objects, one per scenario.

- summary_table:

  Data frame with six rows per scenario, as shown by
  [`print.boin12_oc_multi`](https://gosukehommaex.github.io/simFastBOIN/reference/print.boin12_oc_multi.md).

- scenario_names:

  Character vector of scenario names.

- n_doses:

  Number of doses.

- call:

  The matched call.

## Details

Every scenario is simulated with the same `seed`, so that a scenario
simulated here matches the same scenario simulated on its own with
[`sim_boin12`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin12.md).

## See also

[`sim_boin12`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin12.md)

## Examples

``` r
scenarios <- list(
  "Plateau" = list(p_true = c(0.05, 0.12, 0.27, 0.35, 0.50),
                   p_resp = c(0.20, 0.35, 0.36, 0.37, 0.38)),
  "Increasing" = list(p_true = c(0.05, 0.12, 0.27, 0.35, 0.50),
                      p_resp = c(0.01, 0.05, 0.30, 0.60, 0.60))
)
oc <- sim_boin12_multi(
  target = 0.35, target_resp = 0.25, scenarios = scenarios,
  n_cohort = 12, cohort_size = 3, n_trials = 200, seed = 1
)
oc
#> BOIN12 operating characteristics across 2 scenarios
#>   DLT limit            : 0.35
#>   response limit       : 0.25
#>   trials each          : 200
#>   max sample size      : 36
#>   ties / no candidate  : higher / nearest
#>   OBD selection        : utility by posterior_mean
#> 
#>    Scenario                   Item  DL1  DL2  DL3  DL4  DL5 Total / No OBD
#>     Plateau      True DLT rate (%)  5.0 12.0 27.0 35.0 50.0               
#>             True response rate (%) 20.0 35.0 36.0 37.0 38.0               
#>                   OBD selected (%) 16.5 49.5 22.5 10.5  1.0            0.0
#>                   Patients treated  6.8 12.6  8.7  5.7  2.1           36.0
#>                  Patients with DLT  0.4  1.6  2.4  1.9  1.1            7.5
#>             Patients with response  1.3  4.5  3.1  2.1  0.8           11.8
#>  Increasing      True DLT rate (%)  5.0 12.0 27.0 35.0 50.0               
#>             True response rate (%)  1.0  5.0 30.0 60.0 60.0               
#>                   OBD selected (%)  3.5  7.0 27.0 49.5  5.0            8.0
#>                   Patients treated  5.8  6.9  8.9 10.6  3.0           35.2
#>                  Patients with DLT  0.3  0.9  2.6  3.6  1.6            9.1
#>             Patients with response  0.0  0.4  2.7  6.4  1.8           11.3
```
