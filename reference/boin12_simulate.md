# Simulate BOIN12 Trials

Simulate trials of the BOIN12 design of Lin et al. (2020), which assigns
each cohort to the most desirable dose by a utility that weighs the
response against the DLT, and return the patients, DLTs and responses at
every dose. This is the simulation engine behind
[`sim_boin12`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin12.md).

## Usage

``` r
boin12_simulate(
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
  tie_break = c("higher", "random"),
  no_candidate = c("nearest", "stop"),
  explore_from = c("any", "admissible"),
  resp_cor = 0,
  resp_cor_type = c("copula", "binary"),
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

- explore_from:

  Character scalar. Doses from which the exploration rule (see
  `n_explore`) escalates: `"any"` (the default) from any dose, as in the
  'escalation' package; `"admissible"` only from an admissible dose, as
  the results published by Lin et al. (2020) indicate. See Details.

- resp_cor:

  Numeric scalar between -1 and 1. Correlation between the DLT and the
  response of a patient, see `resp_cor_type`. Defaults to 0.

- resp_cor_type:

  Character scalar. `"copula"` (the default) links the two outcomes by a
  Gaussian copula with correlation `resp_cor`, the bivariate normal
  model of Lin et al. (2020); `"binary"` sets the correlation of the two
  binary outcomes, as in
  [`boin_et_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_et_simulate.md).

- seed:

  Integer scalar or `NULL`. Random seed. The state of the calling
  session is restored on exit. Defaults to 123.

## Value

An object of class `boin12_trials`, which is a list with components

- n_pts, n_tox, n_resp, n_tox_resp:

  Integer matrices, one row per trial and one column per dose, of the
  patients treated, the DLTs, the responses and the patients with both a
  DLT and a response.

- cohorts_used:

  Integer vector of the cohorts enrolled.

- stop_reason:

  Character vector: `"no_admissible_dose"`, `"n_earlystop"`,
  `"max_cohorts"` or `"no_admissible_candidate"` (only with
  `no_candidate = "stop"`).

- early_stop:

  Logical vector, `TRUE` for the trials stopped because no dose was
  admissible or, with `no_candidate = "stop"`, because no dose allowed
  for the next cohort was admissible. These trials select no OBD.

- design:

  `"BOIN12"`.

- boundary:

  List with `lambda_e`, `lambda_d`, `u_low` and `u_benchmark`.

- settings:

  List of the design parameters.

## Details

After each cohort, with `phat` the DLT rate at the current dose: the
trial stops when no dose is admissible; when `n_explore` patients have
been treated at the current dose, `phat` is below `lambda_d` and the
dose above is untried, the dose is escalated (with
`explore_from = "admissible"`, only when the current dose is
admissible); otherwise, when `phat` is at least `lambda_d`, the next
dose is the highest admissible dose below the current one (the current
one if there is none); when `phat` is above `lambda_e` and at least
`n_star` patients have been treated, the next dose is the more desirable
of the current dose and the dose below; otherwise it is the most
desirable of the dose below, the current dose and the dose above. Only
admissible doses are considered. Ties go to the higher dose, as in the
BOIN12 software according to the 'escalation' package, or with
`tie_break = "random"` to one of the tied doses drawn at random. When
none of these doses is admissible, the admissible dose nearest to the
current one is chosen, ties to the lower dose. The desirability is that
of
[`boin12_decision_table`](https://gosukehommaex.github.io/simFastBOIN/reference/boin12_decision_table.md).

With `no_candidate = "stop"`, the trial stops instead, without an OBD
and with the stopping reason `"no_admissible_candidate"`, when none of
these doses is admissible and when `phat` is at least `lambda_d` and no
dose below the current one is admissible. At the lowest dose, a `phat`
of at least `lambda_d` still keeps the next cohort there when that dose
is admissible.

Lin et al. (2020) leave three points open, decided here as follows. A
dose made inadmissible by its DLT data also makes every higher dose
inadmissible, as in BOIN and in Lu et al. (2025), whereas the
'escalation' package only checks the doses next to the current one. A
de-escalation never goes to an inadmissible dose; 'escalation' moves to
the dose below whatever its data. The trial stops when no dose is
admissible; 'escalation' goes on at the lowest dose when its DLT rate is
high. The article breaks ties at random; the software chooses the higher
dose. `tie_break = "random"`, `no_candidate = "stop"` and
`explore_from = "admissible"` are the settings compared with the results
published by Lin et al. (2020) in `inst/reproduce/reproduce_Lin2020.R`;
the last two rules are not stated in the article and were inferred from
those results.

Each patient draws two uniform variates in enrollment order, the first
for the DLT and the second for the response. With
`tie_break = "random"`, a decision among tied doses draws one more
variate, so that the trials depart from those with `"higher"` after the
first tie.

## References

Lin, R., Zhou, Y., Yan, F., Li, D. and Yuan, Y. (2020). BOIN12: Bayesian
optimal interval phase I/II trial design for utility-based dose finding
in immunotherapy and targeted therapies. JCO Precision Oncology, 4,
1393-1402.
[doi:10.1200/PO.20.00257](https://doi.org/10.1200/PO.20.00257)

Zhou, Y., Lin, R., Lee, J. J., Li, D., Wang, L., Li, R. and Yuan, Y.
(2022). TITE-BOIN12: A Bayesian phase I/II trial design to find the
optimal biological dose with late-onset toxicity and efficacy.
Statistics in Medicine, 41(11), 1918-1931.
[doi:10.1002/sim.9337](https://doi.org/10.1002/sim.9337)

Lu, M., Zhang, J., Yuan, Y. and Lin, R. (2025). Comb-BOIN12: A
utility-based Bayesian optimal interval design for dose optimization in
cancer drug-combination trials. Statistics in Biopharmaceutical
Research, 17(2), 266-276.
[doi:10.1080/19466315.2024.2370403](https://doi.org/10.1080/19466315.2024.2370403)

## See also

[`sim_boin12`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin12.md),
[`boin12_select_obd`](https://gosukehommaex.github.io/simFastBOIN/reference/boin12_select_obd.md),
[`boin12_decision_table`](https://gosukehommaex.github.io/simFastBOIN/reference/boin12_decision_table.md)

## Examples

``` r
# Scenario 1 of Lin et al. (2020)
trials <- boin12_simulate(
  target = 0.35, target_resp = 0.25,
  p_true = c(0.05, 0.12, 0.27, 0.35, 0.50),
  p_resp = c(0.20, 0.35, 0.36, 0.37, 0.38),
  n_cohort = 12, cohort_size = 3, n_trials = 200, seed = 123
)
trials
#> Simulated BOIN12 trials
#>   trials               : 200
#>   doses                : 5
#>   DLT limit            : 0.35
#>   response limit       : 0.25
#>   cohorts              : 12 of size 3
#>   max sample size      : 36
#> 
#> Average per trial:
#>                     DL1   DL2  DL3  DL4  DL5
#> Patients treated   8.18 12.38 8.30 5.10 2.06
#> DLTs observed      0.47  1.44 2.19 1.88 1.03
#> Responses observed 1.79  4.23 3.02 1.85 0.81
#> 
#> Stopping reason (%):
#> 
#> max_cohorts 
#>         100 
```
