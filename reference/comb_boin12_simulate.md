# Simulate Comb-BOIN12 Trials

Simulate trials of the Comb-BOIN12 design of Lu et al. (2025), the
extension of BOIN12 to two drugs given together, which assigns each
cohort to the most desirable combination by a utility that weighs the
response against the DLT, and return the patients, DLTs and responses at
every combination. This is the simulation engine behind
[`sim_comb_boin12`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_comb_boin12.md).

## Usage

``` r
comb_boin12_simulate(
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

- seed:

  Integer scalar or `NULL`. Random seed. The state of the calling
  session is restored on exit. Defaults to 123.

## Value

An object of class `comb_boin12_trials`, which is a list with components

- n_pts, n_tox, n_resp, n_tox_resp:

  Integer arrays with dimensions trial by level of drug A by level of
  drug B, of the patients treated, the DLTs, the responses and the
  patients with both a DLT and a response.

- cohorts_used:

  Integer vector of the cohorts enrolled.

- stop_reason:

  Character vector: `"no_admissible_dose"`, `"n_earlystop"`,
  `"max_cohorts"` or `"no_admissible_candidate"` (only with
  `no_candidate = "stop"`).

- early_stop:

  Logical vector, `TRUE` for the trials stopped because no combination
  was admissible or, with `no_candidate = "stop"`, because no
  combination allowed for the next cohort was admissible. These trials
  select no OBD combination.

- design:

  `"Comb-BOIN12"`.

- boundary:

  List with `lambda_e`, `lambda_d`, `u_low` and `u_benchmark`.

- settings:

  List of the design parameters.

## Details

After each cohort, with `phat` the DLT rate at the current combination
(j, k): the trial stops when no combination is admissible; when `phat`
is at least `lambda_d`, the next combination is the most desirable of
the admissible combinations below (j, k), that is with levels at most j
and k, that have the highest sum of the two levels: (j - 1, k) or (j,
k - 1) when either is admissible, as in Lu et al. (2025), and lower ones
otherwise; the current one when there is none. When `phat` is above
`lambda_e` and at least `n_star` patients have been treated, the next
combination is the most desirable of (j, k), (j - 1, k) and (j, k - 1);
otherwise of (j, k) and its four neighbors (j - 1, k), (j, k - 1), (j +
1, k) and (j, k + 1). Only admissible combinations inside the matrix are
considered. When none of these is admissible, the admissible combination
nearest to the current one (by the sum of the differences of the two
levels) is chosen, then the one with the lowest sum of the levels, then
the most desirable. The desirability is that of BOIN12
([`boin12_decision_table`](https://gosukehommaex.github.io/simFastBOIN/reference/boin12_decision_table.md)).
Untried combinations have a DLT rate of 0 and are admissible unless
excluded by the DLT data of a lower combination.

With `no_candidate = "stop"`, the trial stops instead, without an OBD
combination and with the stopping reason `"no_admissible_candidate"`,
when none of these combinations is admissible and when `phat` is at
least `lambda_d` and no combination below the current one is admissible.
At the lowest combination, a `phat` of at least `lambda_d` still keeps
the next cohort there when it is admissible.

Lu et al. (2025) do not state how ties of desirability are broken nor
what happens when no candidate is admissible. With
`tie_break = "higher"`, among tied combinations those with the highest
sum of the levels are kept, so that the design escalates when an untried
combination is as desirable as the current one, and a tie that remains
(between (j + 1, k) and (j, k + 1), for example) is broken at random.
The allocation percentages published by Lu et al. (2025) agree with this
rule, and not with a random choice among all the tied combinations nor
with a fixed preference for one drug (see
`inst/reproduce/reproduce_Lu2025.R`). The two settings of `no_candidate`
could not be told apart from the published results. With one row or one
column, the trials are those of
[`boin12_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/boin12_simulate.md)
without its exploration rule; Lu et al. (2025) have no such rule.

When `psi00 + psi11` is 100, as for the utility of Lu et al. (2025), the
desirability does not depend on the number of patients with both a DLT
and a response, and combinations with the same numbers of patients, DLTs
and responses tie exactly, whatever the rounding of `psi00` and `psi11`.

Each patient draws two uniform variates in enrollment order, the first
for the DLT and the second for the response. A decision between tied
combinations draws one more variate when a tie remains after the rule of
`tie_break`.

## References

Lu, M., Zhang, J., Yuan, Y. and Lin, R. (2025). Comb-BOIN12: A
utility-based Bayesian optimal interval design for dose optimization in
cancer drug-combination trials. Statistics in Biopharmaceutical
Research, 17(2), 266-276.
[doi:10.1080/19466315.2024.2370403](https://doi.org/10.1080/19466315.2024.2370403)

Lin, R., Zhou, Y., Yan, F., Li, D. and Yuan, Y. (2020). BOIN12: Bayesian
optimal interval phase I/II trial design for utility-based dose finding
in immunotherapy and targeted therapies. JCO Precision Oncology, 4,
1393-1402.
[doi:10.1200/PO.20.00257](https://doi.org/10.1200/PO.20.00257)

## See also

[`sim_comb_boin12`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_comb_boin12.md),
[`comb_boin12_select_obd`](https://gosukehommaex.github.io/simFastBOIN/reference/comb_boin12_select_obd.md),
[`boin12_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/boin12_simulate.md)

## Examples

``` r
# Scenario 16 of Lu et al. (2025), with their utility (wT = 0.11)
p_true <- matrix(c(0.05, 0.10, 0.15, 0.15, 0.20, 0.45), nrow = 2)
p_resp <- matrix(c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60), nrow = 2)
trials <- comb_boin12_simulate(
  target = 0.33, target_resp = 0.2, p_true = p_true, p_resp = p_resp,
  n_cohort = 14, cohort_size = 3, n_trials = 200,
  psi00 = 100 * 0.11 / 1.11, psi11 = 100 / 1.11, benchmark_weight = 0.6,
  extrasafe = TRUE, offset = 0.10, seed = 123
)
trials
#> Simulated Comb-BOIN12 trials
#>   trials               : 200
#>   dose levels          : 2 of drug A, 3 of drug B
#>   DLT limit            : 0.33
#>   response limit       : 0.2
#>   cohorts              : 14 of size 3
#>   max sample size      : 42
#> 
#> Patients treated, average per trial:
#>      B1   B2    B3
#> A1 4.16 5.62 12.09
#> A2 3.52 8.23  7.98
#> 
#> DLTs observed, average per trial:
#>      B1   B2   B3
#> A1 0.16 0.88 2.50
#> A2 0.32 1.30 3.56
#> 
#> Responses observed, average per trial:
#>      B1   B2   B3
#> A1 0.45 1.68 6.18
#> A2 0.68 3.18 4.76
#> 
#> Stopping reason (%):
#> 
#>        max_cohorts no_admissible_dose 
#>                 99                  1 
```
