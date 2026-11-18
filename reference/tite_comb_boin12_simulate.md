# Simulate Trials of the TITE-Comb-BOIN12 Design

Simulate trials of the TITE-Comb-BOIN12 design of Lu et al. (2025) for
two drugs given together, in which the DLT and the response may be
observed late: the rules of
[`comb_boin12_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/comb_boin12_simulate.md)
applied at the decision times and with the interim estimates of
[`tite_boin12_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/tite_boin12_simulate.md).
Only the trial conduct is simulated; the OBD combination of each trial
is selected by
[`comb_boin12_select_obd`](https://gosukehommaex.github.io/simFastBOIN/reference/comb_boin12_select_obd.md)
(both are done by
[`sim_tite_comb_boin12`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_tite_comb_boin12.md)).

## Usage

``` r
tite_comb_boin12_simulate(
  target,
  target_resp,
  p_true,
  p_resp,
  n_cohort,
  cohort_size,
  window,
  accrual_rate,
  resp_window = window,
  accrual = c("uniform", "exponential", "fixed"),
  event_time = c("weibull", "uniform"),
  late_fraction = 0.5,
  resp_late_fraction = 0.5,
  min_completed_share = 0.5,
  pending_estimate = c("conditional", "marginal"),
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

- window:

  Numeric scalar. Length of the DLT assessment window.

- accrual_rate:

  Numeric scalar. Average number of patients arriving per unit of time,
  in the same unit as `window`.

- resp_window:

  Numeric scalar. Length of the response assessment window. Defaults to
  `window`.

- accrual:

  Character scalar. Distribution of the times between arrivals:
  `"uniform"` (the default) on `(0, 2 / accrual_rate)`, `"exponential"`
  for a Poisson process, or `"fixed"` for arrivals exactly
  `1 / accrual_rate` apart.

- event_time:

  Character scalar. Distribution of the times to DLT and to response:
  `"weibull"` (the default), with the share of the events in the second
  half of the window set by `late_fraction` and `resp_late_fraction`, or
  `"uniform"` over the window. With `"weibull"` the probabilities must
  be below 1.

- late_fraction:

  Numeric scalar strictly between 0 and 1. Share of the DLTs that occur
  in the second half of the DLT window under the Weibull distribution.
  Defaults to 0.5.

- resp_late_fraction:

  Numeric scalar strictly between 0 and 1. The same for the responses.
  Defaults to 0.5.

- min_completed_share:

  Numeric scalar from 0 to 1. Share of the patients at the current dose
  who must have finished both assessments before the next decision; the
  decision, and the accrual, wait until then. Defaults to 0.5, the rule
  of Zhou et al. (2022), who suspend the accrual when more than half of
  the patients at the current dose have a pending outcome. With 1 every
  decision waits for all the outcomes, as in BOIN12.

- pending_estimate:

  Character scalar. Value that replaces a pending outcome Y in the
  quasi-binomial likelihood: `"conditional"` (the default), the
  probability of the event given that it has not occurred during the
  time followed, `pi_hat * (1 - w) / (1 - pi_hat * w)`, with `w` the
  fraction of the window followed, as in Equation (11) of Zhou et al.
  (2022); or `"marginal"`, the estimate `pi_hat` itself, which
  reproduces the desirabilities of Table S.3 of their Supplement. See
  Details.

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

An object of class `comb_boin12_trials`, as returned by
[`comb_boin12_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/comb_boin12_simulate.md),
with `design` set to `"TITE-Comb-BOIN12"`, the counts of the final
outcomes of all patients, and one more component

- duration:

  Numeric vector of the trial durations, in the unit of `window`: the
  time of the last decision of a trial stopped because no combination
  was admissible, and otherwise the time at which every patient has
  finished both assessments.

The timing arguments are added to `settings`.

## Details

Patients arrive, are followed and are counted as in
[`tite_boin12_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/tite_boin12_simulate.md):
the first patient of a cohort arrives at the time of the previous
decision; the decision waits while fewer than a share
`min_completed_share` of the patients at the current combination have
finished both assessments; a pending outcome counts in the estimates of
the DLT and response rates, in the admissibility and in the quasi-events
of the desirability by the rules of Zhou et al. (2022). The decision
follows the rules of
[`comb_boin12_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/comb_boin12_simulate.md),
with the estimated DLT rate of the current combination, the
admissibility computed from the interim posteriors (the DLT rule of a
combination also excluding every combination with both levels at least
its own, and `extrasafe` lowering the cutoff at the lowest combination)
and the patients treated.

With every outcome observed at every decision these are the quantities
of Comb-BOIN12: with `accrual = "fixed"` and `min_completed_share = 1`
no variate is drawn for the times and the trials are exactly those of
[`comb_boin12_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/comb_boin12_simulate.md)
with the same seed. With one row or one column of combinations the
trials are those of
[`tite_boin12_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/tite_boin12_simulate.md)
without exploration.

Lu et al. (2025) state the time-to-event rules of TITE-BOIN12 for
TITE-Comb-BOIN12 without further detail; the accrual suspension and the
estimates are those of
[`tite_boin12_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/tite_boin12_simulate.md),
as in `inst/reproduce/reproduce_Lu2025.R`.

## References

Lu, M., Zhang, J., Yuan, Y. and Lin, R. (2025). Comb-BOIN12: A
utility-based Bayesian optimal interval design for dose optimization in
cancer drug-combination trials. Statistics in Biopharmaceutical
Research, 17(2), 266-276.
[doi:10.1080/19466315.2024.2370403](https://doi.org/10.1080/19466315.2024.2370403)

Zhou, Y., Lin, R., Lee, J. J., Li, D., Wang, L., Li, R. and Yuan, Y.
(2022). TITE-BOIN12: A Bayesian phase I/II trial design to find the
optimal biological dose with late-onset toxicity and efficacy.
Statistics in Medicine, 41(11), 1918-1931.
[doi:10.1002/sim.9337](https://doi.org/10.1002/sim.9337)

## See also

[`sim_tite_comb_boin12`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_tite_comb_boin12.md),
[`comb_boin12_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/comb_boin12_simulate.md),
[`tite_boin12_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/tite_boin12_simulate.md)

## Examples

``` r
# Scenario 16 of Lu et al. (2025), with their utility (wT = 0.11), windows
# of 1 and 2 months and two patients a month
p_true <- matrix(c(0.05, 0.10, 0.15, 0.15, 0.20, 0.45), nrow = 2)
p_resp <- matrix(c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60), nrow = 2)
trials <- tite_comb_boin12_simulate(
  target = 0.33, target_resp = 0.2, p_true = p_true, p_resp = p_resp,
  n_cohort = 14, cohort_size = 3, window = 1, resp_window = 2,
  accrual_rate = 2, event_time = "uniform", n_trials = 200,
  psi00 = 100 * 0.11 / 1.11, psi11 = 100 / 1.11, benchmark_weight = 0.6,
  extrasafe = TRUE, offset = 0.10, seed = 123
)
trials
#> Simulated TITE-Comb-BOIN12 trials
#>   trials               : 200
#>   dose levels          : 2 of drug A, 3 of drug B
#>   DLT limit            : 0.33
#>   response limit       : 0.2
#>   cohorts              : 14 of size 3
#>   max sample size      : 42
#> 
#> Patients treated, average per trial:
#>      B1   B2    B3
#> A1 3.80 5.79 10.66
#> A2 3.42 9.60  8.73
#> 
#> DLTs observed, average per trial:
#>      B1   B2   B3
#> A1 0.20 0.85 2.14
#> A2 0.29 1.40 3.90
#> 
#> Responses observed, average per trial:
#>      B1   B2   B3
#> A1 0.32 1.73 5.36
#> A2 0.72 3.90 5.32
#> 
#> Stopping reason (%):
#> 
#> max_cohorts 
#>         100 
#> 
#> Average trial duration: 26.3
```
