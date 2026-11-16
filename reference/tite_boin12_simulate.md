# Simulate TITE-BOIN12 Trials

Simulate trials of the TITE-BOIN12 design of Zhou et al. (2022), the
time-to-event version of BOIN12 in which patients arrive over time, the
DLT and the response are assessed over windows and the decisions use the
outcomes still pending through the approximated likelihood. Return the
patients, DLTs and responses at every dose. This is the simulation
engine behind
[`sim_tite_boin12`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_tite_boin12.md).

## Usage

``` r
tite_boin12_simulate(
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
  response of a patient, see `resp_cor_type`; with `"copula"` it also
  links the time to DLT and the time to response. Defaults to 0.

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

An object of class `boin12_trials`, as returned by
[`boin12_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/boin12_simulate.md),
with `design` set to `"TITE-BOIN12"` and one more component

- duration:

  Numeric vector of the duration of each trial: the time at which every
  patient has finished both assessments, or the time of the decision
  that stopped the trial because no dose was admissible.

The DLTs and responses are the final outcomes of all the patients.

## Details

The first patient of a cohort arrives at the time of the previous
decision, and the next decision is due when the next patient arrives.
While fewer than a share `min_completed_share` of the patients at the
current dose have finished both assessments, the decision waits for the
next assessment to end. After the last cohort, and when the trial stops
because `n_earlystop` is reached, the trial waits until every patient
has finished both assessments, and the OBD is selected from the complete
data, as in the example of Zhou et al. (2022).

At a decision, an outcome is observed once its event has occurred or its
window has ended, and pending otherwise. For each dose, with `w` the
fraction of the window followed by a pending patient, the DLT
probability is estimated by `pi_hat` = DLTs observed / (DLTs observed +
patients who completed the window without DLT + sum of `w`) (Equation 14
of the Supplement), and the response probability in the same way. The
quasi-events of BOIN12 count each patient with a pending outcome by the
expected utility, the pending outcome being replaced as set by
`pending_estimate` and the two outcomes taken as independent for this
purpose only. The desirability is that of BOIN12 with these quasi-events
and the number of patients treated, and the rules of
[`boin12_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/boin12_simulate.md)
are applied with `pi_hat` of the DLTs as the DLT rate of the current
dose. With every outcome observed, these are exactly the quantities of
BOIN12, so that `min_completed_share = 1` with `accrual = "fixed"` gives
the trials of
[`boin12_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/boin12_simulate.md)
with the same seed.

Two points are not settled by Zhou et al. (2022) and are decided here as
follows. Equation (11) replaces a pending outcome by its conditional
expectation, but the desirabilities of their Table S.3 follow from
`pi_hat` itself; both are offered by `pending_estimate`. The
admissibility of a dose, which they update at every decision, uses the
posterior of the approximated likelihood, Beta(1 + DLTs observed, 1 +
patients who completed the window without DLT + sum of `w`), and the
same for the responses, so that a response that is merely late does not
exclude a dose. TITE-BOIN-ET
([`tite_boin_et_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/tite_boin_et_simulate.md))
instead counts the pending patients as patients without the event, as
the 'boinet' package does.

Each patient draws two uniform variates as in
[`boin12_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/boin12_simulate.md);
the time to DLT and the time to response are quantiles of their
distributions at these variates, so that a DLT and a response occur
within their windows with the probabilities `p_true` and `p_resp`. The
times between arrivals draw further variates unless `accrual = "fixed"`.

## References

Zhou, Y., Lin, R., Lee, J. J., Li, D., Wang, L., Li, R. and Yuan, Y.
(2022). TITE-BOIN12: A Bayesian phase I/II trial design to find the
optimal biological dose with late-onset toxicity and efficacy.
Statistics in Medicine, 41(11), 1918-1931.
[doi:10.1002/sim.9337](https://doi.org/10.1002/sim.9337)

Lin, R., Zhou, Y., Yan, F., Li, D. and Yuan, Y. (2020). BOIN12: Bayesian
optimal interval phase I/II trial design for utility-based dose finding
in immunotherapy and targeted therapies. JCO Precision Oncology, 4,
1393-1402.
[doi:10.1200/PO.20.00257](https://doi.org/10.1200/PO.20.00257)

## See also

[`sim_tite_boin12`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_tite_boin12.md),
[`boin12_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/boin12_simulate.md),
[`boin12_select_obd`](https://gosukehommaex.github.io/simFastBOIN/reference/boin12_select_obd.md)

## Examples

``` r
# Scenario 1 of Zhou et al. (2022), windows of 45 and 60 days and three
# patients a month
trials <- tite_boin12_simulate(
  target = 0.35, target_resp = 0.25,
  p_true = c(0.05, 0.12, 0.27, 0.35, 0.50),
  p_resp = c(0.20, 0.35, 0.36, 0.37, 0.38),
  n_cohort = 12, cohort_size = 3, window = 45, resp_window = 60,
  accrual_rate = 3 / 30, n_trials = 200, seed = 123
)
trials
#> Simulated TITE-BOIN12 trials
#>   trials               : 200
#>   doses                : 5
#>   DLT limit            : 0.35
#>   response limit       : 0.25
#>   cohorts              : 12 of size 3
#>   max sample size      : 36
#> 
#> Average per trial:
#>                     DL1   DL2  DL3  DL4  DL5
#> Patients treated   7.17 12.18 8.30 5.70 2.65
#> DLTs observed      0.41  1.33 2.19 2.02 1.35
#> Responses observed 1.31  4.36 2.93 2.15 1.04
#> 
#> Stopping reason (%):
#> 
#> max_cohorts 
#>         100 
#> 
#> Average trial duration: 565.4
summary(trials$duration)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>   454.4   538.7   565.7   565.4   592.8   667.0 
```
