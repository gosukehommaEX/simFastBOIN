# Simulate TITE-BOIN-ET Trials

Run the TITE-BOIN-ET dose-finding algorithm of Takeda, Morita and Taguri
(2020) over many simulated trials, in which patients arrive over time
and the DLT and the response are assessed over windows, and return the
patient, DLT and response counts at every dose. This is the simulation
engine behind
[`sim_tite_boin_et`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_tite_boin_et.md),
exposed for users who want the raw trial data.

## Usage

``` r
tite_boin_et_simulate(
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
  resp_cor = 0.2,
  n_trials = 10000,
  start_dose = 1,
  n_earlystop = NULL,
  p_saf = NULL,
  p_tox = NULL,
  p_saf_resp = NULL,
  cutoff_eli = 0.95,
  cutoff_eli_resp = 0.99,
  p_eli = NULL,
  p_eli_resp = NULL,
  boinet_compatible = FALSE,
  prior = rep(1/6, 6),
  seed = 123
)
```

## Arguments

- target:

  Numeric scalar. Target DLT probability, for example 0.30.

- target_resp:

  Numeric scalar. Target response probability, for example 0.60.

- p_true:

  Numeric vector. True DLT probability at each dose level, in increasing
  dose order.

- p_resp:

  Numeric vector of the same length as `p_true`. True response
  probability at each dose level.

- n_cohort:

  Integer scalar. Number of cohorts in a trial.

- cohort_size:

  Integer scalar or vector. Number of patients per cohort. A scalar is
  used for every cohort. A vector shorter than `n_cohort` is padded with
  its last element and a longer one is truncated.

- window:

  Numeric scalar. Length of the DLT assessment window.

- accrual_rate:

  Numeric scalar. Average number of patients arriving per unit of time,
  in the same unit as `window`. The `accrual` argument of the 'boinet'
  package is the average time between arrivals, its inverse.

- resp_window:

  Numeric scalar. Length of the response assessment window. Defaults to
  `window`.

- accrual:

  Character scalar. Distribution of the times between arrivals:
  `"uniform"` (the default, as in 'boinet') on `(0, 2 / accrual_rate)`,
  `"exponential"` for a Poisson process, or `"fixed"` for arrivals
  exactly `1 / accrual_rate` apart.

- event_time:

  Character scalar. Distribution of the times to DLT and to response:
  `"weibull"` (the default), with the share of the events in the second
  half of the window set by `late_fraction` and `resp_late_fraction`, or
  `"uniform"` on `(0, window / p)` for an event of probability `p`. With
  `"weibull"` the probabilities must be below 1.

- late_fraction:

  Numeric scalar strictly between 0 and 1. Share of the DLTs that occur
  in the second half of the DLT window under the Weibull distribution.
  Defaults to 0.5.

- resp_late_fraction:

  Numeric scalar strictly between 0 and 1. The same for the responses.
  Defaults to 0.5.

- resp_cor:

  Numeric scalar strictly between -1 and 1. Correlation of the Gaussian
  copula that links the time to DLT and the time to response of a
  patient. Defaults to 0.2, the default of the 'boinet' package.

- n_trials:

  Integer scalar. Number of trials to simulate. Defaults to 10000.

- start_dose:

  Integer scalar. Dose level for the first cohort. Defaults to 1.

- n_earlystop:

  Integer scalar or `NULL`. The trial stops once this many patients have
  been treated at any dose. Defaults to `NULL`, the maximum sample size,
  which is the default of the 'boinet' package.

- p_saf:

  Numeric scalar. Highest DLT probability deemed subtherapeutic.
  Defaults to `0.1 * target`, the value recommended by Takeda et al.
  (2018), which is far below the `0.6 * target` of the BOIN design.

- p_tox:

  Numeric scalar. Lowest DLT probability deemed overly toxic. Defaults
  to `1.4 * target`.

- p_saf_resp:

  Numeric scalar. Highest response probability deemed subtherapeutic.
  Defaults to `0.6 * target_resp`.

- cutoff_eli:

  Numeric scalar. A dose is no longer admissible when the posterior
  probability that its DLT rate is below `p_eli` falls below
  `1 - cutoff_eli`. Defaults to 0.95.

- cutoff_eli_resp:

  Numeric scalar. A dose is no longer admissible when the posterior
  probability that its response rate is above `p_eli_resp` falls below
  `1 - cutoff_eli_resp`. Defaults to 0.99.

- p_eli:

  Numeric scalar or `NULL`. DLT probability of the rule above. Defaults
  to `NULL`, which uses `target` as the 'boinet' package does. Takeda et
  al. (2018) used 0.40 with a cutoff of 0.90.

- p_eli_resp:

  Numeric scalar or `NULL`. Response probability of the rule above.
  Defaults to `NULL`, which uses `p_saf_resp` as the 'boinet' package
  does. Takeda et al. (2018) used 0.20 with a cutoff of 0.90.

- boinet_compatible:

  Logical scalar. Apply the rules of the 'boinet' package (version
  1.6.0) exactly. Defaults to `FALSE`, which keeps to Takeda et al.
  (2018), who do not allow dose skipping during escalation; see Details.

- prior:

  Numeric vector of length 6. Prior probabilities of the six hypotheses.
  Defaults to `rep(1 / 6, 6)`, the noninformative choice of the article.
  The third hypothesis, under which every decision is correct, does not
  enter the criterion.

- seed:

  Integer scalar or `NULL`. Random seed. When `NULL` the current state
  of the random number generator is used. The state of the calling
  session is restored on exit. Defaults to 123.

## Value

An object of class `boin_et_trials`, as returned by
[`boin_et_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_et_simulate.md),
with `design` set to `"TITE-BOIN-ET"` and one more component

- duration:

  Numeric vector of the time of the last decision of each trial, which
  is the trial duration when all cohorts are enrolled.

The DLTs and responses are those observed at the last decision.

## Details

The first patient of a cohort arrives at the time of the previous
decision, and the next decision is due when the next patient arrives;
after the last cohort it is due once every patient has been followed
over both windows. While fewer than half of the patients at the current
dose have finished both assessments, the decision waits for the next
assessment to end. The DLT and response rates of a dose then count a
patient still under assessment by the fraction of the window followed so
far (Takeda et al., 2020), and the decision rule of BOIN-ET is applied
to these rates. The admissibility of the doses, which counts the
patients still under assessment as patients without the event, and the
moves between doses follow
[`boin_et_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_et_simulate.md),
including the changes from 'boinet' unless `boinet_compatible` is
`TRUE`. These are the rules of
[`boinet::tite.boinet()`](https://rdrr.io/pkg/boinet/man/tite.boinet.html).

The two event times of a patient come from their marginal distributions
linked by a Gaussian copula, so that a DLT and a response occur within
their windows with the probabilities `p_true` and `p_resp`. Each patient
draws two uniform variates, so the trials agree with those of 'boinet'
in distribution, not one by one.

## References

Takeda, K., Morita, S. and Taguri, M. (2020). TITE-BOIN-ET:
Time-to-event Bayesian optimal interval design to accelerate
dose-finding based on both efficacy and toxicity outcomes.
Pharmaceutical Statistics, 19, 335-349.
[doi:10.1002/pst.1995](https://doi.org/10.1002/pst.1995)

Yamaguchi, Y., Takeda, K. and Maruo, K. (2025). boinet: An R Package for
Bayesian Optimal Interval Design for Dose Finding Based on Both Efficacy
and Toxicity Outcomes. The R Journal, 17(2), 15-34.
[doi:10.32614/RJ-2025-012](https://doi.org/10.32614/RJ-2025-012)

## See also

[`sim_tite_boin_et`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_tite_boin_et.md),
[`boin_et_select_obd`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_et_select_obd.md)

## Examples

``` r
trials <- tite_boin_et_simulate(
  target = 0.30, target_resp = 0.60,
  p_true = c(0.05, 0.15, 0.25, 0.30, 0.35, 0.40),
  p_resp = c(0.05, 0.35, 0.55, 0.65, 0.75, 0.70),
  n_cohort = 12, cohort_size = 3, window = 30, resp_window = 45,
  accrual_rate = 1 / 10, n_trials = 200, seed = 123
)
trials
#> Simulated TITE-BOIN-ET trials
#>   trials               : 200
#>   doses                : 6
#>   target DLT rate      : 0.3
#>   target response rate : 0.6
#>   cohorts              : 12 of size 3
#>   max sample size      : 36
#> 
#> Average per trial:
#>                     DL1  DL2   DL3  DL4  DL5  DL6
#> Patients treated   4.28 9.80 13.32 6.38 1.64 0.14
#> DLTs observed      0.24 1.39  3.48 2.02 0.58 0.06
#> Responses observed 0.22 3.49  7.50 4.18 1.24 0.10
#> 
#> Stopping reason (%):
#> 
#>              max_cohorts       no_admissible_dose no_lower_admissible_dose 
#>                       96                        1                        3 
#> 
#> Average trial duration: 462.5
summary(trials$duration)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>   219.4   446.4   472.0   462.5   490.8   558.7 
```
