# Simulate TITE-gBOIN-ET Trials

Run the TITE-gBOIN-ET dose-finding algorithm of Takeda, Yamaguchi,
Taguri and Morita (2023) over many simulated trials, for toxicity and
efficacy outcomes with several categories that each carry a score, in
which patients arrive over time and the outcomes are assessed over
windows, and return the patients, their categories and the total scores
at every dose. This is the simulation engine behind
[`sim_tite_gboin_et`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_tite_gboin_et.md),
exposed for users who want the raw trial data.

## Usage

``` r
tite_gboin_et_simulate(
  target,
  target_resp,
  p_grade,
  p_resp_grade,
  score,
  score_resp,
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
  min_completed_share = 0.5,
  n_trials = 10000,
  start_dose = 1,
  n_earlystop = NULL,
  p_saf = NULL,
  p_tox = NULL,
  p_saf_resp = NULL,
  lambda1 = NULL,
  lambda2 = NULL,
  eta1 = NULL,
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

  Numeric scalar. Target value of the mean normalized toxicity score,
  for example `0.47 / 1.5` for a target equivalent toxicity score of
  0.47 with a largest score of 1.5.

- target_resp:

  Numeric scalar. Target value of the mean normalized efficacy score,
  for example `0.95 / 1.5` for a target equivalent efficacy score of
  0.95 with a largest score of 1.5.

- p_grade:

  Numeric matrix. True probability of each toxicity category (rows) at
  each dose (columns, in increasing dose order). Every column sums to 1.
  The first row is the category without a toxicity event.

- p_resp_grade:

  Numeric matrix with one column per dose. True probability of each
  efficacy category (rows). Every column sums to 1. The first row is the
  category without an efficacy event.

- score:

  Numeric vector, one value per row of `p_grade`. Score of each toxicity
  category, for example `c(0, 0.5, 1, 1.5)` for grades 0 or 1, 2, 3
  and 4. The scores are divided by their maximum.

- score_resp:

  Numeric vector, one value per row of `p_resp_grade`. Score of each
  efficacy category, for example `c(0, 1, 1.5)` for stable or
  progressive disease, partial response and complete response. The
  scores are divided by their maximum.

- n_cohort:

  Integer scalar. Number of cohorts in a trial.

- cohort_size:

  Integer scalar or vector. Number of patients per cohort. A scalar is
  used for every cohort. A vector shorter than `n_cohort` is padded with
  its last element and a longer one is truncated.

- window:

  Numeric scalar. Length of the toxicity assessment window.

- accrual_rate:

  Numeric scalar. Average number of patients arriving per unit of time,
  in the same unit as `window`. The `accrual` argument of the 'boinet'
  package is the average time between arrivals, its inverse.

- resp_window:

  Numeric scalar. Length of the efficacy assessment window. Defaults to
  `window`.

- accrual:

  Character scalar. Distribution of the times between arrivals:
  `"uniform"` (the default, as in 'boinet') on `(0, 2 / accrual_rate)`,
  `"exponential"` for a Poisson process, or `"fixed"` for arrivals
  exactly `1 / accrual_rate` apart.

- event_time:

  Character scalar. Distribution of the times to the toxicity event and
  to the efficacy event: `"weibull"` (the default), with the share of
  the events in the second half of the window set by `late_fraction` and
  `resp_late_fraction`, or `"uniform"` on `(0, window / p)` for an event
  of probability `p`. With `"weibull"` the first category of both
  outcomes, the one without an event, must have a positive probability
  at every dose.

- late_fraction:

  Numeric scalar strictly between 0 and 1. Share of the toxicity events
  that occur in the second half of the toxicity window under the Weibull
  distribution. Defaults to 0.5.

- resp_late_fraction:

  Numeric scalar strictly between 0 and 1. The same for the efficacy
  events. Defaults to 0.5.

- resp_cor:

  Numeric scalar strictly between -1 and 1. Correlation of the Gaussian
  copula that links the time to the toxicity event and the time to the
  efficacy event of a patient. Defaults to 0.2, the default of the
  'boinet' package.

- min_completed_share:

  Numeric scalar from 0 to 1. Share of the patients at the current dose
  who must have finished both assessments before the next decision; the
  decision waits until then. Defaults to 0.5, the rule of Takeda et
  al. (2023) and of the 'boinet' package.

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

- lambda1, lambda2, eta1:

  Numeric scalars or `NULL`. Boundaries of the mean toxicity score and
  of the mean efficacy score. Defaults to `NULL`, the boundaries of
  [`boin_et_boundary`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/boin_et_boundary.md)
  for the targets. Give all three to use other values, such as those
  printed in an article.

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
[`gboin_et_simulate`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/gboin_et_simulate.md),
with `design` set to `"TITE-gBOIN-ET"` and one more component

- duration:

  Numeric vector of the time of the last decision of each trial, which
  is the trial duration when all cohorts are enrolled.

The categories and the total scores are those observed at the last
decision; a patient still under assessment is counted in the first
category.

## Details

The TITE-gBOIN-ET design is the gBOIN-ET design
([`gboin_et_simulate`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/gboin_et_simulate.md))
for outcomes that are observed late: the next cohort is enrolled while
earlier patients are still assessed. The first patient of a cohort
arrives at the time of the previous decision, and the next decision is
due when the next patient arrives; after the last cohort it is due once
every patient has been followed over both windows. While fewer than a
share `min_completed_share` of the patients at the current dose have
finished both assessments, the decision waits for the next assessment to
end.

At a decision, the mean toxicity score of a dose is the total of the
normalized scores observed so far divided by the effective number of
patients: the patients whose toxicity assessment has ended plus, for
every patient still under assessment, the fraction of the window
followed so far (Takeda et al., 2023). The mean efficacy score is formed
in the same way, and the decision rule of gBOIN-ET is applied to the two
means of the current dose. The admissibility of a dose uses the total
scores observed so far and the number of patients treated at the dose,
so that a patient still under assessment counts as a patient without an
event, and the moves between doses follow
[`boin_et_simulate`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/boin_et_simulate.md),
including the changes from 'boinet' unless `boinet_compatible` is
`TRUE`. These are the rules of
[`boinet::tite.gboinet()`](https://rdrr.io/pkg/boinet/man/tite.gboinet.html).

A patient has a toxicity event when the time to the event falls within
the window, which happens with the probability of the categories above
the first, and an efficacy event likewise. The two event times come from
their marginal distributions linked by a Gaussian copula, and the
category of a patient with an event is drawn from the probabilities of
the categories above the first, independently of everything else; it is
observed at the time of the event, and a patient without an event is
known to have none when the window ends. Each patient draws one uniform
variate for the toxicity event, one for the efficacy event, and one more
for each outcome with more than two categories. With two categories of
scores 0 and 1 for both outcomes the trials are, for the same seed,
those of
[`tite_boin_et_simulate`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/tite_boin_et_simulate.md);
with arrivals a fixed time apart that is longer than both windows, they
are those of
[`gboin_et_simulate`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/gboin_et_simulate.md).
The 'boinet' package draws its random numbers differently, so the trials
agree with it in distribution, not one by one.

The program of the supporting information of Takeda et al. (2023)
differs from these rules in two points. Where one dose among the current
dose and its neighbors has the highest mean efficacy score, the program
passes that dose to [`sample()`](https://rdrr.io/r/base/sample.html),
which then draws from all doses up to it, so that the trial can move to
any lower dose; and the program excludes doses with the effective
numbers of patients in the place of the numbers of patients. Neither is
described in the article, and 'boinet' has neither. The first changes
the published operating characteristics substantially; see
`inst/reproduce/reproduce_Takeda2023.R`.

The mean score of a dose is computed as the total of the scores as
given, divided by the largest score times the effective number of
patients, as in
[`gboin_et_simulate`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/gboin_et_simulate.md);
when no patient of the dose is still under assessment, doses whose mean
scores are equal as fractions tie on every platform.

Possible values of `stop_reason` are those of
[`boin_et_simulate`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/boin_et_simulate.md).

## References

Takeda, K., Yamaguchi, Y., Taguri, M. and Morita, S. (2023).
TITE-gBOIN-ET: Time-to-event generalized Bayesian optimal interval
design to accelerate dose-finding accounting for ordinal graded efficacy
and toxicity outcomes. Biometrical Journal, 65(7), 2200265.
[doi:10.1002/bimj.202200265](https://doi.org/10.1002/bimj.202200265)

Yamaguchi, Y., Takeda, K. and Maruo, K. (2025). boinet: An R Package for
Bayesian Optimal Interval Design for Dose Finding Based on Both Efficacy
and Toxicity Outcomes. The R Journal, 17(2), 15-34.
[doi:10.32614/RJ-2025-012](https://doi.org/10.32614/RJ-2025-012)

## See also

[`sim_tite_gboin_et`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_tite_gboin_et.md),
[`gboin_et_select_obd`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/gboin_et_select_obd.md),
[`gboin_et_simulate`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/gboin_et_simulate.md),
[`tite_boin_et_simulate`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/tite_boin_et_simulate.md)

## Examples

``` r
# Toxicity grades 0 or 1, 2, 3 and 4, and progressive disease, stable
# disease, partial response and complete response, at four doses, with a
# 30-day toxicity window, a 60-day efficacy window and one patient every
# 10 days on average
p_grade <- rbind(
  c(0.94, 0.79, 0.62, 0.50),
  c(0.05, 0.15, 0.20, 0.20),
  c(0.01, 0.05, 0.15, 0.25),
  c(0.00, 0.01, 0.03, 0.05)
)
p_resp_grade <- rbind(
  c(0.64, 0.45, 0.20, 0.05),
  c(0.30, 0.40, 0.40, 0.15),
  c(0.05, 0.10, 0.20, 0.35),
  c(0.01, 0.05, 0.20, 0.45)
)
trials <- tite_gboin_et_simulate(
  target = 0.47 / 1.5, target_resp = 1.75 / 3,
  p_grade = p_grade, p_resp_grade = p_resp_grade,
  score = c(0, 0.5, 1, 1.5), score_resp = c(0, 0.25, 1, 3),
  n_cohort = 12, cohort_size = 3, window = 30, resp_window = 60,
  accrual_rate = 1 / 10, n_trials = 200, seed = 123
)
trials
#> Simulated TITE-gBOIN-ET trials
#>   trials               : 200
#>   doses                : 4
#>   target toxicity score: 0.3133333
#>   target efficacy score: 0.5833333
#>   cohorts              : 12 of size 3
#>   max sample size      : 36
#> 
#> Average per trial:
#>                   DL1  DL2   DL3   DL4
#> Patients treated 3.15 4.16 10.05 18.64
#> Toxicity score   0.07 0.39  1.94  5.56
#> Efficacy score   0.16 0.53  3.06 10.77
#> 
#> Stopping reason (%):
#> 
#> max_cohorts 
#>         100 
#> 
#> Average trial duration: 499.9
summary(trials$duration)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>   405.0   476.0   502.4   499.9   525.8   596.3 
```
