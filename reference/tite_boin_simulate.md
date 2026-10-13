# Simulate TITE-BOIN Trials

Run the time-to-event BOIN (TITE-BOIN) design over many simulated
trials, in which patients arrive over time and new patients are treated
while the DLT assessment of earlier patients is still pending. Return
the patient and DLT counts at every dose together with the duration of
each trial and the time spent with accrual suspended.

## Usage

``` r
tite_boin_simulate(
  target,
  p_true,
  n_cohort,
  cohort_size,
  window,
  accrual_rate,
  method = c("imputation", "ess"),
  accrual = c("exponential", "uniform", "fixed"),
  dlt_time = c("weibull", "uniform"),
  late_fraction = 0.5,
  prior_weights = c(1, 1, 1)/3,
  max_pending_ratio = NULL,
  min_completed = NULL,
  min_follow_up = 0,
  n_trials = 10000,
  start_dose = 1,
  n_earlystop = 18,
  p_saf = NULL,
  p_tox = NULL,
  cutoff_eli = 0.95,
  extrasafe = FALSE,
  offset = 0.05,
  n_earlystop_rule = c("with_stay", "simple"),
  seed = 123
)
```

## Arguments

- target:

  Numeric scalar. Target DLT probability, for example 0.30.

- p_true:

  Numeric vector. True DLT probability at each dose level, in increasing
  dose order. The probability of a DLT within the assessment window.

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
  in the same unit as `window`.

- method:

  Character scalar, `"imputation"` (the default) or `"ess"`. How the
  pending patients enter the estimate of the DLT probability. See
  [`tite_boin_decision_table`](https://gosukehommaex.github.io/simFastBOIN/reference/tite_boin_decision_table.md).

- accrual:

  Character scalar. Distribution of the times between arrivals:
  `"exponential"` (the default) for a Poisson process, `"uniform"` for a
  uniform distribution on `(0, 2 / accrual_rate)`, as in the code of Lin
  and Yuan (2020), or `"fixed"` for arrivals exactly `1 / accrual_rate`
  apart.

- dlt_time:

  Character scalar. Distribution of the time to DLT within the window:
  `"weibull"` (the default), with the share of DLTs in the second half
  of the window set by `late_fraction`, or `"uniform"`.

- late_fraction:

  Numeric scalar strictly between 0 and 1. Share of the DLTs that occur
  in the second half of the assessment window under the Weibull
  distribution. Defaults to 0.5, the setting of Yuan et al. (2018) and
  Lin and Yuan (2020).

- prior_weights:

  Numeric vector of three non-negative values. Prior probabilities that
  a DLT occurs in the first, second and last third of the assessment
  window, used to weight the follow-up of the pending patients in STFT
  (Yuan et al., 2018, Supplementary Appendix D) and ESS (Lin and Yuan,
  2020, Supplementary S1). They are rescaled to sum to one. Defaults to
  equal weights, which give the unweighted STFT. The decision table does
  not change with the weights; only the value compared with its
  boundaries does. This describes the design and is separate from
  `dlt_time` and `late_fraction`, which describe how the simulated DLTs
  actually occur.

- max_pending_ratio:

  Numeric scalar or `NULL`. Suspension rule on the proportion of pending
  patients. See
  [`tite_boin_decision_table`](https://gosukehommaex.github.io/simFastBOIN/reference/tite_boin_decision_table.md).

- min_completed:

  Integer scalar or `NULL`. Number of completed assessments required for
  escalation. See
  [`tite_boin_decision_table`](https://gosukehommaex.github.io/simFastBOIN/reference/tite_boin_decision_table.md).

- min_follow_up:

  Numeric scalar between 0 and 1. Shortest follow-up, as a fraction of
  `window`, that every pending patient at the current dose must have
  reached before the dose is escalated; accrual is suspended until then.
  Defaults to 0, which never suspends on this ground. See
  [`tite_boin_decision_table`](https://gosukehommaex.github.io/simFastBOIN/reference/tite_boin_decision_table.md).

- n_trials:

  Integer scalar. Number of trials to simulate. Defaults to 10000.

- start_dose:

  Integer scalar. Dose level for the first cohort. Defaults to 1.

- n_earlystop:

  Integer scalar. The trial stops once this many patients have been
  treated at the current dose and the design would stay there. Defaults
  to 18.

- p_saf:

  Numeric scalar. Highest DLT probability deemed subtherapeutic.
  Defaults to `0.6 * target`.

- p_tox:

  Numeric scalar. Lowest DLT probability deemed overly toxic. Defaults
  to `1.4 * target`.

- cutoff_eli:

  Numeric scalar. Posterior probability cutoff for dose elimination.
  Defaults to 0.95.

- extrasafe:

  Logical scalar. Apply the stricter safety stopping rule at the lowest
  dose. Defaults to `FALSE`.

- offset:

  Numeric scalar between 0 and 0.5. Amount by which `cutoff_eli` is
  relaxed for the safety stopping rule. Defaults to 0.05.

- n_earlystop_rule:

  Character scalar, either `"with_stay"` or `"simple"`. See
  [`boin_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_simulate.md).

- seed:

  Integer scalar or `NULL`. Random seed. When `NULL` the current state
  of the random number generator is used. The state of the calling
  session is restored on exit. Defaults to 123.

## Value

An object of class `tite_boin_trials`, which is a list with components

- n_pts:

  Integer matrix, `n_trials` by number of doses, of patients treated.

- n_tox:

  Integer matrix of the same shape, of DLTs, counted once every patient
  has completed the assessment.

- eliminated:

  Logical matrix of the same shape, marking doses eliminated during the
  trial.

- cohorts_used:

  Integer vector of cohorts treated in each trial.

- stop_reason:

  Character vector giving why each trial ended, with the values of
  [`boin_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_simulate.md).

- duration:

  Numeric vector. Time from the first arrival until every enrolled
  patient has completed the assessment.

- n_suspensions:

  Integer vector. Number of decisions at which accrual had to be
  suspended.

- time_suspended:

  Numeric vector. Total time with accrual suspended.

- boundary:

  The `boin_boundary` object used.

- settings:

  List of the design parameters.

## Details

Patients arrive one at a time. The dose for a cohort is decided when its
first patient arrives, from the data available at that moment, and the
other patients of the cohort receive the same dose. The decision follows
[`tite_boin_decision_table`](https://gosukehommaex.github.io/simFastBOIN/reference/tite_boin_decision_table.md).
A patient with a DLT completes the assessment when the DLT occurs, and a
patient without one at the end of the window. When the decision is to
suspend accrual, the arriving patient waits until the next pending
patient at the current dose completes the assessment, or, with a
positive `min_follow_up`, until the last pending patient to arrive at
the current dose has been followed for that fraction of the window,
whichever comes first, and the decision is taken again at that moment.
Later patients arrive after the one who waited.

At every decision the rules are applied in the order used by
[`boin_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_simulate.md):
elimination at the current dose and the extra safety rule, early
stopping on `n_earlystop`, then the dose transition. Because DLTs can be
observed after a dose has been left, the elimination rule is also
checked at the other doses; a dose that meets it is eliminated together
with all higher doses, and the trial moves below it. Elimination counts
every treated patient, the pending ones as without DLT.

One uniform variate is drawn from R's random number stream per patient,
in enrollment order and in the cohort blocks of
[`boin_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_simulate.md).
It decides whether the patient has a DLT and, under the Weibull
distribution, gives the time to DLT through the inverse distribution
function, which falls inside the window exactly when there is a DLT. The
times between arrivals come from a second generator (xoshiro256\*\*)
seeded with `seed`, or with one draw from R's stream when `seed` is
`NULL`. As a result, whenever no patient is pending at any decision, for
example with `accrual = "fixed"` and `window` shorter than
`1 / accrual_rate`, the trials are identical to those of
[`boin_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_simulate.md)
with the same seed.

The follow-up of a pending patient counts towards STFT as the share of
the window already covered, or, with `prior_weights`, as the prior
probability that a DLT would have occurred by then. Both articles
evaluate their designs with equal weights, the default.

The titration phase and the `stay_on_1_of_3` option of
[`boin_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_simulate.md)
are not available, as neither article defines them for the time-to-event
design.

## References

Yuan, Y., Lin, R., Li, D., Nie, L. and Warren, K. E. (2018).
Time-to-Event Bayesian Optimal Interval Design to Accelerate Phase I
Trials. Clinical Cancer Research, 24(20), 4921-4930.

Lin, R. and Yuan, Y. (2020). Time-to-Event Model-Assisted Designs for
Dose-Finding Trials with Delayed Toxicity. Biostatistics, 21(4),
807-824.

## See also

[`tite_boin_decision_table`](https://gosukehommaex.github.io/simFastBOIN/reference/tite_boin_decision_table.md),
[`boin_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_simulate.md)

## Examples

``` r
# A three month window with two patients a month
trials <- tite_boin_simulate(
  target = 0.30,
  p_true = c(0.05, 0.15, 0.30, 0.45, 0.60),
  n_cohort = 10,
  cohort_size = 3,
  window = 3,
  accrual_rate = 2,
  n_trials = 200,
  seed = 123
)
trials
#> Simulated TITE-BOIN trials
#>   method          : single mean imputation (Yuan et al., 2018)
#>   trials          : 200
#>   doses           : 5
#>   target DLT rate : 0.3
#>   cohorts         : 10 of size 3
#>   max sample size : 30
#>   window          : 3
#>   accrual         : 2 patients per unit of time (exponential)
#> 
#> Average per trial:
#>                   DL1   DL2  DL3  DL4  DL5
#> Patients treated 5.80 10.61 9.12 3.24 0.44
#> DLTs observed    0.32  1.62 2.80 1.46 0.32
#> 
#> Trial duration        : mean 23.55, median 23.44
#> Trials with suspension: 100%
#> Time suspended        : mean 6.19 per trial
#> 
#> Stopping reason (%):
#> 
#> max_sample_size     n_earlystop 
#>            84.5            15.5 

summary(trials$duration)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>   14.66   21.40   23.44   23.55   25.93   32.79 
```
