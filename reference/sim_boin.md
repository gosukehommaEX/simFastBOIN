# Operating Characteristics of a BOIN Design

Simulate a BOIN dose-finding trial many times under one dose-toxicity
scenario and summarize how often each dose is selected as the MTD, how
many patients are treated at each dose and how many DLTs are observed.

## Usage

``` r
sim_boin(
  target,
  p_true,
  n_cohort,
  cohort_size,
  n_trials = 10000,
  start_dose = 1,
  n_earlystop = 18,
  p_saf = NULL,
  p_tox = NULL,
  cutoff_eli = 0.95,
  extrasafe = FALSE,
  offset = 0.05,
  titration = FALSE,
  stay_on_1_of_3 = FALSE,
  bound_mtd = FALSE,
  mtd_max_estimate = NULL,
  min_mtd_sample = 1,
  overdose_cutoff = NULL,
  n_earlystop_rule = c("with_stay", "simple"),
  keep_trials = FALSE,
  verbose = FALSE,
  seed = 123
)
```

## Arguments

- target:

  Numeric scalar. Target DLT probability, for example 0.30.

- p_true:

  Numeric vector. True DLT probability at each dose level, in increasing
  dose order.

- n_cohort:

  Integer scalar. Number of cohorts in a trial.

- cohort_size:

  Integer scalar or vector. Number of patients per cohort. A scalar is
  used for every cohort. A vector shorter than `n_cohort` is padded with
  its last element and a longer one is truncated.

- n_trials:

  Integer scalar. Number of trials to simulate. Defaults to 10000.

- start_dose:

  Integer scalar. Dose level for the first cohort. Defaults to 1. It is
  ignored when `titration` is `TRUE`, because the titration phase always
  begins at the lowest dose.

- n_earlystop:

  Integer scalar. The trial stops once this many patients have been
  treated at the current dose and the design would stay there. Defaults
  to 18. Set it to a value above the maximum sample size to switch this
  rule off.

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

- titration:

  Logical scalar. Start with single patient cohorts until the first DLT
  is seen. Ignored when the first cohort size is one. Defaults to
  `FALSE`.

- stay_on_1_of_3:

  Logical scalar. When `TRUE`, one DLT out of three patients leads to
  staying at the current dose rather than de-escalating. Defaults to
  `FALSE`. See
  [`boin_boundary`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_boundary.md).

- bound_mtd:

  Logical scalar. Require the isotonic estimate at the selected dose to
  be at or below the de-escalation boundary. Defaults to `FALSE`.

- mtd_max_estimate:

  Numeric scalar or `NULL`. Largest isotonic estimate a dose may have
  and still be selected as the MTD. Supplying it bounds the selection
  whatever `bound_mtd` says, and unlike `bound_mtd` it can be set at or
  below the target rate. It changes only the selection, never the
  dose-finding itself. Defaults to `NULL`.

- min_mtd_sample:

  Integer scalar. Smallest number of patients a dose must have received
  to be eligible as the MTD. Defaults to 1.

- overdose_cutoff:

  Numeric scalar or `NULL`. Doses whose true DLT probability exceeds
  this value count as overdoses in the `overdose` component of the
  result. Defaults to `NULL`, which uses `target`.

- n_earlystop_rule:

  Character scalar, either `"with_stay"` or `"simple"`. See
  [`boin_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_simulate.md).

- keep_trials:

  Logical scalar. Keep the full trial by trial data in the result.
  Defaults to `FALSE`.

- verbose:

  Logical scalar. Report progress while the simulation runs. Defaults to
  `FALSE`.

- seed:

  Integer scalar or `NULL`. Random seed. The state of the calling
  session is restored on exit. Defaults to 123.

## Value

An object of class `boin_oc`, which is a list with components

- sel_percent:

  Percentage of trials selecting each dose as the MTD.

- percent_no_mtd:

  Percentage of trials ending without an MTD.

- n_pts_dose:

  Average number of patients treated at each dose.

- n_tox_dose:

  Average number of DLTs observed at each dose.

- total_n_pts:

  Average total number of patients per trial.

- total_n_tox:

  Average total number of DLTs per trial.

- overdose:

  List describing the design's relationship with doses above
  `overdose_cutoff`. Two questions are answered separately. How far were
  patients exposed during the trial: `pct_patients` (the percentage of
  all simulated patients treated at such a dose, that is the probability
  that a patient is dosed above the cutoff), `pct_patients_by_trial`
  (the same percentage computed within each trial and then averaged),
  `avg_n_patients`, `pct_trials_any`, `pct_trials_over_60` and
  `pct_trials_over_80`. And how often did the design end up recommending
  such a dose: `pct_trials_mtd_above` (the percentage of all trials
  whose selected MTD is above the cutoff) and
  `pct_mtd_above_when_selected` (the same among the trials that selected
  an MTD at all). The `cutoff` and the dose levels `doses` that exceed
  it are also returned.

- stop_reason_percent:

  Percentage of trials by reason for stopping.

- trials:

  Trial level data when `keep_trials` is `TRUE`, otherwise `NULL`.

together with the design parameters and the call.

## Details

The engine draws one uniform variate per patient, in enrollment order,
and applies the decision rules in the order used by
[`BOIN::get.oc()`](https://rdrr.io/pkg/BOIN/man/get.oc.html). See
[`boin_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_simulate.md)
for the two places where a variate is drawn but not used. With the same
seed and matching arguments the two implementations agree trial by
trial, not merely on average.

Note that `n_earlystop` defaults to 18 here, whereas the reference
implementation defaults to 100, which in practice switches the rule off.

## References

Liu S. and Yuan, Y. (2015). Bayesian Optimal Interval Designs for Phase
I Clinical Trials. Journal of the Royal Statistical Society: Series C,
64, 507-523.

## See also

[`sim_boin_multi`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin_multi.md),
[`boin_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_simulate.md)

## Examples

``` r
oc <- sim_boin(
  target = 0.30,
  p_true = c(0.05, 0.15, 0.30, 0.45, 0.60),
  n_cohort = 10,
  cohort_size = 3,
  n_trials = 500,
  seed = 123
)
oc
#> BOIN operating characteristics
#>   target DLT rate : 0.3
#>   trials          : 500
#>   cohorts         : 10 of size 3
#>   max sample size : 30
#>   lambda_e / _d   : 0.2365 / 0.3585
#> 
#>                   DL1  DL2  DL3  DL4  DL5 Total / No MTD
#> True DLT rate (%) 5.0 15.0 30.0 45.0 60.0               
#> MTD selected (%)  0.8 24.8 54.2 19.0  1.2            0.0
#> Patients treated  4.0  8.9 10.9  4.6  0.7           29.1
#> Patients with DLT 0.2  1.4  3.2  2.1  0.4            7.3
#> 
#> Doses above a true DLT rate of 0.3: DL4, DL5
#>   Patients treated there              : 18.1%
#>   Trials treating any patient there   : 56%
#>   Trials treating over 60% there      : 1.4%
#>   Trials treating over 80% there      : 0%
#>   Trials selecting an MTD there       : 20.2%

# \donttest{
# A larger run with the safety options switched on
oc_safe <- sim_boin(
  target = 0.30,
  p_true = c(0.30, 0.40, 0.50, 0.60, 0.70),
  n_cohort = 20,
  cohort_size = 3,
  n_trials = 10000,
  extrasafe = TRUE,
  bound_mtd = TRUE,
  titration = TRUE,
  seed = 123
)
oc_safe
#> BOIN operating characteristics
#>   target DLT rate : 0.3
#>   trials          : 10000
#>   cohorts         : 20 of size 3
#>   max sample size : 60
#>   lambda_e / _d   : 0.2365 / 0.3585
#> 
#>                    DL1  DL2  DL3  DL4  DL5 Total / No MTD
#> True DLT rate (%) 30.0 40.0 50.0 60.0 70.0               
#> MTD selected (%)  40.1 20.8  2.7  0.1  0.0           36.2
#> Patients treated  12.3  9.8  3.9  1.2  0.3           27.5
#> Patients with DLT  3.7  3.9  2.0  0.7  0.2           10.5
#> 
#> Doses above a true DLT rate of 0.3: DL2, DL3, DL4, DL5
#>   Patients treated there              : 55.4%
#>   Trials treating any patient there   : 78.7%
#>   Trials treating over 60% there      : 7.6%
#>   Trials treating over 80% there      : 1.2%
#>   Trials selecting an MTD there       : 23.7%

# How often are patients dosed above a true DLT rate of 0.33?
oc_cut <- sim_boin(
  target = 0.30,
  p_true = c(0.10, 0.20, 0.30, 0.42, 0.55),
  n_cohort = 20,
  cohort_size = 3,
  n_trials = 10000,
  overdose_cutoff = 0.33,
  seed = 123
)
oc_cut$overdose$pct_patients            # patients dosed above 0.33
#> [1] 21.38077
oc_cut$overdose$pct_trials_mtd_above     # trials recommending a dose above 0.33
#> [1] 15.89
# }
```
