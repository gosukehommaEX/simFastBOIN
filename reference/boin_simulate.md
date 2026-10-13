# Simulate BOIN Trials

Run the BOIN dose-finding algorithm over many simulated trials and
return the patient and DLT counts at every dose. This is the simulation
engine behind
[`sim_boin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin.md),
exposed for users who want the raw trial data.

## Usage

``` r
boin_simulate(
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
  n_earlystop_rule = c("with_stay", "simple"),
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

- titration:

  Logical scalar. Start with single patient cohorts until the first DLT
  is seen. Ignored when the first cohort size is one. Defaults to
  `FALSE`.

- stay_on_1_of_3:

  Logical scalar. When `TRUE`, one DLT out of three patients leads to
  staying at the current dose rather than de-escalating. Defaults to
  `FALSE`. See
  [`boin_boundary`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_boundary.md).

- n_earlystop_rule:

  Character scalar, either `"with_stay"` or `"simple"`. Under
  `"with_stay"` the trial stops at `n_earlystop` only when the design
  would remain at the current dose, which is the rule used by the
  reference implementation. Under `"simple"` it stops as soon as the
  threshold is reached.

- seed:

  Integer scalar or `NULL`. Random seed. When `NULL` the current state
  of the random number generator is used. The state of the calling
  session is restored on exit. Defaults to 123.

## Value

An object of class `boin_trials`, which is a list with components

- n_pts:

  Integer matrix, `n_trials` by number of doses, of patients treated.

- n_tox:

  Integer matrix of the same shape, of DLTs observed.

- eliminated:

  Logical matrix of the same shape, marking doses eliminated during the
  trial.

- cohorts_used:

  Integer vector of cohorts completed in each trial.

- stop_reason:

  Character vector giving why each trial ended.

- boundary:

  The `boin_boundary` object used.

- settings:

  List of the design parameters.

## Details

One uniform variate is drawn per patient, in enrollment order. Two
details matter for reproducibility: the titration phase draws one
variate per dose level whether or not it is used, and a cohort cut short
by the maximum sample size still draws a full cohort's worth. This is
the random number consumption of
[`BOIN::get.oc()`](https://rdrr.io/pkg/BOIN/man/get.oc.html), so with
the same seed and matching arguments a run reproduces the reference
implementation trial by trial.

Possible values of `stop_reason` are `"lowest_dose_eliminated"`,
`"lowest_dose_too_toxic"`, `"n_earlystop"`, `"max_sample_size"` and
`"max_cohorts"`.

## References

Liu S. and Yuan, Y. (2015). Bayesian Optimal Interval Designs for Phase
I Clinical Trials. Journal of the Royal Statistical Society: Series C,
64, 507-523.

## See also

[`sim_boin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin.md),
[`boin_select_mtd`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_select_mtd.md)

## Examples

``` r
trials <- boin_simulate(
  target = 0.30,
  p_true = c(0.05, 0.15, 0.30, 0.45, 0.60),
  n_cohort = 10,
  cohort_size = 3,
  n_trials = 200,
  seed = 123
)
trials
#> Simulated BOIN trials
#>   trials          : 200
#>   doses           : 5
#>   target DLT rate : 0.3
#>   cohorts         : 10 of size 3
#>   max sample size : 30
#> 
#> Average per trial:
#>                   DL1  DL2   DL3  DL4  DL5
#> Patients treated 4.18 9.12 10.61 4.53 0.64
#> DLTs observed    0.23 1.40  3.18 2.06 0.40
#> 
#> Stopping reason (%):
#> 
#> max_sample_size     n_earlystop 
#>              79              21 

head(trials$n_pts)
#>      DL1 DL2 DL3 DL4 DL5
#> [1,]   3   6   6  15   0
#> [2,]   3  18   6   0   0
#> [3,]   3   6   9  12   0
#> [4,]   3   6   6   6   9
#> [5,]   3  18   6   0   0
#> [6,]   3  15  12   0   0
table(trials$stop_reason)
#> 
#> max_sample_size     n_earlystop 
#>             158              42 
```
