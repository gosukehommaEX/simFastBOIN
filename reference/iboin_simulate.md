# Simulate iBOIN Trials

Run the iBOIN dose-finding algorithm of Zhou et al. (2021) over many
simulated trials and return the patient and DLT counts at every dose.
This is the simulation engine behind
[`sim_iboin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_iboin.md),
exposed for users who want the raw trial data.

## Usage

``` r
iboin_simulate(
  target,
  p_true,
  skeleton,
  prior_ess,
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
  robust = FALSE,
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

- skeleton:

  Numeric vector, one value per dose. Prior DLT probability at each
  dose, strictly between 0 and 1.

- prior_ess:

  Whole number or vector of whole numbers, one per dose. Prior effective
  sample size. See
  [`iboin_boundary`](https://gosukehommaex.github.io/simFastBOIN/reference/iboin_boundary.md).

- n_cohort:

  Integer scalar. Number of cohorts in a trial.

- cohort_size:

  Integer scalar or vector. Number of patients per cohort. A scalar is
  used for every cohort. A vector shorter than `n_cohort` is padded with
  its last element and a longer one is truncated.

- n_trials:

  Integer scalar. Number of trials to simulate. Defaults to 10000.

- start_dose:

  Integer scalar. Dose level for the first cohort. Defaults to 1.

- n_earlystop:

  Integer scalar. The trial stops once this many patients have been
  treated at the current dose and the design would stay there. Defaults
  to 18. The simulations of the article do not stop early; a value above
  the maximum sample size switches the rule off.

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

- robust:

  Logical scalar. Use the robust prior, which drops the prior
  information above the prior MTD when the prior MTD lies in the upper
  half of the doses. Defaults to `FALSE`. See
  [`iboin_boundary`](https://gosukehommaex.github.io/simFastBOIN/reference/iboin_boundary.md).

- n_earlystop_rule:

  Character scalar, either `"with_stay"` or `"simple"`. Under
  `"with_stay"` the trial stops at `n_earlystop` only when the design
  would remain at the current dose. Under `"simple"` it stops as soon as
  the threshold is reached.

- seed:

  Integer scalar or `NULL`. Random seed. When `NULL` the current state
  of the random number generator is used. The state of the calling
  session is restored on exit. Defaults to 123.

## Value

An object of class `c("iboin_trials", "boin_trials")`, which is a list
with components

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

  The `iboin_boundary` object used.

- design:

  The name of the design, `"iBOIN"`.

- settings:

  List of the design parameters.

## Details

The trial runs as a BOIN trial, see
[`boin_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_simulate.md),
with the boundaries of
[`iboin_boundary`](https://gosukehommaex.github.io/simFastBOIN/reference/iboin_boundary.md),
which depend on the dose and on the number of patients treated there.
Dose elimination and the safety stopping rule use a uniform prior, as in
the BOIN design. The titration of the BOIN design is not available.

One uniform variate is drawn per patient, in enrollment order, and a
cohort cut short by the maximum sample size still draws a full cohort's
worth. With a prior effective sample size of 0 at every dose the trials
are therefore, for the same seed, those of
[`boin_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_simulate.md).

Possible values of `stop_reason` are `"lowest_dose_eliminated"`,
`"lowest_dose_too_toxic"`, `"n_earlystop"`, `"max_sample_size"` and
`"max_cohorts"`.

## References

Zhou, Y., Lee, J. J., Wang, S., Bailey, S. and Yuan, Y. (2021).
Incorporating historical information to improve phase I clinical trials.
Pharmaceutical Statistics, 20, 1017-1034.
[doi:10.1002/pst.2121](https://doi.org/10.1002/pst.2121)

## See also

[`sim_iboin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_iboin.md),
[`iboin_boundary`](https://gosukehommaex.github.io/simFastBOIN/reference/iboin_boundary.md),
[`boin_select_mtd`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_select_mtd.md)

## Examples

``` r
trials <- iboin_simulate(
  target = 0.30,
  p_true = c(0.08, 0.15, 0.31, 0.45, 0.55),
  skeleton = c(0.10, 0.19, 0.30, 0.42, 0.54),
  prior_ess = 3,
  n_cohort = 10,
  cohort_size = 3,
  n_trials = 200,
  seed = 123
)
trials
#> Simulated iBOIN trials
#>   trials          : 200
#>   doses           : 5
#>   target DLT rate : 0.3
#>   cohorts         : 10 of size 3
#>   max sample size : 30
#> 
#> Average per trial:
#>                   DL1  DL2   DL3  DL4  DL5
#> Patients treated 3.52 8.84 12.63 3.86 0.42
#> DLTs observed    0.31 1.30  3.86 1.81 0.22
#> 
#> Stopping reason (%):
#> 
#> max_sample_size     n_earlystop 
#>              82              18 

head(trials$n_pts)
#>      DL1 DL2 DL3 DL4 DL5
#> [1,]   3   6   9   9   3
#> [2,]   3  18   6   0   0
#> [3,]   3   6  15   6   0
#> [4,]   3   6   9   6   6
#> [5,]   3  18   6   0   0
#> [6,]   3  15  12   0   0
table(trials$stop_reason)
#> 
#> max_sample_size     n_earlystop 
#>             164              36 
```
