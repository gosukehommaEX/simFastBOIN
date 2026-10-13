# Select the MTD from Completed Trials

Choose the maximum tolerated dose for each of a set of completed trials
from the final patient and DLT counts.

## Usage

``` r
boin_select_mtd(
  n_pts,
  n_tox,
  target,
  cutoff_eli = 0.95,
  extrasafe = FALSE,
  offset = 0.05,
  bound_mtd = FALSE,
  p_tox = NULL,
  mtd_max_estimate = NULL,
  min_mtd_sample = 1
)
```

## Arguments

- n_pts:

  Integer matrix with one row per trial and one column per dose, or a
  vector for a single trial, giving the number of patients treated.

- n_tox:

  Integer matrix or vector of the same shape as `n_pts`, giving the
  number of DLTs observed.

- target:

  Numeric scalar. Target DLT probability.

- cutoff_eli:

  Numeric scalar. Posterior probability cutoff for dose elimination.
  Defaults to 0.95.

- extrasafe:

  Logical scalar. Apply the stricter safety rule at the lowest dose, in
  which case no dose is selected when the lowest dose is judged too
  toxic. Defaults to `FALSE`.

- offset:

  Numeric scalar between 0 and 0.5. Amount by which `cutoff_eli` is
  relaxed for the safety rule. Defaults to 0.05.

- bound_mtd:

  Logical scalar. Require the isotonic estimate at the selected dose to
  be at or below the de-escalation boundary. Defaults to `FALSE`.

- p_tox:

  Numeric scalar. Lowest DLT probability deemed overly toxic, used to
  derive the de-escalation boundary when `bound_mtd` is `TRUE`. Defaults
  to `1.4 * target`.

- mtd_max_estimate:

  Numeric scalar or `NULL`. Largest isotonic estimate a dose may have
  and still be selected. Supplying it bounds the selection whatever
  `bound_mtd` says, and replaces the de-escalation boundary that
  `bound_mtd` would otherwise use. Defaults to `NULL`.

- min_mtd_sample:

  Integer scalar. Smallest number of patients a dose must have received
  to be eligible. Defaults to 1, which admits every dose that treated a
  patient.

## Value

A data frame with one row per trial and columns

- trial:

  Trial index.

- mtd:

  Selected dose level, or `NA` when no dose is selected.

- reason:

  Why the dose was or was not selected.

## Details

Dose elimination is re-derived from the final data of each trial rather
than carried over from the dose-finding stage, which is what the
reference implementation does and which matters when a trial ends at its
maximum sample size. The isotonic fit is computed over the admissible
doses only, so that eliminated doses do not influence the estimates of
the doses that remain. The dose whose estimate is closest to `target` is
selected, with ties broken by a small increasing perturbation across the
admissible doses.

`bound_mtd` caps the estimate at the de-escalation boundary, which
always lies above `target`. A cap at or below the target rate therefore
cannot be expressed that way, and needs `mtd_max_estimate`. The two caps
are independent of the decision table, so bounding the selection does
not change how the trial was run.

Values of `reason` are `"selected"`, `"lowest_dose_eliminated"`,
`"no_admissible_dose"` and `"no_dose_below_bound"`.

## References

Liu S. and Yuan, Y. (2015). Bayesian Optimal Interval Designs for Phase
I Clinical Trials. Journal of the Royal Statistical Society: Series C,
64, 507-523.

## See also

[`boin_isotonic`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_isotonic.md),
[`boin_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_simulate.md)

## Examples

``` r
n_pts <- matrix(c(3, 6, 9, 3,
                  3, 6, 9, 3), nrow = 2, byrow = TRUE)
n_tox <- matrix(c(0, 1, 3, 2,
                  0, 1, 2, 1), nrow = 2, byrow = TRUE)

boin_select_mtd(n_pts, n_tox, target = 0.30)
#>   trial mtd   reason
#> 1     1   3 selected
#> 2     2   4 selected

boin_select_mtd(n_pts, n_tox, target = 0.30, bound_mtd = TRUE)
#>   trial mtd   reason
#> 1     1   3 selected
#> 2     2   4 selected

# Require the isotonic estimate at the selected dose to be at most 0.30,
# which the de-escalation boundary cannot express
boin_select_mtd(n_pts, n_tox, target = 0.30, mtd_max_estimate = 0.30)
#>   trial mtd   reason
#> 1     1   2 selected
#> 2     2   3 selected
```
