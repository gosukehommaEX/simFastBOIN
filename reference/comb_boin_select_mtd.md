# Select the MTD Combination from Completed Drug-Combination Trials

Choose the maximum tolerated dose combination for each of a set of
completed trials of the BOIN drug-combination design from the final
patient and DLT counts.

## Usage

``` r
comb_boin_select_mtd(
  n_pts,
  n_tox,
  target,
  cutoff_eli = 0.95,
  extrasafe = FALSE,
  offset = 0.05,
  bound_mtd = FALSE,
  p_tox = NULL
)
```

## Arguments

- n_pts:

  Integer matrix giving the number of patients treated at each
  combination of one trial, with one row per dose level of drug A and
  one column per dose level of drug B, or an integer array with
  dimensions trial by level of drug A by level of drug B for several
  trials, as returned by
  [`comb_boin_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/comb_boin_simulate.md).

- n_tox:

  Integer matrix or array of the same shape as `n_pts`, giving the
  number of DLTs observed.

- target:

  Numeric scalar. Target DLT probability.

- cutoff_eli:

  Numeric scalar. Posterior probability cutoff for eliminating a
  combination. Defaults to 0.95.

- extrasafe:

  Logical scalar. Apply the stricter safety rule at the lowest
  combination, in which case no combination is selected when it is
  judged too toxic. Defaults to `FALSE`.

- offset:

  Numeric scalar between 0 and 0.5. Amount by which `cutoff_eli` is
  relaxed for the safety rule. Defaults to 0.05.

- bound_mtd:

  Logical scalar. Require the isotonic estimate at the selected
  combination to be below the de-escalation boundary. Defaults to
  `FALSE`.

- p_tox:

  Numeric scalar. Lowest DLT probability deemed overly toxic, used to
  derive the de-escalation boundary when `bound_mtd` is `TRUE`. Defaults
  to `1.4 * target`.

## Value

A data frame with one row per trial and columns

- trial:

  Trial index.

- mtd_a:

  Dose level of drug A in the selected combination, or `NA` when none is
  selected.

- mtd_b:

  Dose level of drug B in the selected combination, or `NA`.

- reason:

  Why a combination was or was not selected.

## Details

The selection follows the reference implementation used by
[`BOIN::get.oc.comb()`](https://rdrr.io/pkg/BOIN/man/get.oc.comb.html).
Elimination is re-derived from the final data: within each row, the
first combination with at least three patients whose posterior
probability of a DLT rate above `target` exceeds `cutoff_eli` is
eliminated together with the rest of its row and the rest of its column.
If the lowest combination is eliminated, no combination is selected.
Otherwise the estimates `(y + 0.05) / (n + 0.1)` are smoothed by
bivariate isotonic regression with weights `n + 0.1`, so that they
increase along every row and column, and the treated, non-eliminated
combination whose estimate is closest to `target` is selected. A small
amount, `1e-5` times the sum of the two dose levels, is added to every
estimate first, so that ties go to the lower combination.

The isotonic regression is
[`Iso::biviso()`](https://rdrr.io/pkg/Iso/man/biviso.html), the routine
the reference implementation calls. Combinations that the regression
pools together have equal estimates in exact arithmetic but can differ
in the last digits, and with the small amount above, these digits decide
between combinations on the same anti-diagonal. Another algorithm for
the same regression would therefore select differently in a few trials
out of a thousand.

The article of Lin and Yin (2017) estimates the DLT rates as `y / n` and
resolves ties differently; the package follows the reference
implementation so that the two can be compared trial by trial.

Values of `reason` are `"selected"`, `"lowest_dose_eliminated"` and
`"no_dose_below_bound"`.

## References

Lin, R. and Yin, G. (2017). Bayesian Optimal Interval Design for Dose
Finding in Drug-Combination Trials. Statistical Methods in Medical
Research, 26(5), 2155-2167.

## See also

[`comb_boin_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/comb_boin_simulate.md),
[`sim_comb_boin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_comb_boin.md)

## Examples

``` r
n_pts <- matrix(c(3, 6, 3,
                  3, 3, 0), nrow = 2, byrow = TRUE)
n_tox <- matrix(c(0, 1, 1,
                  0, 2, 0), nrow = 2, byrow = TRUE)

comb_boin_select_mtd(n_pts, n_tox, target = 0.30)
#>   trial mtd_a mtd_b   reason
#> 1     1     1     3 selected
```
