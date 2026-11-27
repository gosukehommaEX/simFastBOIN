# Select the MTD Contour from Completed Waterfall Trials

Choose the maximum tolerated dose at every dose level of drug A for each
of a set of completed trials of the waterfall design, from the final
patient and DLT counts and the combinations excluded during the trial.

## Usage

``` r
waterfall_select_mtd(n_pts, n_tox, eliminated, target)
```

## Arguments

- n_pts:

  Integer matrix giving the number of patients treated at each
  combination of one trial, with one row per dose level of drug A and
  one column per dose level of drug B, or an integer array with
  dimensions trial by level of drug A by level of drug B for several
  trials, as returned by
  [`waterfall_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/waterfall_simulate.md).

- n_tox:

  Integer matrix or array of the same shape as `n_pts`, giving the
  number of DLTs observed.

- eliminated:

  Logical matrix or array of the same shape as `n_pts`: whether the
  combination is excluded from the selection, as returned by
  [`waterfall_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/waterfall_simulate.md).

- target:

  Numeric scalar. Target DLT probability.

## Value

A data frame with one row per trial and dose level of drug A and columns

- trial:

  Trial index.

- dose_a:

  Dose level of drug A.

- dose_b:

  Dose level of drug B of the MTD at that level of drug A, or `NA` when
  none is selected.

## Details

The selection follows
[`BOIN::get.oc.comb()`](https://rdrr.io/pkg/BOIN/man/get.oc.comb.html)
with `mtd.contour = TRUE`. The estimates `(y + 0.05) / (n + 0.1)`, set
to 1.1 at the excluded combinations, are smoothed by bivariate isotonic
regression with weights `n + 0.1`, so that they increase along every row
and column, and `1e-5` times the sum of the two dose levels is added.
From the highest level of drug A to the lowest, the treated combination
of the row that is not excluded and whose estimate is closest to
`target` is selected. A row whose first combination is excluded, or
without a treated combination that is not excluded, has no MTD. When the
MTD of a row is not at a higher level of drug B than the MTD of the row
above it, it is moved to the level of drug B of that MTD, so that the
contour does not turn back; the combination it is moved to need not have
been treated.

The isotonic regression is
[`Iso::biviso()`](https://rdrr.io/pkg/Iso/man/biviso.html), the routine
the reference implementation calls.

[`BOIN::select.mtd.comb()`](https://rdrr.io/pkg/BOIN/man/select.mtd.comb.html)
with `mtd.contour = TRUE` ends with a different rule (a row takes the
last level of drug B when the row above has its MTD there, and has no
MTD when both are there); the operating characteristics of
[`BOIN::get.oc.comb()`](https://rdrr.io/pkg/BOIN/man/get.oc.comb.html)
come from the rule described here.

## References

Zhang, L. and Yuan, Y. (2016). A Practical Bayesian Design to Identify
the Maximum Tolerated Dose Contour for Drug Combination Trials.
Statistics in Medicine, 35(27), 4924-4936.

## See also

[`waterfall_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/waterfall_simulate.md),
[`sim_waterfall`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_waterfall.md)

## Examples

``` r
n_pts <- matrix(c(3, 3, 9,
                  6, 9, 3), nrow = 2, byrow = TRUE)
n_tox <- matrix(c(0, 0, 3,
                  1, 3, 2), nrow = 2, byrow = TRUE)
eliminated <- matrix(FALSE, nrow = 2, ncol = 3)
eliminated[2, 3] <- TRUE

waterfall_select_mtd(n_pts, n_tox, eliminated, target = 0.30)
#>   trial dose_a dose_b
#> 1     1      1      3
#> 2     1      2      2
```
