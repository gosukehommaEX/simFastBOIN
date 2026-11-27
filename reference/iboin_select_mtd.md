# Select the MTD of an iBOIN Trial

Select the maximum tolerated dose from the final data of one or more
trials of the iBOIN design of Zhou et al. (2021), with or without the
prior information in the estimates of the DLT probabilities.

## Usage

``` r
iboin_select_mtd(
  n_pts,
  n_tox,
  target,
  skeleton,
  prior_ess,
  with_prior = TRUE,
  cutoff_eli = 0.95,
  extrasafe = FALSE,
  offset = 0.05,
  robust = FALSE,
  min_mtd_sample = 1
)
```

## Arguments

- n_pts:

  Integer vector or matrix. Number of patients treated at each dose. A
  matrix holds one trial per row.

- n_tox:

  Integer vector or matrix of the same shape. Number of DLTs observed at
  each dose.

- target:

  Numeric scalar. Target DLT probability.

- skeleton:

  Numeric vector, one value per dose. Prior DLT probability at each
  dose.

- prior_ess:

  Whole number or vector of whole numbers, one per dose. Prior effective
  sample size.

- with_prior:

  Logical scalar. When `TRUE`, the default, the prior information is
  added to the data before the isotonic regression. When `FALSE` the
  selection is that of
  [`boin_select_mtd`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_select_mtd.md).

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

  Logical scalar. Use the robust prior, see
  [`iboin_boundary`](https://gosukehommaex.github.io/simFastBOIN/reference/iboin_boundary.md).
  Defaults to `FALSE`.

- min_mtd_sample:

  Integer scalar. Smallest number of patients a dose must have been
  given to be selected. Defaults to 1.

## Value

A data frame with one row per trial and the columns `trial`, `mtd` (the
selected dose, `NA` when none is selected) and `reason` (`"selected"`,
`"lowest_dose_eliminated"` or `"no_admissible_dose"`).

## Details

The article gives two estimates of the DLT probability of a dose for the
selection of the MTD: the observed rate `n_tox / n_pts`, and the rate
with the prior information,
`(n_tox + prior_ess * skeleton) / (n_pts + prior_ess)`. Either is made
monotone by isotonic regression, and the dose whose estimate is closest
to the target is selected. With `with_prior = TRUE` the numbers
`n_tox + prior_ess * skeleton` and `n_pts + prior_ess` take the place of
the DLTs and the patients in the estimate of
[`boin_select_mtd`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_select_mtd.md);
the doses that can be selected are still those that have been given to
patients and are not eliminated, and the elimination uses the data
alone, under a uniform prior.

The article does not say which estimate its simulation study uses. Its
operating characteristics (Table S.1 of its supporting information) are
reproduced with the prior information and not without it, see
`inst/reproduce/reproduce_Zhou2021.R`. The estimate with the prior
selects the true MTD more often when the prior MTD is right and less
often when it is wrong.

## References

Zhou, Y., Lee, J. J., Wang, S., Bailey, S. and Yuan, Y. (2021).
Incorporating historical information to improve phase I clinical trials.
Pharmaceutical Statistics, 20, 1017-1034.
[doi:10.1002/pst.2121](https://doi.org/10.1002/pst.2121)

## See also

[`sim_iboin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_iboin.md),
[`boin_select_mtd`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_select_mtd.md)

## Examples

``` r
skeleton <- c(0.10, 0.19, 0.30, 0.42, 0.54)
n_pts <- c(3, 6, 12, 9, 0)
n_tox <- c(0, 1, 3, 4, 0)

iboin_select_mtd(n_pts, n_tox, target = 0.30, skeleton = skeleton,
                 prior_ess = 3)
#>   trial mtd   reason
#> 1     1   3 selected
iboin_select_mtd(n_pts, n_tox, target = 0.30, skeleton = skeleton,
                 prior_ess = 3, with_prior = FALSE)
#>   trial mtd   reason
#> 1     1   3 selected
```
