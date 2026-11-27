# Select the MTD of a gBOIN Trial

Select the maximum tolerated dose from the final data of one or more
trials of the gBOIN design of Mu et al. (2019), in which the toxicity of
a dose is the total of the normalized scores of its patients.

## Usage

``` r
gboin_select_mtd(
  n_pts,
  score_sum,
  target,
  cutoff_eli = 0.95,
  extrasafe = FALSE,
  offset = 0.05,
  min_mtd_sample = 1
)
```

## Arguments

- n_pts:

  Integer vector or matrix. Number of patients treated at each dose. A
  matrix holds one trial per row.

- score_sum:

  Numeric vector or matrix of the same shape. Total of the normalized
  toxicity scores, each between 0 and 1, of the patients treated at each
  dose.

- target:

  Numeric scalar. Target value of the mean normalized toxicity score.

- cutoff_eli:

  Numeric scalar. Posterior probability cutoff for dose elimination.
  Defaults to 0.95.

- extrasafe:

  Logical scalar. Apply the stricter safety stopping rule at the lowest
  dose. Defaults to `FALSE`.

- offset:

  Numeric scalar between 0 and 0.5. Amount by which `cutoff_eli` is
  relaxed for the safety stopping rule. Defaults to 0.05.

- min_mtd_sample:

  Integer scalar. Smallest number of patients a dose must have been
  given to be selected. Defaults to 1.

## Value

A data frame with one row per trial and the columns `trial`, `mtd` (the
selected dose, `NA` when none is selected) and `reason` (`"selected"`,
`"lowest_dose_eliminated"` or `"no_admissible_dose"`).

## Details

The selection is that of
[`boin_select_mtd`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_select_mtd.md)
with the total score in the place of the number of DLTs: doses whose
posterior probability of a mean score above the target exceeds
`cutoff_eli` are excluded together with the doses above them, the mean
scores of the remaining doses that were given to patients are made
monotone by isotonic regression, and the dose whose estimate is closest
to the target is selected. This is the selection of the
'UnifiedDoseFinding' and 'TITEgBOIN' packages. With whole numbers of
DLTs as `score_sum` the result is that of
[`boin_select_mtd`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_select_mtd.md).

## References

Mu, R., Yuan, Y., Xu, J., Mandrekar, S. J. and Yin, J. (2019). gBOIN: a
unified model-assisted phase I trial design accounting for toxicity
grades, and binary or continuous end points. Journal of the Royal
Statistical Society: Series C, 68, 289-308.
[doi:10.1111/rssc.12263](https://doi.org/10.1111/rssc.12263)

## See also

[`sim_gboin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_gboin.md),
[`boin_select_mtd`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_select_mtd.md)

## Examples

``` r
# Equivalent toxicity scores divided by their maximum of 1.5
n_pts <- c(3, 6, 12, 9, 0, 0)
score_sum <- c(0.5, 1.5, 5.5, 6.0, 0, 0) / 1.5
gboin_select_mtd(n_pts, score_sum, target = 0.47 / 1.5)
#>   trial mtd   reason
#> 1     1   3 selected
```
