# Isotonic Estimates of the DLT Probability

Estimate the DLT probability at each dose under the constraint that
toxicity does not decrease with dose, using the pool adjacent violators
algorithm.

## Usage

``` r
boin_isotonic(n_pts, n_tox, admissible = NULL)
```

## Arguments

- n_pts:

  Integer matrix with one row per trial and one column per dose, or a
  vector for a single trial, giving the number of patients treated.

- n_tox:

  Integer matrix or vector of the same shape as `n_pts`, giving the
  number of DLTs observed.

- admissible:

  Logical matrix of the same shape as `n_pts`, or `NULL`. Doses that are
  `FALSE` are left out of the estimation and returned as `NA`. Defaults
  to all doses that have treated at least one patient.

## Value

A numeric matrix with one row per trial and one column per dose. Doses
that did not enter the estimation are `NA`.

## Details

Pseudo-counts of 0.05 DLTs and 0.1 patients are added before the fit,
and doses are pooled with inverse variance weights. Doses left out
through `admissible` do not influence the estimates of the remaining
doses, which matters when eliminated doses must be excluded before
selecting the MTD.

## References

Liu S. and Yuan, Y. (2015). Bayesian Optimal Interval Designs for Phase
I Clinical Trials. Journal of the Royal Statistical Society: Series C,
64, 507-523.

## See also

[`boin_select_mtd`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_select_mtd.md)

## Examples

``` r
# A single trial
boin_isotonic(n_pts = c(3, 6, 9, 12), n_tox = c(0, 1, 3, 4))
#>            [,1]      [,2]      [,3]      [,4]
#> [1,] 0.01612903 0.1721311 0.3349084 0.3349084

# Several trials at once
n_pts <- matrix(c(3, 6, 9, 12,
                  3, 6, 9, 12,
                  3, 6, 9, 12), nrow = 3, byrow = TRUE)
n_tox <- matrix(c(0, 1, 3, 4,
                  0, 0, 2, 3,
                  1, 2, 4, 6), nrow = 3, byrow = TRUE)
boin_isotonic(n_pts, n_tox)
#>            [,1]       [,2]      [,3]      [,4]
#> [1,] 0.01612903 0.17213115 0.3349084 0.3349084
#> [2,] 0.01000763 0.01000763 0.2252747 0.2520661
#> [3,] 0.33703115 0.33703115 0.4450549 0.5000000
```
