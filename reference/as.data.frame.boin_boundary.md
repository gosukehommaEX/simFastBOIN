# Coerce BOIN Decision Boundaries to a Data Frame

Return the integer decision boundaries produced by
[`boin_boundary`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_boundary.md)
as one row per sample size.

## Usage

``` r
# S3 method for class 'boin_boundary'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x:

  An object of class `boin_boundary`.

- row.names:

  Passed on for compatibility with the generic, normally unused.

- optional:

  Passed on for compatibility with the generic, normally unused.

- ...:

  Further arguments, currently ignored.

## Value

A data frame with columns `n_pts`, `escalate_if_dlt_leq`,
`deescalate_if_dlt_geq` and `eliminate_if_dlt_geq`, plus
`stop_if_dlt_geq` when the boundaries were built with
`extrasafe = TRUE`.

## Examples

``` r
bd <- boin_boundary(target = 0.30, max_n = 18)
as.data.frame(bd)[seq(3, 18, by = 3), ]
#>    n_pts escalate_if_dlt_leq deescalate_if_dlt_geq eliminate_if_dlt_geq
#> 3      3                   0                     2                    3
#> 6      6                   1                     3                    4
#> 9      9                   2                     4                    5
#> 12    12                   2                     5                    7
#> 15    15                   3                     6                    8
#> 18    18                   4                     7                    9
```
