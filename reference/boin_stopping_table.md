# Safety Stopping Boundary as a Two-Row Table

Present the safety stopping boundary at the lowest dose as a compact
table with one column per sample size, which is the shape usually wanted
in a protocol or a report.

## Usage

``` r
boin_stopping_table(x, cohort_size = NULL)
```

## Arguments

- x:

  An object of class `boin_boundary`, built with `extrasafe = TRUE`.

- cohort_size:

  Integer scalar or `NULL`. When supplied, only the sample sizes that
  are multiples of `cohort_size` are kept.

## Value

A data frame with two rows, labeled by the quantity they hold, and one
column per sample size. Sample sizes at which no number of DLTs triggers
the rule are dropped.

## Details

The same numbers are available one row per sample size from
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) on the
boundary object. This function transposes them and labels the rows,
which suits a table pasted into a document or rendered with
[`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html) or
`DT::datatable()`.

## See also

[`boin_boundary`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_boundary.md)

## Examples

``` r
bd <- boin_boundary(target = 0.30, max_n = 18, extrasafe = TRUE)

boin_stopping_table(bd)
#>                                                               3 4 5 6 7 8 9 10
#> Number of evaluable patients treated at the lowest dose level 3 4 5 6 7 8 9 10
#> Stop the trial if # of DLT >=                                 2 3 3 4 4 4 5  5
#>                                                               11 12 13 14 15 16
#> Number of evaluable patients treated at the lowest dose level 11 12 13 14 15 16
#> Stop the trial if # of DLT >=                                  6  6  6  7  7  8
#>                                                               17 18
#> Number of evaluable patients treated at the lowest dose level 17 18
#> Stop the trial if # of DLT >=                                  8  8

# Only the sample sizes reached at the end of a cohort of three
boin_stopping_table(bd, cohort_size = 3)
#>                                                               3 6 9 12 15 18
#> Number of evaluable patients treated at the lowest dose level 3 6 9 12 15 18
#> Stop the trial if # of DLT >=                                 2 4 5  6  7  8
```
