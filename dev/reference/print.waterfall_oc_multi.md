# Print Operating Characteristics of the Waterfall Design Across Scenarios

Print the design and the cross-scenario summary table of a
`waterfall_oc_multi` object.

## Usage

``` r
# S3 method for class 'waterfall_oc_multi'
print(x, digits = 1, kable = FALSE, kable_format = "pipe", ...)
```

## Arguments

- x:

  An object of class `waterfall_oc_multi`, as returned by
  [`sim_waterfall_multi`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_waterfall_multi.md).

- digits:

  Integer scalar. Number of decimal places. Defaults to 1.

- kable:

  Logical scalar. Print the table with
  [`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html) instead
  of the plain console layout. Defaults to `FALSE`.

- kable_format:

  Character scalar. Format passed to
  [`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html). Defaults
  to `"pipe"`.

- ...:

  Not used.

## Value

`x`, invisibly.

## Examples

``` r
scenarios <- list(
  "Contour A1B3, A2B2" = matrix(c(0.03, 0.10, 0.28,
                                  0.10, 0.30, 0.50), nrow = 2, byrow = TRUE),
  "Contour A1B1" = matrix(c(0.30, 0.40, 0.50,
                            0.42, 0.49, 0.55), nrow = 2, byrow = TRUE)
)

oc <- sim_waterfall_multi(
  target = 0.30,
  scenarios = scenarios,
  cohort_size = 3,
  n_trials = 200,
  seed = 123
)

print(oc)
#> Waterfall drug-combination operating characteristics across 2 scenarios
#>   target DLT rate : 0.3
#>   dose levels     : 2 of drug A, 3 of drug B
#>   trials each     : 200
#>   max sample size : 27
#> 
#>            Scenario MTD contour Contour selected (%) Patients on contour (%)
#>  Contour A1B3, A2B2  A1B3, A2B2                 50.0                    50.5
#>        Contour A1B1        A1B1                 55.5                    52.0
#>  Patients above contour (%) No MTD (%) Patients DLTs
#>                         8.8          0     26.9  6.0
#>                        48.0         11     18.2  6.6
```
