# Print Operating Characteristics of the BOIN Drug-Combination Design Across Scenarios

Display the table produced by
[`sim_comb_boin_multi`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_comb_boin_multi.md),
one row per scenario. The full results of a scenario are printed by
printing its element of `x$results`.

## Usage

``` r
# S3 method for class 'comb_boin_oc_multi'
print(x, digits = 1, kable = FALSE, kable_format = "pipe", ...)
```

## Arguments

- x:

  An object of class `comb_boin_oc_multi`.

- digits:

  Integer scalar. Number of decimal places. Defaults to 1.

- kable:

  Logical scalar. Return the table through
  [`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html) rather
  than printing it directly. Requires the knitr package. Defaults to
  `FALSE`.

- kable_format:

  Character scalar passed to
  [`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html), for
  example `"pipe"`, `"html"` or `"latex"`. Defaults to `"pipe"`.

- ...:

  Further arguments, currently ignored.

## Value

The object `x`, invisibly.

## Examples

``` r
oc <- sim_comb_boin_multi(
  target = 0.30,
  scenarios = list(
    "MTD in the middle" = matrix(c(0.10, 0.20, 0.30,
                                   0.20, 0.30, 0.45), 2, byrow = TRUE),
    "All toxic"         = matrix(c(0.40, 0.50, 0.60,
                                   0.50, 0.60, 0.70), 2, byrow = TRUE)
  ),
  n_cohort = 10,
  cohort_size = 3,
  n_trials = 200,
  seed = 123
)

print(oc)
#> BOIN drug-combination operating characteristics across 2 scenarios
#>   target DLT rate : 0.3
#>   dose levels     : 2 of drug A, 3 of drug B
#>   trials each     : 200
#>   max sample size : 30
#>   overdose cutoff : 0.3
#> 
#>           Scenario   True MTD MTD selected (%) Patients at MTD (%) No MTD (%)
#>  MTD in the middle A2B2, A1B3               59                35.9        0.0
#>          All toxic       A1B1               53                76.4       39.5
#>  Above cutoff selected (%) Patients above cutoff (%) Patients DLTs
#>                       11.0                      12.1     29.1  7.4
#>                       60.5                     100.0     18.0  7.7
```
