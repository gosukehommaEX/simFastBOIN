# Print the Operating Characteristics of the BOIN-ETC Design Across Scenarios

Print the Operating Characteristics of the BOIN-ETC Design Across
Scenarios

## Usage

``` r
# S3 method for class 'boin_etc_oc_multi'
print(x, digits = 1, kable = FALSE, kable_format = "pipe", ...)
```

## Arguments

- x:

  An object of class `boin_etc_oc_multi` from
  [`sim_boin_etc_multi`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_boin_etc_multi.md).

- digits:

  Integer scalar. Decimal places. Defaults to 1.

- kable:

  Logical scalar. Print the table with
  [`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html). Defaults
  to `FALSE`.

- kable_format:

  Character scalar. Format passed to
  [`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html). Defaults
  to `"pipe"`.

- ...:

  Not used.

## Value

`x`, invisibly.

## See also

[`sim_boin_etc_multi`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_boin_etc_multi.md)

## Examples

``` r
scenarios <- list(
  "Scenario 1" = list(
    p_true = matrix(c(0.05, 0.10, 0.10, 0.30, 0.30, 0.50), nrow = 2),
    p_resp = matrix(c(0.10, 0.25, 0.25, 0.60, 0.60, 0.65), nrow = 2)),
  "Scenario 2" = list(
    p_true = matrix(c(0.05, 0.15, 0.15, 0.35, 0.30, 0.55), nrow = 2),
    p_resp = matrix(c(0.25, 0.60, 0.60, 0.65, 0.65, 0.70), nrow = 2))
)
res <- sim_boin_etc_multi(target = 0.30, target_resp = 0.60,
                          scenarios = scenarios, n_cohort = 8,
                          cohort_size = 3, n_trials = 100)
print(res)
#> BOIN-ETC1 operating characteristics across 2 scenarios
#>   target DLT rate      : 0.3
#>   target response rate : 0.6
#>   dose levels          : 2 of drug A, 3 of drug B
#>   trials each          : 100
#>   max sample size      : 24
#>   subtrial end (n_stop): 6
#>   start of a row       : diagonal
#>   rule 4               : paper
#>   utility scores       : 100, 40, 60, 0; ties lowest
#> 
#>    Scenario   True OBD OBD selected (%) Patients at OBD (%) No OBD (%)
#>  Scenario 1 A2B2, A1B3               74                43.0          0
#>  Scenario 2 A2B1, A1B2               84                61.9          0
#>  Above limit selected (%) Patients above limit (%) Patients DLTs Responses
#>                         1                      5.0     22.6  4.7       8.9
#>                         6                      9.3     17.3  2.9       9.3
```
