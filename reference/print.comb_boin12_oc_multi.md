# Print the Operating Characteristics of Comb-BOIN12 Across Scenarios

Print the Operating Characteristics of Comb-BOIN12 Across Scenarios

## Usage

``` r
# S3 method for class 'comb_boin12_oc_multi'
print(x, digits = 1, kable = FALSE, kable_format = "pipe", ...)
```

## Arguments

- x:

  An object of class `comb_boin12_oc_multi` from
  [`sim_comb_boin12_multi`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_comb_boin12_multi.md).

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

## Details

One row per scenario: the true OBD combinations, the percentage of
trials selecting one of them and of the patients treated there, the
percentage of trials selecting no OBD combination, the percentages
selecting and treating combinations above the DLT limit, and the average
numbers of patients, DLTs and responses.

## See also

[`sim_comb_boin12_multi`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_comb_boin12_multi.md)

## Examples

``` r
scenarios <- list(
  "Scenario 16" = list(
    p_true = matrix(c(0.05, 0.10, 0.15, 0.15, 0.20, 0.45), nrow = 2),
    p_resp = matrix(c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60), nrow = 2)),
  "Scenario 17" = list(
    p_true = matrix(c(0.05, 0.10, 0.15, 0.20, 0.40, 0.50), nrow = 2),
    p_resp = matrix(c(0.10, 0.20, 0.30, 0.40, 0.50, 0.55), nrow = 2))
)
oc <- sim_comb_boin12_multi(target = 0.33, target_resp = 0.2,
                            scenarios = scenarios, n_cohort = 14,
                            cohort_size = 3, n_trials = 100,
                            est_method = "observed", seed = 1)
print(oc)
#> Comb-BOIN12 operating characteristics across 2 scenarios
#>   DLT limit            : 0.33
#>   response limit       : 0.2
#>   dose levels          : 2 of drug A, 3 of drug B
#>   trials each          : 100
#>   max sample size      : 42
#>   ties / no candidate  : higher / nearest
#>   OBD selection        : observed response, rectangle bound
#> 
#>     Scenario True OBD OBD selected (%) Patients at OBD (%) No OBD (%)
#>  Scenario 16     A1B3               31                28.6          0
#>  Scenario 17     A2B2               18                23.4          1
#>  Above limit selected (%) Patients above limit (%) Patients DLTs Responses
#>                         6                     16.1       42  8.7      16.8
#>                        26                     30.9       42 10.2      15.2
```
