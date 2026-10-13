# Print Operating Characteristics of a Backfill Design Across Scenarios

Display the result of
[`sim_bf_boin_multi`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_bf_boin_multi.md)
or
[`sim_be_boin_multi`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_be_boin_multi.md)
as one table with seven rows per scenario.

## Usage

``` r
# S3 method for class 'backfill_oc_multi'
print(
  x,
  digits = 1,
  percent = FALSE,
  kable = FALSE,
  kable_format = "pipe",
  ...
)
```

## Arguments

- x:

  An object of class `backfill_oc_multi`.

- digits:

  Integer scalar. Number of decimal places. Defaults to 1.

- percent:

  Logical scalar. Show patients, backfilled patients and DLTs as
  percentages of their totals. Defaults to `FALSE`.

- kable:

  Logical scalar. Print only the table, formatted by
  [`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html). Defaults
  to `FALSE`.

- kable_format:

  Character scalar passed to
  [`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html) as
  `format`. Defaults to `"pipe"`.

- ...:

  Further arguments, currently ignored.

## Value

The object `x`, invisibly.

## See also

[`sim_bf_boin_multi`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_bf_boin_multi.md),
[`sim_be_boin_multi`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_be_boin_multi.md)

## Examples

``` r
scenarios <- list(
  list(name = "Rising", p_true = c(0.04, 0.12, 0.25, 0.43, 0.63),
       p_resp = c(0.10, 0.20, 0.30, 0.45, 0.58)),
  list(name = "Plateau", p_true = c(0.04, 0.12, 0.25, 0.43, 0.63),
       p_resp = c(0.30, 0.32, 0.35, 0.36, 0.36))
)
oc <- sim_bf_boin_multi(
  target = 0.25, scenarios = scenarios, n_cohort = 10, cohort_size = 3,
  window = 1, accrual_rate = 3, n_trials = 100, seed = 1
)
print(oc, percent = TRUE)
#> BF-BOIN operating characteristics across 2 scenarios
#>   target DLT rate : 0.25
#>   trials each     : 100
#>   escalation size : 30
#>   backfill cap    : 12 patients per dose
#>   window          : 1 (DLT), 1 (response)
#>   accrual         : 3 patients per unit of time (exponential)
#> 
#>  Scenario                    Item  DL1  DL2  DL3  DL4  DL5 Total / No MTD
#>    Rising       True DLT rate (%)  4.0 12.0 25.0 43.0 63.0               
#>            True response rate (%) 10.0 20.0 30.0 45.0 58.0               
#>                  MTD selected (%)  2.0 25.0 64.0  9.0  0.0            0.0
#>              Patients treated (%) 23.0 36.2 31.4  8.6  0.8           38.2
#>           Backfilled patients (%)  7.5  9.7  5.0  0.6  0.0           22.8
#>             Patients with DLT (%)  6.4 25.5 42.9 22.6  2.6            6.9
#>             Trial duration (mean)                                    17.1
#>   Plateau       True DLT rate (%)  4.0 12.0 25.0 43.0 63.0               
#>            True response rate (%) 30.0 32.0 35.0 36.0 36.0               
#>                  MTD selected (%)  2.0 30.0 56.0 12.0  0.0            0.0
#>              Patients treated (%) 25.6 35.5 29.2  8.6  1.1           40.2
#>           Backfilled patients (%) 11.0 11.2  4.6  0.7  0.0           27.5
#>             Patients with DLT (%)  6.6 26.5 40.7 22.6  3.6            7.0
#>             Trial duration (mean)                                    17.2
```
