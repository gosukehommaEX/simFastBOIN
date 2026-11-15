# Print the Operating Characteristics of BF-BOIN-ET Across Scenarios

Display the table produced by
[`sim_bf_boin_et_multi`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_bf_boin_et_multi.md),
with seven rows per scenario.

## Usage

``` r
# S3 method for class 'bf_boin_et_oc_multi'
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

  An object of class `bf_boin_et_oc_multi`.

- digits:

  Integer scalar. Number of decimal places. Defaults to 1.

- percent:

  Logical scalar. Express the average numbers of patients, backfilled
  patients, DLTs and responses at each dose as percentages of the trial
  totals. Defaults to `FALSE`.

- kable:

  Logical scalar. Print the table with
  [`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html), for use
  in R Markdown. Defaults to `FALSE`.

- kable_format:

  Character scalar. Format passed to
  [`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html). Defaults
  to `"pipe"`.

- ...:

  Further arguments, currently ignored.

## Value

The object `x`, invisibly.

## Examples

``` r
scenarios <- list(
  "Example" = list(p_true = c(0.03, 0.05, 0.20, 0.22, 0.45),
                   p_resp = c(0.05, 0.10, 0.50, 0.68, 0.70)),
  "Plateau" = list(p_true = c(0.05, 0.10, 0.15, 0.30, 0.45),
                   p_resp = c(0.10, 0.30, 0.45, 0.50, 0.52))
)
res <- sim_bf_boin_et_multi(
  target = 0.30, target_resp = 0.25, scenarios = scenarios,
  n_cohort = 10, cohort_size = 3, window = 1, accrual_rate = 3,
  n_trials = 100, seed = 1
)
print(res, digits = 2)
#> BF-BOIN-ET operating characteristics across 2 scenarios
#>   target DLT rate      : 0.3
#>   target response rate : 0.25
#>   trials each          : 100
#>   cohorts              : 10 of size 3
#> 
#>  Scenario                   Item   DL1   DL2   DL3   DL4   DL5 Total / No OBD
#>   Example      True DLT rate (%)  3.00  5.00 20.00 22.00 45.00               
#>           True response rate (%)  5.00 10.00 50.00 68.00 70.00               
#>                 OBD selected (%)  3.00  3.00 31.00 58.00  5.00           0.00
#>                 Patients treated  4.22  5.88 11.25  9.83  5.11          36.29
#>              Patients backfilled  0.80  1.02  4.05  2.90  0.28           9.05
#>                Patients with DLT  0.08  0.35  2.38  2.12  2.19           7.12
#>           Patients with response  0.28  0.48  5.39  6.90  3.65          16.70
#>   Plateau      True DLT rate (%)  5.00 10.00 15.00 30.00 45.00               
#>           True response rate (%) 10.00 30.00 45.00 50.00 52.00               
#>                 OBD selected (%)  2.00 22.00 51.00 25.00  0.00           0.00
#>                 Patients treated  4.78  8.68 12.07  9.94  3.78          39.25
#>              Patients backfilled  1.36  4.48  5.02  1.42  0.18          12.46
#>                Patients with DLT  0.14  0.71  1.66  3.30  1.76           7.57
#>           Patients with response  0.50  2.65  5.31  4.98  1.84          15.28
```
