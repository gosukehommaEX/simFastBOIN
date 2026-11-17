# Print the Operating Characteristics of the BOIN12 Design

Display the summary produced by
[`sim_boin12`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin12.md)
as a table with one column per dose.

## Usage

``` r
# S3 method for class 'boin12_oc'
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

  An object of class `boin12_oc`.

- digits:

  Integer scalar. Number of decimal places. Defaults to 1.

- percent:

  Logical scalar. Express the average numbers of patients, DLTs and
  responses at each dose as percentages of the trial totals. Defaults to
  `FALSE`.

- kable:

  Logical scalar. Print the table with
  [`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html), for use
  in R Markdown. Defaults to `FALSE`.

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
oc <- sim_boin12(
  target = 0.35, target_resp = 0.25,
  p_true = c(0.05, 0.12, 0.27, 0.35, 0.50),
  p_resp = c(0.20, 0.35, 0.36, 0.37, 0.38),
  n_cohort = 12, cohort_size = 3, n_trials = 200, seed = 1
)
print(oc)
#> BOIN12 operating characteristics
#>   DLT limit            : 0.35
#>   response limit       : 0.25
#>   trials               : 200
#>   cohorts              : 12 of size 3
#>   max sample size      : 36
#>   lambda_e / lambda_d  : 0.2763 / 0.4189
#>   utility scores       : 100, 40, 60, 0; benchmark 70.5
#>   ties / no candidate  : higher / nearest
#>   OBD selection        : utility by posterior_mean
#> 
#>                         DL1  DL2  DL3  DL4  DL5 Total / No OBD
#> True DLT rate (%)       5.0 12.0 27.0 35.0 50.0               
#> True response rate (%) 20.0 35.0 36.0 37.0 38.0               
#> OBD selected (%)       16.5 49.5 22.5 10.5  1.0            0.0
#> Patients treated        6.8 12.6  8.7  5.7  2.1           36.0
#> Patients with DLT       0.4  1.6  2.4  1.9  1.1            7.5
#> Patients with response  1.3  4.5  3.1  2.1  0.8           11.8
#> 
#> Trials stopped early: 0%
print(oc, percent = TRUE)
#> BOIN12 operating characteristics
#>   DLT limit            : 0.35
#>   response limit       : 0.25
#>   trials               : 200
#>   cohorts              : 12 of size 3
#>   max sample size      : 36
#>   lambda_e / lambda_d  : 0.2763 / 0.4189
#>   utility scores       : 100, 40, 60, 0; benchmark 70.5
#>   ties / no candidate  : higher / nearest
#>   OBD selection        : utility by posterior_mean
#> 
#>                             DL1  DL2  DL3  DL4  DL5 Total / No OBD
#> True DLT rate (%)           5.0 12.0 27.0 35.0 50.0               
#> True response rate (%)     20.0 35.0 36.0 37.0 38.0               
#> OBD selected (%)           16.5 49.5 22.5 10.5  1.0            0.0
#> Patients treated (%)       18.9 35.1 24.3 16.0  5.7           36.0
#> Patients with DLT (%)       4.9 21.6 32.6 26.1 14.8            7.5
#> Patients with response (%) 10.7 38.5 26.7 17.8  6.4           11.8
#> 
#> Trials stopped early: 0%
```
