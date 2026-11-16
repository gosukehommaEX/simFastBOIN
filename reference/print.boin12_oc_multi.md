# Print the Operating Characteristics of the BOIN12 Design Across Scenarios

Display the summary produced by
[`sim_boin12_multi`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin12_multi.md)
as one table with a block of rows per scenario.

## Usage

``` r
# S3 method for class 'boin12_oc_multi'
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

  An object of class `boin12_oc_multi`.

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
scenarios <- list(
  "Plateau" = list(p_true = c(0.05, 0.12, 0.27, 0.35, 0.50),
                   p_resp = c(0.20, 0.35, 0.36, 0.37, 0.38)),
  "Increasing" = list(p_true = c(0.05, 0.12, 0.27, 0.35, 0.50),
                      p_resp = c(0.01, 0.05, 0.30, 0.60, 0.60))
)
oc <- sim_boin12_multi(
  target = 0.35, target_resp = 0.25, scenarios = scenarios,
  n_cohort = 12, cohort_size = 3, n_trials = 200, seed = 1
)
print(oc, percent = TRUE)
#> BOIN12 operating characteristics across 2 scenarios
#>   DLT limit            : 0.35
#>   response limit       : 0.25
#>   trials each          : 200
#>   max sample size      : 36
#>   OBD selection        : utility by posterior_mean
#> 
#>    Scenario                       Item  DL1  DL2  DL3  DL4  DL5 Total / No OBD
#>     Plateau          True DLT rate (%)  5.0 12.0 27.0 35.0 50.0               
#>                 True response rate (%) 20.0 35.0 36.0 37.0 38.0               
#>                       OBD selected (%) 16.5 49.5 22.5 10.5  1.0            0.0
#>                   Patients treated (%) 18.9 35.1 24.3 16.0  5.7           36.0
#>                  Patients with DLT (%)  4.9 21.6 32.6 26.1 14.8            7.5
#>             Patients with response (%) 10.7 38.5 26.7 17.8  6.4           11.8
#>  Increasing          True DLT rate (%)  5.0 12.0 27.0 35.0 50.0               
#>                 True response rate (%)  1.0  5.0 30.0 60.0 60.0               
#>                       OBD selected (%)  3.5  7.0 27.0 49.5  5.0            8.0
#>                   Patients treated (%) 16.4 19.5 25.3 30.2  8.6           35.2
#>                  Patients with DLT (%)  3.9  9.9 28.8 39.7 17.8            9.1
#>             Patients with response (%)  0.4  3.2 23.6 56.6 16.3           11.3
```
