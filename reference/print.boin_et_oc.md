# Print the Operating Characteristics of the BOIN-ET Design

Display the summary produced by
[`sim_boin_et`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin_et.md)
as a table with one column per dose.

## Usage

``` r
# S3 method for class 'boin_et_oc'
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

  An object of class `boin_et_oc`.

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

  Further arguments, currently ignored.

## Value

The object `x`, invisibly.

## See also

[`sim_boin_et`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin_et.md)

## Examples

``` r
oc <- sim_boin_et(
  target = 0.30, target_resp = 0.60,
  p_true = c(0.05, 0.15, 0.25, 0.30, 0.35, 0.40),
  p_resp = c(0.05, 0.35, 0.55, 0.65, 0.75, 0.70),
  n_cohort = 12, cohort_size = 3, n_trials = 500, seed = 1
)
print(oc)
#> BOIN-ET operating characteristics
#>   target DLT rate      : 0.3
#>   target response rate : 0.6
#>   trials               : 500
#>   cohorts              : 12 of size 3
#>   max sample size      : 36
#>   lambda1 / lambda2    : 0.14 / 0.35
#>   eta1                 : 0.48
#>   OBD selection        : max_resp with observed response rates
#> 
#>                        DL1  DL2  DL3  DL4  DL5  DL6 Total / No OBD
#> True DLT rate (%)      5.0 15.0 25.0 30.0 35.0 40.0               
#> True response rate (%) 5.0 35.0 55.0 65.0 75.0 70.0               
#> OBD selected (%)       1.0 13.6 47.6 27.6  8.2  0.6            1.4
#> Patients treated       3.5  8.6 14.9  6.7  1.8  0.2           35.8
#> Patients with DLT      0.2  1.3  3.8  2.0  0.6  0.1            8.0
#> Patients with response 0.2  3.0  8.2  4.4  1.4  0.1           17.3
#> 
#> Trials stopped early: 1.4%
print(oc, percent = TRUE)
#> BOIN-ET operating characteristics
#>   target DLT rate      : 0.3
#>   target response rate : 0.6
#>   trials               : 500
#>   cohorts              : 12 of size 3
#>   max sample size      : 36
#>   lambda1 / lambda2    : 0.14 / 0.35
#>   eta1                 : 0.48
#>   OBD selection        : max_resp with observed response rates
#> 
#>                            DL1  DL2  DL3  DL4  DL5  DL6 Total / No OBD
#> True DLT rate (%)          5.0 15.0 25.0 30.0 35.0 40.0               
#> True response rate (%)     5.0 35.0 55.0 65.0 75.0 70.0               
#> OBD selected (%)           1.0 13.6 47.6 27.6  8.2  0.6            1.4
#> Patients treated (%)       9.9 24.1 41.6 18.8  5.0  0.6           35.8
#> Patients with DLT (%)      2.0 16.0 47.8 25.6  7.5  1.0            8.0
#> Patients with response (%) 1.0 17.5 47.6 25.2  7.8  0.8           17.3
#> 
#> Trials stopped early: 1.4%
```
