# Print the Operating Characteristics of the BF-BOIN-ET Design

Display the summary produced by
[`sim_bf_boin_et`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_bf_boin_et.md)
as a table with one column per dose.

## Usage

``` r
# S3 method for class 'bf_boin_et_oc'
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

  An object of class `bf_boin_et_oc`.

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
oc <- sim_bf_boin_et(
  target = 0.30, target_resp = 0.25,
  p_true = c(0.03, 0.05, 0.20, 0.22, 0.45),
  p_resp = c(0.05, 0.10, 0.50, 0.68, 0.70),
  n_cohort = 10, cohort_size = 3, window = 1, accrual_rate = 3,
  n_trials = 200, seed = 1
)
print(oc, digits = 2)
#> BF-BOIN-ET operating characteristics
#>   target DLT rate      : 0.3
#>   target response rate : 0.25
#>   trials               : 200
#>   cohorts              : 10 of size 3
#>   patients per dose    : 12 at most with backfill
#>   lambda1 / lambda2    : 0.2365 / 0.3585
#>   eta1                 : 0.1968
#> 
#>                         DL1   DL2   DL3   DL4   DL5 Total / No OBD
#> True DLT rate (%)      3.00  5.00 20.00 22.00 45.00               
#> True response rate (%) 5.00 10.00 50.00 68.00 70.00               
#> OBD selected (%)       3.00  4.50 28.50 58.50  5.50           0.00
#> Patients treated       4.10  5.80 11.21 10.27  5.47          36.85
#> Patients backfilled    0.64  1.27  4.57  2.87  0.21           9.57
#> Patients with DLT      0.12  0.30  2.27  2.20  2.42           7.32
#> Patients with response 0.22  0.58  5.39  7.06  3.78          17.03
#> 
#> Trials stopped early: 0%
#> Average trial duration: 17.09
print(oc, percent = TRUE)
#> BF-BOIN-ET operating characteristics
#>   target DLT rate      : 0.3
#>   target response rate : 0.25
#>   trials               : 200
#>   cohorts              : 10 of size 3
#>   patients per dose    : 12 at most with backfill
#>   lambda1 / lambda2    : 0.2365 / 0.3585
#>   eta1                 : 0.1968
#> 
#>                             DL1  DL2  DL3  DL4  DL5 Total / No OBD
#> True DLT rate (%)           3.0  5.0 20.0 22.0 45.0               
#> True response rate (%)      5.0 10.0 50.0 68.0 70.0               
#> OBD selected (%)            3.0  4.5 28.5 58.5  5.5            0.0
#> Patients treated (%)       11.1 15.7 30.4 27.9 14.9           36.9
#> Patients backfilled (%)     6.7 13.3 47.8 30.0  2.2            9.6
#> Patients with DLT (%)       1.6  4.2 31.0 30.1 33.1            7.3
#> Patients with response (%)  1.3  3.4 31.7 41.4 22.2           17.0
#> 
#> Trials stopped early: 0%
#> Average trial duration: 17.1
```
