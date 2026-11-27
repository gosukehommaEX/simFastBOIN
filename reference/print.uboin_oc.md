# Print Operating Characteristics of the U-BOIN Design

Print the design, a table of the true rates and utilities, the selection
percentages and the average numbers of patients, DLTs and responses of
every dose, and the reasons the trials ended, for a `uboin_oc` object.

## Usage

``` r
# S3 method for class 'uboin_oc'
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

  An object of class `uboin_oc`, as returned by
  [`sim_uboin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_uboin.md).

- digits:

  Integer scalar. Number of decimal places. Defaults to 1.

- percent:

  Logical scalar. Show the numbers of patients, DLTs and responses as
  percentages of their totals. Defaults to `FALSE`.

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
oc <- sim_uboin(
  tox_upper = 0.30,
  resp_lower = 0.20,
  p_true = c(0.02, 0.15, 0.30, 0.45, 0.60),
  p_resp = c(0.20, 0.65, 0.65, 0.65, 0.65),
  n_cohort = 18,
  cohort_size = 3,
  resp_assoc = 0.2,
  n_trials = 200,
  seed = 123
)

print(oc)
#> U-BOIN operating characteristics
#>   DLT upper limit      : 0.3
#>   response lower limit : 0.2
#>   stage I target       : 0.25
#>   trials               : 200
#>   cohorts              : 18 of size 3
#>   max sample size      : 54
#>   stage II from / stop : 12 / 54 patients at a dose
#>   lambda_e / lambda_d  : 0.1968 / 0.2984
#>   utilities            : 0, 30, 50, 100
#>   cutoffs (DLT / resp) : 0.95 / 0.9
#> 
#>                         DL1  DL2  DL3  DL4  DL5 Total / No OBD
#> True DLT rate (%)       2.0 15.0 30.0 45.0 60.0               
#> True response rate (%) 20.0 65.0 65.0 65.0 65.0               
#> True utility           43.3 69.0 62.5 56.0 49.6               
#> OBD selected (%)        1.0 82.0 15.0  2.0  0.0            0.0
#> Patients treated        6.0 31.2 12.7  3.6  0.5           54.0
#> Patients with DLT       0.1  4.7  3.9  1.6  0.3           10.7
#> Patients with response  1.2 20.5  8.2  2.3  0.3           32.5
#> 
#> Trials reaching stage II: 100%
#> Reason for ending (%):
#>        max_sample_size                 n_stop lowest_dose_eliminated 
#>                    100                      0                      0 
#>     no_admissible_dose 
#>                      0 
print(oc, percent = TRUE)
#> U-BOIN operating characteristics
#>   DLT upper limit      : 0.3
#>   response lower limit : 0.2
#>   stage I target       : 0.25
#>   trials               : 200
#>   cohorts              : 18 of size 3
#>   max sample size      : 54
#>   stage II from / stop : 12 / 54 patients at a dose
#>   lambda_e / lambda_d  : 0.1968 / 0.2984
#>   utilities            : 0, 30, 50, 100
#>   cutoffs (DLT / resp) : 0.95 / 0.9
#> 
#>                             DL1  DL2  DL3  DL4  DL5 Total / No OBD
#> True DLT rate (%)           2.0 15.0 30.0 45.0 60.0               
#> True response rate (%)     20.0 65.0 65.0 65.0 65.0               
#> True utility               43.3 69.0 62.5 56.0 49.6               
#> OBD selected (%)            1.0 82.0 15.0  2.0  0.0            0.0
#> Patients treated (%)       11.1 57.8 23.5  6.6  0.9           54.0
#> Patients with DLT (%)       1.3 44.3 36.7 14.9  2.8           10.7
#> Patients with response (%)  3.6 63.1 25.2  7.2  0.9           32.5
#> 
#> Trials reaching stage II: 100%
#> Reason for ending (%):
#>        max_sample_size                 n_stop lowest_dose_eliminated 
#>                    100                      0                      0 
#>     no_admissible_dose 
#>                      0 
```
