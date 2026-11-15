# Print the Operating Characteristics of BOIN-ET Across Scenarios

Display the table produced by
[`sim_boin_et_multi`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin_et_multi.md),
with six rows per scenario.

## Usage

``` r
# S3 method for class 'boin_et_oc_multi'
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

  An object of class `boin_et_oc_multi`.

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

[`sim_boin_et_multi`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin_et_multi.md)

## Examples

``` r
scenarios <- list(
  list(name = "Monotone",
       p_true = c(0.05, 0.10, 0.15, 0.20, 0.25, 0.40),
       p_resp = c(0.05, 0.10, 0.25, 0.40, 0.65, 0.70)),
  list(name = "Peaked",
       p_true = c(0.05, 0.15, 0.25, 0.35, 0.45, 0.55),
       p_resp = c(0.30, 0.70, 0.60, 0.50, 0.40, 0.30))
)
res <- sim_boin_et_multi(
  target = 0.30, target_resp = 0.60, scenarios = scenarios,
  n_cohort = 12, cohort_size = 3, n_trials = 300, seed = 1
)
print(res)
#> BOIN-ET operating characteristics across 2 scenarios
#>   target DLT rate      : 0.3
#>   target response rate : 0.6
#>   trials each          : 300
#>   max sample size      : 36
#>   OBD selection        : max_resp with observed response rates
#> 
#>  Scenario                   Item  DL1  DL2  DL3  DL4  DL5  DL6 Total / No OBD
#>  Monotone      True DLT rate (%)  5.0 10.0 15.0 20.0 25.0 40.0               
#>           True response rate (%)  5.0 10.0 25.0 40.0 65.0 70.0               
#>                 OBD selected (%)  1.7  2.7  6.7 20.3 57.3 11.3            0.0
#>                 Patients treated  3.2  3.7  5.0  8.3 12.7  3.1           36.0
#>                Patients with DLT  0.1  0.4  0.7  1.6  3.2  1.2            7.3
#>           Patients with response  0.2  0.4  1.2  3.3  8.1  2.2           15.5
#>    Peaked      True DLT rate (%)  5.0 15.0 25.0 35.0 45.0 55.0               
#>           True response rate (%) 30.0 70.0 60.0 50.0 40.0 30.0               
#>                 OBD selected (%)  3.0 77.3 16.3  3.0  0.3  0.0            0.0
#>                 Patients treated  5.7 24.6  4.0  1.2  0.4  0.2           36.0
#>                Patients with DLT  0.3  3.8  1.0  0.4  0.2  0.1            5.8
#>           Patients with response  1.7 17.3  2.4  0.6  0.1  0.0           22.1
print(res, percent = TRUE)
#> BOIN-ET operating characteristics across 2 scenarios
#>   target DLT rate      : 0.3
#>   target response rate : 0.6
#>   trials each          : 300
#>   max sample size      : 36
#>   OBD selection        : max_resp with observed response rates
#> 
#>  Scenario                       Item  DL1  DL2  DL3  DL4  DL5  DL6
#>  Monotone          True DLT rate (%)  5.0 10.0 15.0 20.0 25.0 40.0
#>               True response rate (%)  5.0 10.0 25.0 40.0 65.0 70.0
#>                     OBD selected (%)  1.7  2.7  6.7 20.3 57.3 11.3
#>                 Patients treated (%)  8.9 10.3 13.9 22.9 35.2  8.7
#>                Patients with DLT (%)  2.0  5.2  9.9 22.0 44.5 16.5
#>           Patients with response (%)  1.1  2.6  8.1 21.4 52.5 14.3
#>    Peaked          True DLT rate (%)  5.0 15.0 25.0 35.0 45.0 55.0
#>               True response rate (%) 30.0 70.0 60.0 50.0 40.0 30.0
#>                     OBD selected (%)  3.0 77.3 16.3  3.0  0.3  0.0
#>                 Patients treated (%) 15.7 68.2 11.0  3.4  1.1  0.5
#>                Patients with DLT (%)  4.4 65.3 17.9  6.8  3.5  2.1
#>           Patients with response (%)  7.5 78.3 10.7  2.6  0.6  0.2
#>  Total / No OBD
#>                
#>                
#>             0.0
#>            36.0
#>             7.3
#>            15.5
#>                
#>                
#>             0.0
#>            36.0
#>             5.8
#>            22.1
```
