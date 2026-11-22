# Print the Operating Characteristics of BOIN-ET Across Scenarios

Display the table produced by
[`sim_boin_et_multi`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_boin_et_multi.md),
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

[`sim_boin_et_multi`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_boin_et_multi.md)

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
#>                 OBD selected (%)  3.3  2.3  7.3 19.0 56.3 11.3            0.3
#>                 Patients treated  3.2  3.9  5.9  8.6 11.8  2.5           35.9
#>                Patients with DLT  0.1  0.4  0.9  1.8  3.0  0.9            7.2
#>           Patients with response  0.1  0.4  1.5  3.4  7.8  1.8           15.0
#>    Peaked      True DLT rate (%)  5.0 15.0 25.0 35.0 45.0 55.0               
#>           True response rate (%) 30.0 70.0 60.0 50.0 40.0 30.0               
#>                 OBD selected (%)  3.3 77.3 16.3  2.7  0.3  0.0            0.0
#>                 Patients treated  5.8 24.6  4.0  1.2  0.3  0.1           36.0
#>                Patients with DLT  0.3  3.8  1.0  0.4  0.2  0.1            5.7
#>           Patients with response  1.7 17.3  2.4  0.6  0.1  0.0           22.0
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
#>                     OBD selected (%)  3.3  2.3  7.3 19.0 56.3 11.3
#>                 Patients treated (%)  9.0 11.0 16.3 23.8 32.9  7.0
#>                Patients with DLT (%)  2.1  6.2 12.0 25.6 41.8 12.4
#>           Patients with response (%)  1.0  2.6  9.9 22.9 51.7 11.9
#>    Peaked          True DLT rate (%)  5.0 15.0 25.0 35.0 45.0 55.0
#>               True response rate (%) 30.0 70.0 60.0 50.0 40.0 30.0
#>                     OBD selected (%)  3.3 77.3 16.3  2.7  0.3  0.0
#>                 Patients treated (%) 16.0 68.3 11.0  3.4  1.0  0.3
#>                Patients with DLT (%)  4.8 66.2 17.5  7.0  3.2  1.4
#>           Patients with response (%)  7.6 78.5 10.7  2.6  0.6  0.1
#>  Total / No OBD
#>                
#>                
#>             0.3
#>            35.9
#>             7.2
#>            15.0
#>                
#>                
#>             0.0
#>            36.0
#>             5.7
#>            22.0
```
