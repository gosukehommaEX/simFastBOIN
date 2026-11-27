# Print Operating Characteristics of the U-BOIN Design Across Scenarios

Print the design and the cross-scenario summary table of a
`uboin_oc_multi` object.

## Usage

``` r
# S3 method for class 'uboin_oc_multi'
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

  An object of class `uboin_oc_multi`, as returned by
  [`sim_uboin_multi`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_uboin_multi.md).

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
scenarios <- list(
  "Scenario 1" = list(p_true = c(0.02, 0.15, 0.30, 0.45, 0.60),
                      p_resp = c(0.20, 0.65, 0.65, 0.65, 0.65)),
  "Scenario 3" = list(p_true = c(0.05, 0.15, 0.30, 0.45, 0.60),
                      p_resp = c(0.08, 0.46, 0.25, 0.20, 0.10))
)
oc <- sim_uboin_multi(
  tox_upper = 0.30,
  resp_lower = 0.20,
  scenarios = scenarios,
  n_cohort = 18,
  cohort_size = 3,
  n_trials = 200,
  seed = 123
)

print(oc)
#> U-BOIN operating characteristics across 2 scenarios
#>   DLT upper limit      : 0.3
#>   response lower limit : 0.2
#>   trials each          : 200
#>   max sample size      : 54
#>   utilities            : 0, 30, 50, 100
#> 
#>    Scenario                   Item  DL1  DL2  DL3  DL4  DL5 Total / No OBD
#>  Scenario 1      True DLT rate (%)  2.0 15.0 30.0 45.0 60.0               
#>             True response rate (%) 20.0 65.0 65.0 65.0 65.0               
#>                   OBD selected (%)  0.5 81.5 15.0  3.0  0.0            0.0
#>                   Patients treated  6.0 31.0 12.6  3.8  0.6           54.0
#>                  Patients with DLT  0.1  4.7  3.9  1.7  0.3           10.8
#>             Patients with response  1.2 20.4  8.0  2.5  0.3           32.5
#>  Scenario 3      True DLT rate (%)  5.0 15.0 30.0 45.0 60.0               
#>             True response rate (%)  8.0 46.0 25.0 20.0 10.0               
#>                   OBD selected (%)  0.5 90.5  6.0  1.0  0.0            2.0
#>                   Patients treated  6.8 34.2  9.3  2.8  0.3           53.5
#>                  Patients with DLT  0.3  5.3  2.7  1.3  0.2            9.7
#>             Patients with response  0.5 15.6  2.4  0.6  0.0           19.0
```
