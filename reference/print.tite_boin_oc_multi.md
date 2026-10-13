# Print Operating Characteristics of a TITE-BOIN Design Across Scenarios

Display the table produced by
[`sim_tite_boin_multi`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_tite_boin_multi.md),
one block of rows per scenario.

## Usage

``` r
# S3 method for class 'tite_boin_oc_multi'
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

  An object of class `tite_boin_oc_multi`.

- digits:

  Integer scalar. Number of decimal places. Defaults to 1.

- percent:

  Logical scalar. Express the average number of patients and DLTs at
  each dose as a percentage of the trial total instead of a count.
  Defaults to `FALSE`.

- kable:

  Logical scalar. Return the table through
  [`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html). Requires
  the knitr package. Defaults to `FALSE`.

- kable_format:

  Character scalar passed to
  [`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html). Defaults
  to `"pipe"`.

- ...:

  Further arguments, currently ignored.

## Value

The object `x`, invisibly.

## Examples

``` r
scenarios <- list(
  list(name = "MTD at dose 3", p_true = c(0.05, 0.15, 0.30, 0.45, 0.60)),
  list(name = "All toxic",     p_true = c(0.35, 0.45, 0.55, 0.65, 0.75))
)

oc <- sim_tite_boin_multi(
  target = 0.30,
  scenarios = scenarios,
  n_cohort = 10,
  cohort_size = 3,
  window = 3,
  accrual_rate = 2,
  n_trials = 200,
  seed = 123
)

print(oc)
#> TITE-BOIN operating characteristics across 2 scenarios
#>   method          : single mean imputation (Yuan et al., 2018)
#>   target DLT rate : 0.3
#>   trials each     : 200
#>   max sample size : 30
#>   window          : 3
#>   accrual         : 2 patients per unit of time (exponential)
#> 
#>       Scenario                         Item  DL1  DL2  DL3  DL4  DL5
#>  MTD at dose 3            True DLT rate (%)  5.0 15.0 30.0 45.0 60.0
#>                            MTD selected (%)  2.0 30.5 56.0 11.5  0.0
#>                            Patients treated  5.8 10.6  9.1  3.2  0.4
#>                           Patients with DLT  0.3  1.6  2.8  1.5  0.3
#>                       Trial duration (mean)                         
#>                Accrual suspended (% trials)                         
#>      All toxic            True DLT rate (%) 35.0 45.0 55.0 65.0 75.0
#>                            MTD selected (%) 75.0  8.0  1.0  0.0  0.0
#>                            Patients treated 16.4  4.1  0.4  0.0  0.0
#>                           Patients with DLT  5.7  1.9  0.2  0.0  0.0
#>                       Trial duration (mean)                         
#>                Accrual suspended (% trials)                         
#>  Total / No MTD
#>                
#>             0.0
#>            29.2
#>             6.5
#>            23.5
#>           100.0
#>                
#>            16.0
#>            21.0
#>             7.8
#>            15.6
#>            98.0
print(oc, percent = TRUE)
#> TITE-BOIN operating characteristics across 2 scenarios
#>   method          : single mean imputation (Yuan et al., 2018)
#>   target DLT rate : 0.3
#>   trials each     : 200
#>   max sample size : 30
#>   window          : 3
#>   accrual         : 2 patients per unit of time (exponential)
#> 
#>       Scenario                         Item  DL1  DL2  DL3  DL4  DL5
#>  MTD at dose 3            True DLT rate (%)  5.0 15.0 30.0 45.0 60.0
#>                            MTD selected (%)  2.0 30.5 56.0 11.5  0.0
#>                        Patients treated (%) 19.9 36.3 31.2 11.1  1.5
#>                       Patients with DLT (%)  4.9 24.9 43.0 22.3  4.8
#>                       Trial duration (mean)                         
#>                Accrual suspended (% trials)                         
#>      All toxic            True DLT rate (%) 35.0 45.0 55.0 65.0 75.0
#>                            MTD selected (%) 75.0  8.0  1.0  0.0  0.0
#>                        Patients treated (%) 78.2 19.6  2.0  0.2  0.0
#>                       Patients with DLT (%) 72.8 24.1  2.7  0.4  0.0
#>                       Trial duration (mean)                         
#>                Accrual suspended (% trials)                         
#>  Total / No MTD
#>                
#>             0.0
#>            29.2
#>             6.5
#>            23.5
#>           100.0
#>                
#>            16.0
#>            21.0
#>             7.8
#>            15.6
#>            98.0
```
