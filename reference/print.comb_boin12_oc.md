# Print the Operating Characteristics of the Comb-BOIN12 Design

Print the Operating Characteristics of the Comb-BOIN12 Design

## Usage

``` r
# S3 method for class 'comb_boin12_oc'
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

  An object of class `comb_boin12_oc` from
  [`sim_comb_boin12`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_comb_boin12.md).

- digits:

  Integer scalar. Decimal places. Defaults to 1.

- percent:

  Logical scalar. Show the patients, DLTs and responses at each
  combination as percentages of their totals. Defaults to `FALSE`.

- kable:

  Logical scalar. Print the table with
  [`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html). Defaults
  to `FALSE`.

- kable_format:

  Character scalar. Format passed to
  [`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html). Defaults
  to `"pipe"`.

- ...:

  Not used.

## Value

`x`, invisibly.

## See also

[`sim_comb_boin12`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_comb_boin12.md)

## Examples

``` r
p_true <- matrix(c(0.05, 0.10, 0.15, 0.15, 0.20, 0.45), nrow = 2)
p_resp <- matrix(c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60), nrow = 2)
oc <- sim_comb_boin12(target = 0.33, target_resp = 0.2, p_true = p_true,
                      p_resp = p_resp, n_cohort = 14, cohort_size = 3,
                      n_trials = 100, est_method = "observed", seed = 1)
print(oc)
#> Comb-BOIN12 operating characteristics
#>   DLT limit            : 0.33
#>   response limit       : 0.2
#>   dose levels          : 2 of drug A, 3 of drug B
#>   trials               : 100
#>   cohorts              : 14 of size 3
#>   max sample size      : 42
#>   lambda_e / lambda_d  : 0.2604 / 0.3947
#>   utility scores       : 100, 40, 60, 0; benchmark 69.4
#>   ties / no candidate  : higher / nearest
#>   OBD selection        : observed response, rectangle bound
#> 
#>                    Item Drug A   B1   B2   B3
#>       True DLT rate (%)     A1  5.0 15.0 20.0
#>                             A2 10.0 15.0 45.0
#>  True response rate (%)     A1 10.0 30.0 50.0
#>                             A2 20.0 40.0 60.0
#>        OBD selected (%)     A1 18.0 21.0 29.0
#>                             A2  8.0 18.0  6.0
#>        Patients treated     A1  4.7  6.2 12.0
#>                             A2  3.9  8.4  6.8
#>       Patients with DLT     A1  0.3  1.0  2.3
#>                             A2  0.5  1.3  3.2
#>  Patients with response     A1  0.5  1.8  6.2
#>                             A2  0.8  3.2  4.3
#> 
#> No OBD selected                      : 0%
#> Trials stopped early                 : 0%
#> Average number of patients           : 42
#> 
#> True OBD (highest true utility)      : A1B3
#>   Selected as the OBD                : 29%
#>   Patients treated there             : 28.6%
#> Combinations above the DLT limit     : 1 of 6
#>   Selected as the OBD                : 6%
#>   Patients treated there             : 16.1%
print(oc, percent = TRUE)
#> Comb-BOIN12 operating characteristics
#>   DLT limit            : 0.33
#>   response limit       : 0.2
#>   dose levels          : 2 of drug A, 3 of drug B
#>   trials               : 100
#>   cohorts              : 14 of size 3
#>   max sample size      : 42
#>   lambda_e / lambda_d  : 0.2604 / 0.3947
#>   utility scores       : 100, 40, 60, 0; benchmark 69.4
#>   ties / no candidate  : higher / nearest
#>   OBD selection        : observed response, rectangle bound
#> 
#>                        Item Drug A   B1   B2   B3
#>           True DLT rate (%)     A1  5.0 15.0 20.0
#>                                 A2 10.0 15.0 45.0
#>      True response rate (%)     A1 10.0 30.0 50.0
#>                                 A2 20.0 40.0 60.0
#>            OBD selected (%)     A1 18.0 21.0 29.0
#>                                 A2  8.0 18.0  6.0
#>        Patients treated (%)     A1 11.2 14.9 28.6
#>                                 A2  9.3 19.9 16.1
#>       Patients with DLT (%)     A1  4.0 11.3 26.9
#>                                 A2  5.4 15.4 37.0
#>  Patients with response (%)     A1  3.0 10.9 36.7
#>                                 A2  4.9 19.3 25.3
#> 
#> No OBD selected                      : 0%
#> Trials stopped early                 : 0%
#> Average number of patients           : 42
#> 
#> True OBD (highest true utility)      : A1B3
#>   Selected as the OBD                : 29%
#>   Patients treated there             : 28.6%
#> Combinations above the DLT limit     : 1 of 6
#>   Selected as the OBD                : 6%
#>   Patients treated there             : 16.1%
```
