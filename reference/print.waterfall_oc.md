# Print Operating Characteristics of the Waterfall Design

Print the design, the true DLT rates, the selection percentages and the
average numbers of patients and DLTs of every combination, and the
summaries of the MTD contour, for a `waterfall_oc` object.

## Usage

``` r
# S3 method for class 'waterfall_oc'
print(x, digits = 1, kable = FALSE, kable_format = "pipe", ...)
```

## Arguments

- x:

  An object of class `waterfall_oc`, as returned by
  [`sim_waterfall`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_waterfall.md).

- digits:

  Integer scalar. Number of decimal places. Defaults to 1.

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
p_true <- matrix(c(0.03, 0.10, 0.28,
                   0.10, 0.30, 0.50), nrow = 2, byrow = TRUE)

oc <- sim_waterfall(
  target = 0.30,
  p_true = p_true,
  cohort_size = 3,
  n_trials = 200,
  seed = 123
)

print(oc)
#> Waterfall drug-combination operating characteristics
#>   target DLT rate : 0.3
#>   dose levels     : 2 of drug A, 3 of drug B
#>   trials          : 200
#>   cohorts         : at most 6/3 per subtrial, of size 3
#>   max sample size : 27
#>   lambda_e / _d   : 0.2365 / 0.3585
#> 
#>               Item Drug A   B1   B2   B3
#>  True DLT rate (%)     A1  3.0 10.0 28.0
#>                        A2 10.0 30.0 50.0
#>   MTD selected (%)     A1  0.0 17.0 83.0
#>                        A2 24.5 59.0 16.5
#>   Patients treated     A1  3.3  2.0  7.0
#>                        A2  5.6  6.6  2.4
#>  Patients with DLT     A1  0.1  0.2  2.0
#>                        A2  0.5  2.1  1.2
#> 
#> No MTD selected at any level of drug A : 0%
#> Average number of patients             : 26.9
#> Average number of DLTs                 : 6
#> Average number of subtrials            : 2
#> 
#> True MTD contour                       : A1B3, A2B2
#>   Trials selecting exactly the contour : 50%
#>   Patients treated on the contour      : 50.5%
#>   Patients treated above the contour   : 8.8%
#>   Patients treated below the contour   : 40.7%
print(oc, digits = 2)
#> Waterfall drug-combination operating characteristics
#>   target DLT rate : 0.3
#>   dose levels     : 2 of drug A, 3 of drug B
#>   trials          : 200
#>   cohorts         : at most 6/3 per subtrial, of size 3
#>   max sample size : 27
#>   lambda_e / _d   : 0.2365 / 0.3585
#> 
#>               Item Drug A    B1    B2    B3
#>  True DLT rate (%)     A1  3.00 10.00 28.00
#>                        A2 10.00 30.00 50.00
#>   MTD selected (%)     A1  0.00 17.00 83.00
#>                        A2 24.50 59.00 16.50
#>   Patients treated     A1  3.35  2.01  6.99
#>                        A2  5.61  6.62  2.37
#>  Patients with DLT     A1  0.09  0.17  1.97
#>                        A2  0.50  2.08  1.20
#> 
#> No MTD selected at any level of drug A : 0%
#> Average number of patients             : 26.94
#> Average number of DLTs                 : 6.01
#> Average number of subtrials            : 2
#> 
#> True MTD contour                       : A1B3, A2B2
#>   Trials selecting exactly the contour : 50%
#>   Patients treated on the contour      : 50.5%
#>   Patients treated above the contour   : 8.8%
#>   Patients treated below the contour   : 40.7%
```
