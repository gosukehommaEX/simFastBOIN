# Print Operating Characteristics of the BOIN Drug-Combination Design

Display the summary produced by
[`sim_comb_boin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_comb_boin.md):
the true DLT rates, the selection percentages and the average numbers of
patients and DLTs at every dose combination, followed by the selection
of the true MTD and the exposure to combinations above the overdose
cutoff.

## Usage

``` r
# S3 method for class 'comb_boin_oc'
print(x, digits = 1, kable = FALSE, kable_format = "pipe", ...)
```

## Arguments

- x:

  An object of class `comb_boin_oc`.

- digits:

  Integer scalar. Number of decimal places. Defaults to 1.

- kable:

  Logical scalar. Return the table through
  [`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html) rather
  than printing it directly, which is convenient inside a report.
  Requires the knitr package. Defaults to `FALSE`.

- kable_format:

  Character scalar passed to
  [`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html), for
  example `"pipe"`, `"html"` or `"latex"`. Defaults to `"pipe"`.

- ...:

  Further arguments, currently ignored.

## Value

The object `x`, invisibly.

## Examples

``` r
p_true <- matrix(c(0.02, 0.04, 0.08, 0.14,
                   0.08, 0.25, 0.42, 0.48,
                   0.25, 0.45, 0.50, 0.60), nrow = 3, byrow = TRUE)

oc <- sim_comb_boin(
  target = 0.25,
  p_true = p_true,
  n_cohort = 16,
  cohort_size = 3,
  n_trials = 200,
  seed = 123
)

print(oc)
#> BOIN drug-combination operating characteristics
#>   target DLT rate : 0.25
#>   dose levels     : 3 of drug A, 4 of drug B
#>   trials          : 200
#>   cohorts         : 16 of size 3
#>   max sample size : 48
#>   lambda_e / _d   : 0.1968 / 0.2984
#> 
#>               Item Drug A   B1   B2   B3   B4
#>  True DLT rate (%)     A1  2.0  4.0  8.0 14.0
#>                        A2  8.0 25.0 42.0 48.0
#>                        A3 25.0 45.0 50.0 60.0
#>   MTD selected (%)     A1  0.0  1.0  5.0 16.0
#>                        A2  3.0 37.0 10.0  0.5
#>                        A3 20.5  7.0  0.0  0.0
#>   Patients treated     A1  4.0  3.7  2.7  2.8
#>                        A2  5.8  9.0  3.9  1.5
#>                        A3  4.7  4.0  1.0  0.3
#>  Patients with DLT     A1  0.1  0.2  0.2  0.4
#>                        A2  0.4  2.1  1.6  0.7
#>                        A3  1.2  1.7  0.5  0.2
#> 
#> No MTD selected                      : 0%
#> Average number of patients           : 43.6
#> Average number of DLTs               : 9.4
#> 
#> True MTD (closest to the target)     : A3B1, A2B2
#>   Selected as the MTD                : 57.5%
#>   Patients treated there             : 31.6%
#> 
#> Combinations above a true DLT rate of 0.25: 5 of 12
#>   Patients treated there             : 24.8%
#>   Trials treating any patient there  : 82%
#>   Trials treating over 60% there     : 2.5%
#>   Trials treating over 80% there     : 0.5%
#>   Trials selecting an MTD there      : 17.5%
print(oc, digits = 2)
#> BOIN drug-combination operating characteristics
#>   target DLT rate : 0.25
#>   dose levels     : 3 of drug A, 4 of drug B
#>   trials          : 200
#>   cohorts         : 16 of size 3
#>   max sample size : 48
#>   lambda_e / _d   : 0.1968 / 0.2984
#> 
#>               Item Drug A    B1    B2    B3    B4
#>  True DLT rate (%)     A1  2.00  4.00  8.00 14.00
#>                        A2  8.00 25.00 42.00 48.00
#>                        A3 25.00 45.00 50.00 60.00
#>   MTD selected (%)     A1  0.00  1.00  5.00 16.00
#>                        A2  3.00 37.00 10.00  0.50
#>                        A3 20.50  7.00  0.00  0.00
#>   Patients treated     A1  3.98  3.70  2.71  2.84
#>                        A2  5.80  9.03  3.93  1.50
#>                        A3  4.74  4.04  1.03  0.30
#>  Patients with DLT     A1  0.07  0.18  0.18  0.41
#>                        A2  0.40  2.12  1.65  0.72
#>                        A3  1.22  1.73  0.50  0.20
#> 
#> No MTD selected                      : 0%
#> Average number of patients           : 43.6
#> Average number of DLTs               : 9.37
#> 
#> True MTD (closest to the target)     : A3B1, A2B2
#>   Selected as the MTD                : 57.5%
#>   Patients treated there             : 31.58%
#> 
#> Combinations above a true DLT rate of 0.25: 5 of 12
#>   Patients treated there             : 24.77%
#>   Trials treating any patient there  : 82%
#>   Trials treating over 60% there     : 2.5%
#>   Trials treating over 80% there     : 0.5%
#>   Trials selecting an MTD there      : 17.5%
```
