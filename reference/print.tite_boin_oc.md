# Print Operating Characteristics of a TITE-BOIN Design

Display the summary produced by
[`sim_tite_boin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_tite_boin.md):
the design, the table of
[`print.boin_oc`](https://gosukehommaex.github.io/simFastBOIN/reference/print.boin_oc.md),
the exposure to doses above the overdose cutoff, and the duration of the
trials and the suspensions of accrual.

## Usage

``` r
# S3 method for class 'tite_boin_oc'
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

  An object of class `tite_boin_oc`.

- digits:

  Integer scalar. Number of decimal places. Defaults to 1.

- percent:

  Logical scalar. Express the average number of patients and DLTs at
  each dose as a percentage of the trial total instead of a count.
  Defaults to `FALSE`.

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
oc <- sim_tite_boin(
  target = 0.30,
  p_true = c(0.05, 0.15, 0.30, 0.45, 0.60),
  n_cohort = 10,
  cohort_size = 3,
  window = 3,
  accrual_rate = 2,
  n_trials = 200,
  seed = 123
)

print(oc)
#> TITE-BOIN operating characteristics
#>   method          : single mean imputation (Yuan et al., 2018)
#>   target DLT rate : 0.3
#>   trials          : 200
#>   cohorts         : 10 of size 3
#>   max sample size : 30
#>   lambda_e / _d   : 0.2365 / 0.3585
#>   window          : 3
#>   accrual         : 2 patients per unit of time (exponential)
#> 
#>                   DL1  DL2  DL3  DL4  DL5 Total / No MTD
#> True DLT rate (%) 5.0 15.0 30.0 45.0 60.0               
#> MTD selected (%)  2.0 30.5 56.0 11.5  0.0            0.0
#> Patients treated  5.8 10.6  9.1  3.2  0.4           29.2
#> Patients with DLT 0.3  1.6  2.8  1.5  0.3            6.5
#> 
#> Doses above a true DLT rate of 0.3: DL4, DL5
#>   Patients treated there              : 12.6%
#>   Trials treating any patient there   : 38.5%
#>   Trials treating over 60% there      : 2.5%
#>   Trials treating over 80% there      : 0%
#>   Trials selecting an MTD there       : 11.5%
#> 
#> Trial duration and suspension of accrual
#>   Average trial duration              : 23.5
#>   Trials with accrual suspended       : 100%
#>   Average number of suspensions       : 3.9
#>   Average time suspended              : 6.2
print(oc, percent = TRUE)
#> TITE-BOIN operating characteristics
#>   method          : single mean imputation (Yuan et al., 2018)
#>   target DLT rate : 0.3
#>   trials          : 200
#>   cohorts         : 10 of size 3
#>   max sample size : 30
#>   lambda_e / _d   : 0.2365 / 0.3585
#>   window          : 3
#>   accrual         : 2 patients per unit of time (exponential)
#> 
#>                        DL1  DL2  DL3  DL4  DL5 Total / No MTD
#> True DLT rate (%)      5.0 15.0 30.0 45.0 60.0               
#> MTD selected (%)       2.0 30.5 56.0 11.5  0.0            0.0
#> Patients treated (%)  19.9 36.3 31.2 11.1  1.5           29.2
#> Patients with DLT (%)  4.9 24.9 43.0 22.3  4.8            6.5
#> 
#> Doses above a true DLT rate of 0.3: DL4, DL5
#>   Patients treated there              : 12.6%
#>   Trials treating any patient there   : 38.5%
#>   Trials treating over 60% there      : 2.5%
#>   Trials treating over 80% there      : 0%
#>   Trials selecting an MTD there       : 11.5%
#> 
#> Trial duration and suspension of accrual
#>   Average trial duration              : 23.5
#>   Trials with accrual suspended       : 100%
#>   Average number of suspensions       : 3.9
#>   Average time suspended              : 6.2
```
