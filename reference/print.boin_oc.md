# Print Operating Characteristics of a BOIN Design

Display the summary table produced by
[`sim_boin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin.md).

## Usage

``` r
# S3 method for class 'boin_oc'
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

  An object of class `boin_oc`.

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
oc <- sim_boin(
  target = 0.30,
  p_true = c(0.05, 0.15, 0.30, 0.45, 0.60),
  n_cohort = 10,
  cohort_size = 3,
  n_trials = 200,
  seed = 123
)

print(oc)
#> BOIN operating characteristics
#>   target DLT rate : 0.3
#>   trials          : 200
#>   cohorts         : 10 of size 3
#>   max sample size : 30
#>   lambda_e / _d   : 0.2365 / 0.3585
#> 
#>                   DL1  DL2  DL3  DL4  DL5 Total / No MTD
#> True DLT rate (%) 5.0 15.0 30.0 45.0 60.0               
#> MTD selected (%)  1.0 26.0 52.0 21.0  0.0            0.0
#> Patients treated  4.2  9.1 10.6  4.5  0.6           29.1
#> Patients with DLT 0.2  1.4  3.2  2.1  0.4            7.3
#> 
#> Doses above a true DLT rate of 0.3: DL4, DL5
#>   Patients treated there              : 17.8%
#>   Trials treating any patient there   : 56%
#>   Trials treating over 60% there      : 1.5%
#>   Trials treating over 80% there      : 0%
#>   Trials selecting an MTD there       : 21%
print(oc, percent = TRUE)
#> BOIN operating characteristics
#>   target DLT rate : 0.3
#>   trials          : 200
#>   cohorts         : 10 of size 3
#>   max sample size : 30
#>   lambda_e / _d   : 0.2365 / 0.3585
#> 
#>                        DL1  DL2  DL3  DL4  DL5 Total / No MTD
#> True DLT rate (%)      5.0 15.0 30.0 45.0 60.0               
#> MTD selected (%)       1.0 26.0 52.0 21.0  0.0            0.0
#> Patients treated (%)  14.4 31.4 36.5 15.6  2.2           29.1
#> Patients with DLT (%)  3.2 19.2 43.9 28.3  5.4            7.3
#> 
#> Doses above a true DLT rate of 0.3: DL4, DL5
#>   Patients treated there              : 17.8%
#>   Trials treating any patient there   : 56%
#>   Trials treating over 60% there      : 1.5%
#>   Trials treating over 80% there      : 0%
#>   Trials selecting an MTD there       : 21%
```
