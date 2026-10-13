# Print Operating Characteristics of a BOIN Design with Backfilling

Display the result of
[`sim_bf_boin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_bf_boin.md)
or
[`sim_be_boin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_be_boin.md):
the design, a table by dose of the true DLT and response rates, the MTD
selection, the patients treated and backfilled and the DLTs, then the
exposure above the overdose cutoff and the timing.

## Usage

``` r
# S3 method for class 'backfill_oc'
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

  An object of class `backfill_oc`.

- digits:

  Integer scalar. Number of decimal places. Defaults to 1.

- percent:

  Logical scalar. Show patients, backfilled patients and DLTs as
  percentages of their totals rather than as averages per trial.
  Defaults to `FALSE`.

- kable:

  Logical scalar. Print only the table, formatted by
  [`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html). Defaults
  to `FALSE`.

- kable_format:

  Character scalar passed to
  [`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html) as
  `format`. Defaults to `"pipe"`.

- ...:

  Further arguments, currently ignored.

## Value

The object `x`, invisibly.

## See also

[`sim_bf_boin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_bf_boin.md),
[`sim_be_boin`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_be_boin.md)

## Examples

``` r
oc <- sim_bf_boin(
  target = 0.25,
  p_true = c(0.04, 0.12, 0.25, 0.43, 0.63),
  p_resp = c(0.10, 0.20, 0.30, 0.45, 0.58),
  n_cohort = 10, cohort_size = 3, window = 1, accrual_rate = 3,
  n_trials = 200, seed = 1
)
print(oc)
#> BF-BOIN operating characteristics (Zhao et al., 2024)
#>   target DLT rate : 0.25
#>   trials          : 200
#>   cohorts         : 10 of size 3
#>   escalation size : 30
#>   backfill cap    : 12 patients per dose
#>   lambda_e / _d   : 0.1968 / 0.2984
#>   window          : 1 (DLT), 1 (response)
#>   accrual         : 3 patients per unit of time (exponential)
#> 
#>                         DL1  DL2  DL3  DL4  DL5 Total / No MTD
#> True DLT rate (%)       4.0 12.0 25.0 43.0 63.0               
#> True response rate (%) 10.0 20.0 30.0 45.0 58.0               
#> MTD selected (%)        1.5 25.5 63.0  9.5  0.5            0.0
#> Patients treated        8.2 13.6 12.4  3.9  0.3           38.4
#>   of which backfilled   2.7  3.8  2.1  0.2  0.0            8.9
#> Patients with DLT       0.4  1.7  3.0  1.7  0.2            7.0
#> 
#> Doses above a true DLT rate of 0.25: DL4, DL5
#>   Patients treated there              : 11%
#>   Trials treating any patient there   : 56.5%
#>   Trials selecting an MTD there       : 10%
#> 
#> Backfilling and timing
#>   Trials with backfilled patients     : 89%
#>   Average trial duration              : 17
#>   Trials in which the escalation waited: 100%
#>   Average time the escalation waited  : 6.7
print(oc, percent = TRUE)
#> BF-BOIN operating characteristics (Zhao et al., 2024)
#>   target DLT rate : 0.25
#>   trials          : 200
#>   cohorts         : 10 of size 3
#>   escalation size : 30
#>   backfill cap    : 12 patients per dose
#>   lambda_e / _d   : 0.1968 / 0.2984
#>   window          : 1 (DLT), 1 (response)
#>   accrual         : 3 patients per unit of time (exponential)
#> 
#>                            DL1  DL2  DL3  DL4  DL5 Total / No MTD
#> True DLT rate (%)          4.0 12.0 25.0 43.0 63.0               
#> True response rate (%)    10.0 20.0 30.0 45.0 58.0               
#> MTD selected (%)           1.5 25.5 63.0  9.5  0.5            0.0
#> Patients treated (%)      21.2 35.5 32.2 10.1  0.9           38.4
#>   of which backfilled (%)  7.1 10.0  5.4  0.6  0.0           23.1
#> Patients with DLT (%)      5.4 23.6 43.3 24.8  2.9            7.0
#> 
#> Doses above a true DLT rate of 0.25: DL4, DL5
#>   Patients treated there              : 11%
#>   Trials treating any patient there   : 56.5%
#>   Trials selecting an MTD there       : 10%
#> 
#> Backfilling and timing
#>   Trials with backfilled patients     : 89%
#>   Average trial duration              : 17
#>   Trials in which the escalation waited: 100%
#>   Average time the escalation waited  : 6.7
```
