# Print Simulated TITE-BOIN Trials

Show a compact summary of the object returned by
[`tite_boin_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/tite_boin_simulate.md)
instead of printing the full trial matrices.

## Usage

``` r
# S3 method for class 'tite_boin_trials'
print(x, ...)
```

## Arguments

- x:

  An object of class `tite_boin_trials`.

- ...:

  Further arguments, currently ignored.

## Value

The object `x`, invisibly.

## Examples

``` r
trials <- tite_boin_simulate(
  target = 0.30,
  p_true = c(0.05, 0.15, 0.30, 0.45, 0.60),
  n_cohort = 10,
  cohort_size = 3,
  window = 3,
  accrual_rate = 2,
  n_trials = 100,
  seed = 1
)
print(trials)
#> Simulated TITE-BOIN trials
#>   method          : single mean imputation (Yuan et al., 2018)
#>   trials          : 100
#>   doses           : 5
#>   target DLT rate : 0.3
#>   cohorts         : 10 of size 3
#>   max sample size : 30
#>   window          : 3
#>   accrual         : 2 patients per unit of time (exponential)
#> 
#> Average per trial:
#>                   DL1   DL2  DL3  DL4  DL5
#> Patients treated 5.58 10.29 9.90 3.00 0.42
#> DLTs observed    0.31  1.53 3.06 1.45 0.21
#> 
#> Trial duration        : mean 23.56, median 23.48
#> Trials with suspension: 100%
#> Time suspended        : mean 6.5 per trial
#> 
#> Stopping reason (%):
#> 
#> max_sample_size     n_earlystop 
#>              80              20 
```
