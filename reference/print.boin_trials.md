# Print Simulated BOIN Trials

Show a compact summary of the object returned by
[`boin_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_simulate.md)
instead of printing the full trial matrices.

## Usage

``` r
# S3 method for class 'boin_trials'
print(x, ...)
```

## Arguments

- x:

  An object of class `boin_trials`.

- ...:

  Further arguments, currently ignored.

## Value

The object `x`, invisibly.

## Examples

``` r
trials <- boin_simulate(
  target = 0.30,
  p_true = c(0.05, 0.15, 0.30, 0.45, 0.60),
  n_cohort = 10,
  cohort_size = 3,
  n_trials = 100,
  seed = 1
)
print(trials)
#> Simulated BOIN trials
#>   trials          : 100
#>   doses           : 5
#>   target DLT rate : 0.3
#>   cohorts         : 10 of size 3
#>   max sample size : 30
#> 
#> Average per trial:
#>                   DL1  DL2   DL3  DL4  DL5
#> Patients treated 4.32 8.76 10.74 4.32 0.78
#> DLTs observed    0.22 1.41  3.35 2.00 0.51
#> 
#> Stopping reason (%):
#> 
#> max_sample_size     n_earlystop 
#>              75              25 
```
