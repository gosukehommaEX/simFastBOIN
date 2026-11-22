# Print Simulated BOIN-ET Trials

Show a compact summary of the object returned by
[`boin_et_simulate`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/boin_et_simulate.md)
instead of printing the full trial matrices.

## Usage

``` r
# S3 method for class 'boin_et_trials'
print(x, ...)
```

## Arguments

- x:

  An object of class `boin_et_trials`.

- ...:

  Further arguments, currently ignored.

## Value

The object `x`, invisibly.

## See also

[`boin_et_simulate`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/boin_et_simulate.md)

## Examples

``` r
trials <- boin_et_simulate(
  target = 0.30, target_resp = 0.60,
  p_true = c(0.05, 0.15, 0.25, 0.30, 0.35, 0.40),
  p_resp = c(0.05, 0.35, 0.55, 0.65, 0.75, 0.70),
  n_cohort = 12, cohort_size = 3, n_trials = 100, seed = 1
)
print(trials)
#> Simulated BOIN-ET trials
#>   trials               : 100
#>   doses                : 6
#>   target DLT rate      : 0.3
#>   target response rate : 0.6
#>   cohorts              : 12 of size 3
#>   max sample size      : 36
#> 
#> Average per trial:
#>                     DL1  DL2   DL3  DL4  DL5  DL6
#> Patients treated   3.57 8.07 14.82 6.96 1.86 0.18
#> DLTs observed      0.15 1.22  3.69 2.17 0.59 0.08
#> Responses observed 0.13 2.81  8.27 4.45 1.42 0.13
#> 
#> Stopping reason (%):
#> 
#>        max_cohorts no_admissible_dose 
#>                 97                  3 
```
