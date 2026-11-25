# Print Simulated U-BOIN Trials

Print the design and the averages over the trials of a `uboin_trials`
object.

## Usage

``` r
# S3 method for class 'uboin_trials'
print(x, ...)
```

## Arguments

- x:

  An object of class `uboin_trials`, as returned by
  [`uboin_simulate`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/uboin_simulate.md).

- ...:

  Not used.

## Value

`x`, invisibly.

## Examples

``` r
trials <- uboin_simulate(
  tox_upper = 0.30,
  resp_lower = 0.20,
  p_true = c(0.02, 0.15, 0.30, 0.45, 0.60),
  p_resp = c(0.20, 0.65, 0.65, 0.65, 0.65),
  n_cohort = 18,
  cohort_size = 3,
  n_trials = 200,
  seed = 123
)
print(trials)
#> Simulated U-BOIN trials
#>   trials               : 200
#>   doses                : 5
#>   DLT upper limit      : 0.3
#>   response lower limit : 0.2
#>   cohorts              : 18 of size 3
#>   max sample size      : 54
#> 
#> Averages per trial:
#>                     DL1   DL2   DL3  DL4  DL5
#> Patients treated   6.03 31.02 12.59 3.78 0.58
#> DLTs observed      0.15  4.71  3.91 1.69 0.34
#> Responses observed 1.16 20.41  8.03 2.52 0.35
#> 
#> Trials reaching stage II: 100%
#> Reason for ending (%):
#> 
#>        max_sample_size                 n_stop lowest_dose_eliminated 
#>                    100                      0                      0 
#>     no_admissible_dose 
#>                      0 
```
