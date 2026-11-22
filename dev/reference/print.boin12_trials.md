# Print Simulated BOIN12 Trials

Display the settings, the average numbers of patients, DLTs and
responses at each dose and the stopping reasons of trials simulated by
[`boin12_simulate`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/boin12_simulate.md).

## Usage

``` r
# S3 method for class 'boin12_trials'
print(x, ...)
```

## Arguments

- x:

  An object of class `boin12_trials`.

- ...:

  Not used.

## Value

`x`, invisibly.

## Examples

``` r
trials <- boin12_simulate(
  target = 0.35, target_resp = 0.25,
  p_true = c(0.05, 0.12, 0.27, 0.35, 0.50),
  p_resp = c(0.20, 0.35, 0.36, 0.37, 0.38),
  n_cohort = 12, cohort_size = 3, n_trials = 100, seed = 1
)
print(trials)
#> Simulated BOIN12 trials
#>   trials               : 100
#>   doses                : 5
#>   DLT limit            : 0.35
#>   response limit       : 0.25
#>   cohorts              : 12 of size 3
#>   max sample size      : 36
#> 
#> Average per trial:
#>                     DL1   DL2  DL3  DL4  DL5
#> Patients treated   7.32 12.24 8.64 5.76 2.04
#> DLTs observed      0.40  1.63 2.32 1.91 1.05
#> Responses observed 1.41  4.35 3.04 2.10 0.74
#> 
#> Stopping reason (%):
#> 
#> max_cohorts 
#>         100 
```
