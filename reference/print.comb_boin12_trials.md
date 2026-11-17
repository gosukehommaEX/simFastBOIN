# Print Simulated Comb-BOIN12 Trials

Print Simulated Comb-BOIN12 Trials

## Usage

``` r
# S3 method for class 'comb_boin12_trials'
print(x, ...)
```

## Arguments

- x:

  An object of class `comb_boin12_trials` from
  [`comb_boin12_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/comb_boin12_simulate.md).

- ...:

  Not used.

## Value

`x`, invisibly.

## See also

[`comb_boin12_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/comb_boin12_simulate.md)

## Examples

``` r
p_true <- matrix(c(0.05, 0.10, 0.15, 0.15, 0.20, 0.45), nrow = 2)
p_resp <- matrix(c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60), nrow = 2)
trials <- comb_boin12_simulate(target = 0.33, target_resp = 0.2,
                               p_true = p_true, p_resp = p_resp,
                               n_cohort = 14, cohort_size = 3,
                               n_trials = 100, seed = 1)
print(trials)
#> Simulated Comb-BOIN12 trials
#>   trials               : 100
#>   dose levels          : 2 of drug A, 3 of drug B
#>   DLT limit            : 0.33
#>   response limit       : 0.2
#>   cohorts              : 14 of size 3
#>   max sample size      : 42
#> 
#> Patients treated, average per trial:
#>      B1   B2    B3
#> A1 4.71 6.24 12.03
#> A2 3.90 8.37  6.75
#> 
#> DLTs observed, average per trial:
#>      B1   B2   B3
#> A1 0.35 0.98 2.33
#> A2 0.47 1.33 3.20
#> 
#> Responses observed, average per trial:
#>      B1   B2   B3
#> A1 0.50 1.83 6.18
#> A2 0.82 3.25 4.27
#> 
#> Stopping reason (%):
#> 
#> max_cohorts 
#>         100 
```
