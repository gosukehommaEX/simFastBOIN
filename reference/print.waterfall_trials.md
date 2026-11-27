# Print Simulated Waterfall Trials

Print the design and the averages over the trials of a
`waterfall_trials` object.

## Usage

``` r
# S3 method for class 'waterfall_trials'
print(x, ...)
```

## Arguments

- x:

  An object of class `waterfall_trials`, as returned by
  [`waterfall_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/waterfall_simulate.md).

- ...:

  Not used.

## Value

`x`, invisibly.

## Examples

``` r
p_true <- matrix(c(0.03, 0.10, 0.28,
                   0.10, 0.30, 0.50), nrow = 2, byrow = TRUE)

trials <- waterfall_simulate(
  target = 0.30,
  p_true = p_true,
  cohort_size = 3,
  n_trials = 200,
  seed = 123
)
print(trials)
#> Simulated waterfall drug-combination trials
#>   trials          : 200
#>   dose levels     : 2 of drug A, 3 of drug B
#>   target DLT rate : 0.3
#>   cohorts         : at most 6/3 per subtrial, of size 3
#>   max sample size : 27
#> 
#> Average patients treated per trial:
#>      B1   B2   B3
#> A1 3.35 2.01 6.99
#> A2 5.61 6.62 2.37
#> 
#> Average DLTs observed per trial:
#>      B1   B2   B3
#> A1 0.09 0.17 1.97
#> A2 0.50 2.08 1.20
#> 
#> Reason for ending (%):
#> 
#> completed 
#>       100 
```
