# Print Simulated BOIN-ETC Trials

Print Simulated BOIN-ETC Trials

## Usage

``` r
# S3 method for class 'boin_etc_trials'
print(x, ...)
```

## Arguments

- x:

  An object of class `boin_etc_trials` from
  [`boin_etc_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_etc_simulate.md).

- ...:

  Not used.

## Value

`x`, invisibly.

## See also

[`boin_etc_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_etc_simulate.md)

## Examples

``` r
p_true <- matrix(c(0.05, 0.10, 0.10, 0.30, 0.30, 0.50), nrow = 2)
p_resp <- matrix(c(0.10, 0.25, 0.25, 0.60, 0.60, 0.65), nrow = 2)
trials <- boin_etc_simulate(target = 0.30, target_resp = 0.60,
                            p_true = p_true, p_resp = p_resp,
                            n_cohort = 8, cohort_size = 3,
                            n_trials = 100, seed = 1)
print(trials)
#> Simulated BOIN-ETC1 trials
#>   trials               : 100
#>   dose levels          : 2 of drug A, 3 of drug B
#>   target DLT rate      : 0.3
#>   target response rate : 0.6
#>   cohorts              : 8 of size 3
#>   max sample size      : 24
#>   subtrial end (n_stop): 6
#> 
#> Patients treated, average per trial:
#>      B1   B2   B3
#> A1 3.09 4.23 4.14
#> A2 4.35 4.62 1.56
#> 
#> DLTs observed, average per trial:
#>      B1   B2   B3
#> A1 0.17 0.43 1.29
#> A2 0.46 1.42 0.82
#> 
#> Responses observed, average per trial:
#>      B1   B2   B3
#> A1 0.36 1.19 2.59
#> A2 1.28 2.52 1.04
#> 
#> Average number of subtrials: 2
#> 
#> Stopping reason (%):
#> 
#>         max_cohorts subtrials_completed 
#>                  55                  45 
```
