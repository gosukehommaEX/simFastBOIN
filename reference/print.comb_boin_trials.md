# Print Simulated BOIN Drug-Combination Trials

Show a compact summary of the object returned by
[`comb_boin_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/comb_boin_simulate.md)
instead of printing the full trial arrays.

## Usage

``` r
# S3 method for class 'comb_boin_trials'
print(x, ...)
```

## Arguments

- x:

  An object of class `comb_boin_trials`.

- ...:

  Further arguments, currently ignored.

## Value

The object `x`, invisibly.

## Examples

``` r
p_true <- matrix(c(0.05, 0.15, 0.30,
                   0.15, 0.30, 0.45), nrow = 2, byrow = TRUE)

trials <- comb_boin_simulate(
  target = 0.30,
  p_true = p_true,
  n_cohort = 10,
  cohort_size = 3,
  n_trials = 100,
  seed = 1
)
print(trials)
#> Simulated BOIN drug-combination trials
#>   trials          : 100
#>   dose levels     : 2 of drug A, 3 of drug B
#>   target DLT rate : 0.3
#>   cohorts         : 10 of size 3
#>   max sample size : 30
#> 
#> Average patients treated per trial:
#>      B1   B2   B3
#> A1 4.26 5.16 3.48
#> A2 3.60 8.25 4.47
#> 
#> Average DLTs observed per trial:
#>      B1   B2   B3
#> A1 0.27 0.96 1.23
#> A2 0.48 2.49 2.25
#> 
#> Stopping reason (%):
#> 
#> max_cohorts n_earlystop 
#>          71          29 
```
