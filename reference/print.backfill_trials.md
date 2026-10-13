# Print Simulated Trials of a BOIN Design with Backfilling

Summarize the trials returned by
[`bf_boin_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/bf_boin_simulate.md)
or
[`be_boin_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/be_boin_simulate.md):
the design, the average numbers of patients, backfilled patients, DLTs
and responses at each dose, the trial duration, the waiting of the dose
escalation and the reasons for stopping.

## Usage

``` r
# S3 method for class 'backfill_trials'
print(x, ...)
```

## Arguments

- x:

  An object of class `backfill_trials`.

- ...:

  Further arguments, currently ignored.

## Value

The object `x`, invisibly.

## See also

[`bf_boin_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/bf_boin_simulate.md),
[`be_boin_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/be_boin_simulate.md)

## Examples

``` r
trials <- bf_boin_simulate(
  target = 0.25,
  p_true = c(0.04, 0.12, 0.25, 0.43, 0.63),
  p_resp = c(0.10, 0.20, 0.30, 0.45, 0.58),
  n_cohort = 10, cohort_size = 3, window = 1, accrual_rate = 3,
  n_trials = 100, seed = 1
)
print(trials)
#> Simulated BF-BOIN (Zhao et al., 2024) trials
#>   trials          : 100
#>   doses           : 5
#>   target DLT rate : 0.25
#>   cohorts         : 10 of size 3
#>   escalation size : 30
#>   backfill cap    : 12 patients per dose
#>   window          : 1 (DLT), 1 (response)
#>   accrual         : 3 patients per unit of time (exponential)
#> 
#> Average per trial:
#>                        DL1   DL2   DL3  DL4  DL5
#> Patients treated      8.79 13.81 11.97 3.29 0.30
#>   of which backfilled 2.88  3.70  1.89 0.23 0.00
#> DLTs observed         0.44  1.76  2.96 1.56 0.18
#> Responses             1.01  2.90  3.36 1.45 0.17
#> 
#> Patients per trial    : mean 38.16, of which backfilled 8.7
#> Trials with backfill  : 92%
#> Trial duration        : mean 17.11, median 17.21
#> Escalation waiting    : 100% of trials, mean time 6.72 per trial
#> 
#> Stopping reason (%):
#> 
#> max_sample_size     n_earlystop 
#>              87              13 
```
