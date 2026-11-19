# Print the Operating Characteristics of the BOIN-ETC Design

Print the Operating Characteristics of the BOIN-ETC Design

## Usage

``` r
# S3 method for class 'boin_etc_oc'
print(
  x,
  digits = 1,
  percent = FALSE,
  kable = FALSE,
  kable_format = "pipe",
  ...
)
```

## Arguments

- x:

  An object of class `boin_etc_oc` from
  [`sim_boin_etc`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin_etc.md).

- digits:

  Integer scalar. Decimal places. Defaults to 1.

- percent:

  Logical scalar. Show the patients, DLTs and responses at each
  combination as percentages of their totals. Defaults to `FALSE`.

- kable:

  Logical scalar. Print the table with
  [`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html). Defaults
  to `FALSE`.

- kable_format:

  Character scalar. Format passed to
  [`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html). Defaults
  to `"pipe"`.

- ...:

  Not used.

## Value

`x`, invisibly.

## See also

[`sim_boin_etc`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin_etc.md)

## Examples

``` r
p_true <- matrix(c(0.05, 0.10, 0.10, 0.30, 0.30, 0.50), nrow = 2)
p_resp <- matrix(c(0.10, 0.25, 0.25, 0.60, 0.60, 0.65), nrow = 2)
oc <- sim_boin_etc(target = 0.30, target_resp = 0.60, p_true = p_true,
                   p_resp = p_resp, n_cohort = 8, cohort_size = 3,
                   n_trials = 100, seed = 1)
print(oc)
#> BOIN-ETC1 operating characteristics
#>   target DLT rate      : 0.3
#>   target response rate : 0.6
#>   dose levels          : 2 of drug A, 3 of drug B
#>   trials               : 100
#>   cohorts              : 8 of size 3
#>   max sample size      : 24
#>   subtrial end (n_stop): 6
#>   lambda1/lambda2/eta1 : 0.14 / 0.35 / 0.48
#>   start of a row       : diagonal
#>   rule 4               : paper
#>   utility scores       : 100, 40, 60, 0; ties lowest
#> 
#>                    Item Drug A   B1   B2   B3
#>       True DLT rate (%)     A1  5.0 10.0 30.0
#>                             A2 10.0 30.0 50.0
#>  True response rate (%)     A1 10.0 25.0 60.0
#>                             A2 25.0 60.0 65.0
#>        OBD selected (%)     A1  2.0 19.0 30.0
#>                             A2 17.0 29.0  3.0
#>        Patients treated     A1  3.1  4.2  4.1
#>                             A2  4.3  4.6  1.6
#>       Patients with DLT     A1  0.2  0.4  1.3
#>                             A2  0.5  1.4  0.8
#>  Patients with response     A1  0.4  1.2  2.6
#>                             A2  1.3  2.5  1.0
#> 
#> No OBD selected                      : 0%
#> Trials stopped early                 : 0%
#> Average number of patients           : 22
#> Average number of subtrials          : 2
#> 
#> True OBD (highest true utility)      : A2B2, A1B3
#>   Selected as the OBD                : 59%
#>   Patients treated there             : 39.8%
#> Combinations above the target DLT    : 1 of 6
#>   Selected as the OBD                : 3%
#>   Patients treated there             : 7.1%
print(oc, percent = TRUE)
#> BOIN-ETC1 operating characteristics
#>   target DLT rate      : 0.3
#>   target response rate : 0.6
#>   dose levels          : 2 of drug A, 3 of drug B
#>   trials               : 100
#>   cohorts              : 8 of size 3
#>   max sample size      : 24
#>   subtrial end (n_stop): 6
#>   lambda1/lambda2/eta1 : 0.14 / 0.35 / 0.48
#>   start of a row       : diagonal
#>   rule 4               : paper
#>   utility scores       : 100, 40, 60, 0; ties lowest
#> 
#>                        Item Drug A   B1   B2   B3
#>           True DLT rate (%)     A1  5.0 10.0 30.0
#>                                 A2 10.0 30.0 50.0
#>      True response rate (%)     A1 10.0 25.0 60.0
#>                                 A2 25.0 60.0 65.0
#>            OBD selected (%)     A1  2.0 19.0 30.0
#>                                 A2 17.0 29.0  3.0
#>        Patients treated (%)     A1 14.1 19.2 18.8
#>                                 A2 19.8 21.0  7.1
#>       Patients with DLT (%)     A1  3.7  9.4 28.1
#>                                 A2 10.0 30.9 17.9
#>  Patients with response (%)     A1  4.0 13.3 28.8
#>                                 A2 14.3 28.1 11.6
#> 
#> No OBD selected                      : 0%
#> Trials stopped early                 : 0%
#> Average number of patients           : 22
#> Average number of subtrials          : 2
#> 
#> True OBD (highest true utility)      : A2B2, A1B3
#>   Selected as the OBD                : 59%
#>   Patients treated there             : 39.8%
#> Combinations above the target DLT    : 1 of 6
#>   Selected as the OBD                : 3%
#>   Patients treated there             : 7.1%
```
