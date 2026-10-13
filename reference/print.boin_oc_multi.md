# Print Operating Characteristics Across Scenarios

Display the combined summary table produced by
[`sim_boin_multi`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin_multi.md).

## Usage

``` r
# S3 method for class 'boin_oc_multi'
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

  An object of class `boin_oc_multi`.

- digits:

  Integer scalar. Number of decimal places. Defaults to 1.

- percent:

  Logical scalar. Express the average number of patients and DLTs at
  each dose as a percentage of the trial total instead of a count.
  Defaults to `FALSE`.

- kable:

  Logical scalar. Return the table through
  [`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html) rather
  than printing it directly. Requires the knitr package. Defaults to
  `FALSE`.

- kable_format:

  Character scalar passed to
  [`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html), for
  example `"pipe"`, `"html"` or `"latex"`. Defaults to `"pipe"`.

- ...:

  Further arguments, currently ignored.

## Value

The object `x`, invisibly.

## Examples

``` r
oc <- sim_boin_multi(
  target = 0.30,
  scenarios = list(
    list(name = "MTD at dose 3", p_true = c(0.05, 0.15, 0.30, 0.45, 0.60)),
    list(name = "All toxic",     p_true = c(0.35, 0.45, 0.55, 0.65, 0.75))
  ),
  n_cohort = 10,
  cohort_size = 3,
  n_trials = 200,
  seed = 123
)

print(oc)
#> BOIN operating characteristics across 2 scenarios
#>   target DLT rate : 0.3
#>   trials each     : 200
#>   max sample size : 30
#> 
#>       Scenario              Item  DL1  DL2  DL3  DL4  DL5 Total / No MTD
#>  MTD at dose 3 True DLT rate (%)  5.0 15.0 30.0 45.0 60.0               
#>                 MTD selected (%)  1.0 26.0 52.0 21.0  0.0            0.0
#>                 Patients treated  4.2  9.1 10.6  4.5  0.6           29.1
#>                Patients with DLT  0.2  1.4  3.2  2.1  0.4            7.3
#>      All toxic True DLT rate (%) 35.0 45.0 55.0 65.0 75.0               
#>                 MTD selected (%) 67.0  9.0  2.0  0.0  0.0           22.0
#>                 Patients treated 14.5  5.4  0.8  0.1  0.0           20.8
#>                Patients with DLT  5.0  2.5  0.4  0.1  0.0            8.0
print(oc, percent = TRUE)
#> BOIN operating characteristics across 2 scenarios
#>   target DLT rate : 0.3
#>   trials each     : 200
#>   max sample size : 30
#> 
#>       Scenario                  Item  DL1  DL2  DL3  DL4  DL5 Total / No MTD
#>  MTD at dose 3     True DLT rate (%)  5.0 15.0 30.0 45.0 60.0               
#>                     MTD selected (%)  1.0 26.0 52.0 21.0  0.0            0.0
#>                 Patients treated (%) 14.4 31.4 36.5 15.6  2.2           29.1
#>                Patients with DLT (%)  3.2 19.2 43.9 28.3  5.4            7.3
#>      All toxic     True DLT rate (%) 35.0 45.0 55.0 65.0 75.0               
#>                     MTD selected (%) 67.0  9.0  2.0  0.0  0.0           22.0
#>                 Patients treated (%) 69.6 26.1  3.8  0.5  0.0           20.8
#>                Patients with DLT (%) 62.7 31.0  5.4  0.9  0.0            8.0
```
