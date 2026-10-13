# Print Operating Characteristics of the 3+3 Design

Display the summary table produced by
[`oc_3p3`](https://gosukehommaex.github.io/simFastBOIN/reference/oc_3p3.md)
or
[`sim_3p3`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_3p3.md).
The layout matches
[`print.boin_oc`](https://gosukehommaex.github.io/simFastBOIN/reference/print.boin_oc.md),
so a 3+3 table can be read beside a BOIN one.

## Usage

``` r
# S3 method for class 'oc_3p3'
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

  An object of class `oc_3p3`.

- digits:

  Integer scalar. Number of decimal places. Defaults to 1.

- percent:

  Logical scalar. Express the average number of patients and DLTs at
  each dose as a percentage of the trial total instead of a count.
  Defaults to `FALSE`.

- kable:

  Logical scalar. Return the table through
  [`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html) rather
  than printing it directly, which is convenient inside a report.
  Requires the knitr package. Defaults to `FALSE`.

- kable_format:

  Character scalar passed to
  [`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html), for
  example `"pipe"`, `"html"` or `"latex"`. Defaults to `"pipe"`.

- ...:

  Further arguments, currently ignored.

## Value

The object `x`, invisibly.

## See also

[`oc_3p3`](https://gosukehommaex.github.io/simFastBOIN/reference/oc_3p3.md),
[`sim_3p3`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_3p3.md)

## Examples

``` r
oc <- oc_3p3(p_true = c(0.30, 0.48, 0.67))

print(oc)
#> 3+3 operating characteristics
#>   MTD rule        : previous (dose below the toxic one)
#>   obtained by     : exact enumeration
#>   start dose      : DL1
#>   max sample size : 18 (six per dose)
#> 
#>                    DL1  DL2  DL3 Total / No MTD
#> True DLT rate (%) 30.0 48.0 67.0               
#> MTD selected (%)  39.8  9.2  0.4           50.6
#> Patients treated   4.3  2.1  0.4            6.7
#> Patients with DLT  1.3  1.0  0.2            2.5
#> 
#> Doses above a true DLT rate of 0.3333: DL2, DL3
#>   Patients treated there              : 35.8%
#>   Trials treating any patient there   : 49.4%
#>   Trials selecting an MTD there       : 9.7%
print(oc, percent = TRUE)
#> 3+3 operating characteristics
#>   MTD rule        : previous (dose below the toxic one)
#>   obtained by     : exact enumeration
#>   start dose      : DL1
#>   max sample size : 18 (six per dose)
#> 
#>                        DL1  DL2  DL3 Total / No MTD
#> True DLT rate (%)     30.0 48.0 67.0               
#> MTD selected (%)      39.8  9.2  0.4           50.6
#> Patients treated (%)  64.2 30.6  5.2            6.7
#> Patients with DLT (%) 51.4 39.2  9.4            2.5
#> 
#> Doses above a true DLT rate of 0.3333: DL2, DL3
#>   Patients treated there              : 35.8%
#>   Trials treating any patient there   : 49.4%
#>   Trials selecting an MTD there       : 9.7%
```
