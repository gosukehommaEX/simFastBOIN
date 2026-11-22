# Print the Desirability Score Table of the BOIN12 Design

Display the table of
[`boin12_decision_table`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/boin12_decision_table.md)
with the score of every outcome, `"E"` marking the outcomes that make
the dose inadmissible.

## Usage

``` r
# S3 method for class 'boin12_decision_table'
print(x, desirability = FALSE, ...)
```

## Arguments

- x:

  An object of class `boin12_decision_table`.

- desirability:

  Logical scalar. Show the posterior probabilities as well. Defaults to
  `FALSE`.

- ...:

  Not used.

## Value

`x`, invisibly.

## Examples

``` r
tab <- boin12_decision_table(target = 0.35, target_resp = 0.25, n = c(0, 3))
print(tab, desirability = TRUE)
#> BOIN12 desirability score table
#>   DLT limit / response limit : 0.35 / 0.25
#>   lambda_e / lambda_d        : 0.276 / 0.419
#>   utility scores             : 100 (response only), 40 (neither), 60 (both), 0 (DLT only)
#>   utility benchmark          : 70.5
#>   E: inadmissible (DLT cutoff 0.95, response cutoff 0.9)
#> 
#>  n n_tox n_resp score desirability
#>  0     0      0     8       0.2950
#>  3     0      0     4       0.1134
#>  3     0      1     7       0.2691
#>  3     0      2    11       0.5009
#>  3     0      3    13       0.7530
#>  3     1      0     2       0.0545
#>  3     1      1     5       0.1558
#>  3     1      2     9       0.3395
#>  3     1      3    12       0.5868
#>  3     2      0     1       0.0224
#>  3     2      1     3       0.0800
#>  3     2      2     6       0.2076
#>  3     2      3    10       0.4175
#>  3     3      0     E       0.0076
#>  3     3      1     E       0.0357
#>  3     3      2     E       0.1134
#>  3     3      3     E       0.2691
```
