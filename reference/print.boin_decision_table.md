# Print a BOIN Decision Table

Display the decision table produced by
[`boin_decision_table`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_decision_table.md),
optionally restricted to the sample sizes reached at the end of a
cohort.

## Usage

``` r
# S3 method for class 'boin_decision_table'
print(x, cohort_size = NULL, ...)
```

## Arguments

- x:

  An object of class `boin_decision_table`.

- cohort_size:

  Integer scalar or `NULL`. When supplied, only the sample sizes that
  are multiples of `cohort_size` are shown, which is how the table is
  consulted during a trial with equally sized cohorts.

- ...:

  Further arguments, currently ignored.

## Value

The object `x`, invisibly.

## See also

[`boin_decision_table`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_decision_table.md),
[`plot.boin_decision_table`](https://gosukehommaex.github.io/simFastBOIN/reference/plot.boin_decision_table.md)

## Examples

``` r
decisions <- boin_decision_table(target = 0.30, max_n = 18)

print(decisions, cohort_size = 3)
#>     Patients
#> DLTs 3  6  9  12 15 18
#>   0  E  E  E  E  E  E 
#>   1  S  E  E  E  E  E 
#>   2  D  S  E  E  E  E 
#>   3  DE D  S  S  E  E 
#>   4     DE D  S  S  E 
#>   5     DE DE D  S  S 
#>   6     DE DE D  D  S 
#>   7        DE DE D  D 
#>   8        DE DE DE D 
#>   9        DE DE DE DE
#>   10          DE DE DE
#>   11          DE DE DE
#>   12          DE DE DE
#>   13             DE DE
#>   14             DE DE
#>   15             DE DE
#>   16                DE
#>   17                DE
#>   18                DE
#> 
#> E = escalate, S = stay, D = de-escalate, DE = de-escalate and eliminate this dose and above
```
