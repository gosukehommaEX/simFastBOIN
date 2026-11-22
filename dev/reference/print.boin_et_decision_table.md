# Print a BOIN-ET Decision Table

Display the decision table produced by
[`boin_et_decision_table`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/boin_et_decision_table.md)
as one matrix per number of patients, with the DLTs in the rows and the
responses in the columns.

## Usage

``` r
# S3 method for class 'boin_et_decision_table'
print(x, ...)
```

## Arguments

- x:

  An object of class `boin_et_decision_table`.

- ...:

  Further arguments, currently ignored.

## Value

The object `x`, invisibly.

## See also

[`boin_et_decision_table`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/boin_et_decision_table.md)

## Examples

``` r
print(boin_et_decision_table(target = 0.30, target_resp = 0.60, n = 3))
#> BOIN-ET decision table
#>   lambda1 = 0.14, lambda2 = 0.35, eta1 = 0.48
#> 
#> 3 patients
#>     Responses
#> DLTs 0     1     2  3 
#>    0 E     E     S  S 
#>    1 E/S/D E/S/D S  S 
#>    2 D     D     D  D 
#>    3 D*    D*    D* D*
#> 
#> E = escalate, S = stay, D = de-escalate, E/S/D = choose among the current and neighboring doses;
#> * = the dose is no longer admissible
```
