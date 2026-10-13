# Print a 3+3 Decision Table

Display the decision table produced by
[`decision_table_3p3`](https://gosukehommaex.github.io/simFastBOIN/reference/decision_table_3p3.md):
one row per number of patients and range of DLTs sharing a decision, for
dose escalation and, with `mtd_rule = "expand"`, for the search for the
MTD, followed by what happens once escalation stops.

## Usage

``` r
# S3 method for class 'decision_table_3p3'
print(x, ...)
```

## Arguments

- x:

  An object of class `decision_table_3p3`.

- ...:

  Further arguments, currently ignored.

## Value

The object `x`, invisibly.

## Details

A data frame that no longer has the columns or attributes of a decision
table, for example after selecting columns, is printed as a plain data
frame.

## See also

[`decision_table_3p3`](https://gosukehommaex.github.io/simFastBOIN/reference/decision_table_3p3.md),
[`plot.decision_table_3p3`](https://gosukehommaex.github.io/simFastBOIN/reference/plot.decision_table_3p3.md)

## Examples

``` r
print(decision_table_3p3())
#> 3+3 decision table, MTD rule: previous (the dose below the one declared too toxic)
#> 
#> Dose escalation
#>  Patients DLTs Decision                          
#>  3        0    Escalate                          
#>  3        1    Treat 3 more patients at this dose
#>  3        >= 2 Stop escalation                   
#>  6        <= 1 Escalate                          
#>  6        >= 2 Stop escalation                   
#> 
#> Stop escalation means that this dose and all higher doses are too
#> toxic. The MTD is the dose below the one at which escalation stopped,
#> and no MTD is selected when escalation stops at the starting dose. When
#> the highest dose is escalated from, it is the MTD.

print(decision_table_3p3(mtd_rule = "expand"))
#> 3+3 decision table, MTD rule: expand (the highest dose with at most 1 DLT in 6)
#> 
#> Dose escalation
#>  Patients DLTs Decision                          
#>  3        0    Escalate                          
#>  3        1    Treat 3 more patients at this dose
#>  3        >= 2 Stop escalation                   
#>  6        <= 1 Escalate                          
#>  6        >= 2 Stop escalation                   
#> 
#> Search for the MTD
#>  Patients DLTs Decision                          
#>  3        0    Treat 3 more patients at this dose
#>  6        <= 1 Select this dose as the MTD       
#>  6        >= 2 Move to the next lower dose       
#> 
#> Stop escalation means that this dose and all higher doses are too
#> toxic. The search for the MTD starts at the dose below the one at which
#> escalation stopped, or at the highest dose when it is escalated from.
#> No MTD is selected when the search moves below the starting dose.
```
