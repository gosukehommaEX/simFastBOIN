# Decision Table for the 3+3 Design

Tabulate the dose assignment rule of the traditional 3+3 design used by
[`oc_3p3`](https://gosukehommaex.github.io/simFastBOIN/reference/oc_3p3.md)
and
[`sim_3p3`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_3p3.md):
for every number of DLTs among the three or six patients treated at the
current dose, the decision taken during dose escalation and, for
`mtd_rule = "expand"`, during the search for the MTD that follows it.

## Usage

``` r
decision_table_3p3(mtd_rule = c("previous", "expand"))
```

## Arguments

- mtd_rule:

  Character scalar, either `"previous"` or `"expand"`. The definition of
  the MTD, as in
  [`oc_3p3`](https://gosukehommaex.github.io/simFastBOIN/reference/oc_3p3.md).
  Defaults to `"previous"`.

## Value

An object of class `decision_table_3p3`, which is a data frame with one
row per state of the current dose and the columns

- stage:

  `"escalation"` for dose escalation, or `"search"` for the search for
  the MTD, which only `mtd_rule = "expand"` has.

- n:

  Number of patients treated at the current dose, 3 or 6.

- n_tox:

  Number of patients with a DLT, from 0 to `n`.

- decision:

  Decision code. During escalation, `"E"` escalates to the next higher
  dose, `"S"` treats three more patients at the current dose and
  `"STOP"` stops escalation, this dose and all higher doses being too
  toxic. During the search, `"S"` treats three more patients at the
  current dose, `"MTD"` selects it as the MTD and `"D"` moves the search
  to the next lower dose.

The definition of the MTD is stored in the attribute `mtd_rule`. The
object has `print` and `plot` methods.

## Details

During escalation three patients are treated at the current dose. No DLT
escalates, one DLT adds three more patients at the same dose, and two or
more stop escalation. With six patients, at most one DLT escalates and
two or more stop escalation. Escalating from the highest dose also ends
escalation.

With `mtd_rule = "previous"` the MTD is the dose below the one at which
escalation stopped, or the highest dose when it was escalated from, and
no MTD is selected when escalation stops at the starting dose. With
`mtd_rule = "expand"` the search starts at the dose below the one at
which escalation stopped, or at the highest dose when it was escalated
from. A dose with three patients receives three more; a dose with at
most one DLT out of six is selected as the MTD, and one with two or more
passes the search to the next lower dose. No MTD is selected when the
search moves below the starting dose.

Every number of DLTs from 0 to `n` is listed for six patients, although
in a trial the six patients at a dose during escalation always include
one DLT among the first three. The search examines only doses that were
escalated from, so a dose with three patients has no DLT there, and only
that state is listed.

The rules are those of
[`oc_3p3`](https://gosukehommaex.github.io/simFastBOIN/reference/oc_3p3.md):
the tests follow this table through every possible course of a trial and
obtain the operating characteristics of
[`oc_3p3()`](https://gosukehommaex.github.io/simFastBOIN/reference/oc_3p3.md).

## See also

[`print.decision_table_3p3`](https://gosukehommaex.github.io/simFastBOIN/reference/print.decision_table_3p3.md),
[`plot.decision_table_3p3`](https://gosukehommaex.github.io/simFastBOIN/reference/plot.decision_table_3p3.md),
[`oc_3p3`](https://gosukehommaex.github.io/simFastBOIN/reference/oc_3p3.md),
[`boin_decision_table`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_decision_table.md)

## Examples

``` r
decision_table_3p3()
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

# The other definition of the MTD adds the search for the MTD
decision_table_3p3(mtd_rule = "expand")
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
