# Decision Table for the BOIN Design

Build the lookup table that maps a pair of DLT and patient counts at the
current dose to the dose assignment decision.

## Usage

``` r
boin_decision_table(
  target,
  max_n,
  p_saf = NULL,
  p_tox = NULL,
  cutoff_eli = 0.95,
  stay_on_1_of_3 = FALSE
)
```

## Arguments

- target:

  Numeric scalar. Target DLT probability.

- max_n:

  Integer scalar. Largest number of patients treated at a single dose.

- p_saf:

  Numeric scalar. Highest DLT probability deemed subtherapeutic.
  Defaults to `0.6 * target`.

- p_tox:

  Numeric scalar. Lowest DLT probability deemed overly toxic. Defaults
  to `1.4 * target`.

- cutoff_eli:

  Numeric scalar. Posterior probability cutoff for dose elimination.
  Defaults to 0.95.

- stay_on_1_of_3:

  Logical scalar. When `TRUE`, one DLT out of three patients leads to
  staying at the current dose rather than de-escalating. Defaults to
  `FALSE`. See
  [`boin_boundary`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_boundary.md)
  for what the modification does.

## Value

An object of class `boin_decision_table`, which is a character matrix
with `max_n + 1` rows, labeled by the number of DLTs from 0 to `max_n`,
and `max_n` columns, labeled by the number of patients from 1 to
`max_n`. Entries are `"E"` to escalate, `"S"` to stay, `"D"` to
de-escalate and `"DE"` to de-escalate and eliminate the dose together
with all higher doses. Combinations with more DLTs than patients are
`NA`. The object has `print` and `plot` methods; subsetting it returns a
plain matrix.

## Details

The table is derived from
[`boin_boundary`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_boundary.md)
and therefore agrees with the decisions taken by the simulation engine.
Elimination takes precedence over the other rules, then escalation, then
de-escalation.

## See also

[`boin_boundary`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_boundary.md),
[`plot.boin_decision_table`](https://gosukehommaex.github.io/simFastBOIN/reference/plot.boin_decision_table.md)

## Examples

``` r
decisions <- boin_decision_table(target = 0.30, max_n = 12)
decisions
#>     Patients
#> DLTs 1 2 3  4  5  6  7  8  9  10 11 12
#>   0  E E E  E  E  E  E  E  E  E  E  E 
#>   1  D D S  S  E  E  E  E  E  E  E  E 
#>   2    D D  D  D  S  S  S  E  E  E  E 
#>   3      DE DE D  D  D  D  S  S  S  S 
#>   4         DE DE DE D  D  D  D  D  S 
#>   5            DE DE DE DE DE D  D  D 
#>   6               DE DE DE DE DE DE D 
#>   7                  DE DE DE DE DE DE
#>   8                     DE DE DE DE DE
#>   9                        DE DE DE DE
#>   10                          DE DE DE
#>   11                             DE DE
#>   12                                DE
#> 
#> E = escalate, S = stay, D = de-escalate, DE = de-escalate and eliminate this dose and above

# Only the sample sizes reached at the end of a cohort of three
print(decisions, cohort_size = 3)
#>     Patients
#> DLTs 3  6  9  12
#>   0  E  E  E  E 
#>   1  S  E  E  E 
#>   2  D  S  E  E 
#>   3  DE D  S  S 
#>   4     DE D  S 
#>   5     DE DE D 
#>   6     DE DE D 
#>   7        DE DE
#>   8        DE DE
#>   9        DE DE
#>   10          DE
#>   11          DE
#>   12          DE
#> 
#> E = escalate, S = stay, D = de-escalate, DE = de-escalate and eliminate this dose and above

# One DLT out of three de-escalates at a target of 0.25
boin_decision_table(target = 0.25, max_n = 9)["1", "3"]
#> [1] "D"

# and stays once the modification is switched on
boin_decision_table(target = 0.25, max_n = 9, stay_on_1_of_3 = TRUE)["1", "3"]
#> [1] "S"
```
