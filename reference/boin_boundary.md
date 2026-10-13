# Integer Decision Boundaries for the BOIN Design

Tabulate, for every attainable number of patients treated at a dose, the
number of DLTs that triggers escalation, de-escalation, elimination and,
when requested, the stricter safety stopping rule at the lowest dose.

## Usage

``` r
boin_boundary(
  target,
  max_n,
  p_saf = NULL,
  p_tox = NULL,
  cutoff_eli = 0.95,
  extrasafe = FALSE,
  offset = 0.05,
  stay_on_1_of_3 = FALSE
)
```

## Arguments

- target:

  Numeric scalar. Target DLT probability, for example 0.30.

- max_n:

  Integer scalar. Largest number of patients that can be treated at a
  single dose, which is the number of columns of the resulting table.
  For a trial with `n_cohort` cohorts of size `cohort_size` this is
  `n_cohort * cohort_size`.

- p_saf:

  Numeric scalar. Highest DLT probability deemed subtherapeutic.
  Defaults to `0.6 * target`.

- p_tox:

  Numeric scalar. Lowest DLT probability deemed overly toxic. Defaults
  to `1.4 * target`.

- cutoff_eli:

  Numeric scalar. Posterior probability cutoff above which a dose is
  eliminated for toxicity. Defaults to 0.95.

- extrasafe:

  Logical scalar. When `TRUE`, also tabulate the safety stopping
  boundary applied at the lowest dose. Defaults to `FALSE`.

- offset:

  Numeric scalar between 0 and 0.5. Amount by which `cutoff_eli` is
  relaxed for the safety stopping rule. Defaults to 0.05.

- stay_on_1_of_3:

  Logical scalar. When `TRUE`, one DLT out of three patients leads to
  staying at the current dose rather than de-escalating. Defaults to
  `FALSE`. See the details.

## Value

An object of class `boin_boundary`, which is a list with components

- lambda_e:

  Escalation interval boundary.

- lambda_d:

  De-escalation interval boundary.

- n:

  Number of patients, 1 to `max_n`.

- b_esc:

  Escalate when the number of DLTs is at most this value.

- b_deesc:

  De-escalate when the number of DLTs reaches this value.

- b_elim:

  Eliminate the dose when the number of DLTs reaches this value, `NA`
  when no elimination is possible at that sample size.

- b_stop:

  Safety stopping boundary at the lowest dose, all `NA` unless
  `extrasafe` is `TRUE`.

- stay_on_1_of_3:

  Whether the modification was requested.

- stay_on_1_of_3_applied:

  Whether it changed the boundaries.

together with the design parameters used.

## Details

The escalation boundary is `floor(lambda_e * n)` and the de-escalation
boundary is `ceiling(lambda_d * n)`, with the convention that when
`lambda_d * n` is a whole number de-escalation requires one more DLT.
The elimination boundary is the smallest positive number of DLTs for
which `Pr(p > target | data)` exceeds `cutoff_eli` under a uniform
Beta(1, 1) prior, and is not evaluated before three patients have been
treated. The de-escalation boundary is capped at the elimination
boundary.

With the default thresholds these definitions reproduce
[`BOIN::get.boundary()`](https://rdrr.io/pkg/BOIN/man/get.boundary.html)
exactly.

## Staying on one DLT out of three

For some target rates the optimal BOIN decision after one DLT in three
patients is to de-escalate, because the likelihood of overdosing then
exceeds the likelihood of proper dosing. Long practice with the 3+3
design has nonetheless made staying at the dose the widely accepted
choice, and `stay_on_1_of_3 = TRUE` aligns the design with that practice
by raising the de-escalation boundary at three patients from one DLT to
two. Nothing else in the table changes.

The modification is applied only where it is meaningful, that is when
one DLT out of three currently triggers de-escalation. It never
overrides an escalation, and it never overrides an elimination, which is
a safety rule. Whether it took effect is reported in
`stay_on_1_of_3_applied`. With the default thresholds it takes effect
for target rates from about 0.098 to 0.279; below that range one DLT out
of three already eliminates the dose, and above it the design already
stays.

## References

Liu S. and Yuan, Y. (2015). Bayesian Optimal Interval Designs for Phase
I Clinical Trials. Journal of the Royal Statistical Society: Series C,
64, 507-523.

## See also

[`boin_lambda`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_lambda.md),
[`boin_decision_table`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_decision_table.md)

## Examples

``` r
bd <- boin_boundary(target = 0.30, max_n = 18)
bd
#> BOIN decision boundaries
#>   target DLT rate : 0.3
#>   p_saf / p_tox   : 0.18 / 0.42
#>   lambda_e        : 0.2365
#>   lambda_d        : 0.3585
#>   cutoff_eli      : 0.95
#> 
#>                                                                          
#> Number of patients treated  1  2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18
#> Escalate if # of DLT <=     0  0 0 0 1 1 1 1 2  2  2  2  3  3  3  3  4  4
#> Deescalate if # of DLT >=   1  1 2 2 2 3 3 3 4  4  4  5  5  6  6  6  7  7
#> Eliminate if # of DLT >=   NA NA 3 3 4 4 5 5 5  6  6  7  7  8  8  8  9  9

# At a target of 0.25 one DLT out of three de-escalates by default
boin_boundary(target = 0.25, max_n = 18)$b_deesc[3]
#> [1] 1

# and stays once the modification is switched on
modified <- boin_boundary(target = 0.25, max_n = 18, stay_on_1_of_3 = TRUE)
modified$b_deesc[3]
#> [1] 2
modified$stay_on_1_of_3_applied
#> [1] TRUE

# At a target of 0.30 the design already stays, so nothing changes
boin_boundary(target = 0.30, max_n = 18, stay_on_1_of_3 = TRUE)$stay_on_1_of_3_applied
#> [1] FALSE
```
