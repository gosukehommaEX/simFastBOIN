# Print the Decision Boundaries of the iBOIN Design

Print the prior information and, for every dose, the numbers of DLTs
that lead to escalation and de-escalation, in the layout of Table 1 of
Zhou et al. (2021).

## Usage

``` r
# S3 method for class 'iboin_boundary'
print(x, cohort_size = NULL, ...)
```

## Arguments

- x:

  An object of class `iboin_boundary`, from
  [`iboin_boundary`](https://gosukehommaex.github.io/simFastBOIN/reference/iboin_boundary.md).

- cohort_size:

  Integer scalar or `NULL`. When given, only the numbers of patients
  that are multiples of the cohort size are shown. Defaults to `NULL`,
  every number of patients.

- ...:

  Unused.

## Value

The object `x`, invisibly.

## Details

An escalation boundary printed as `NA` means that no number of DLTs
leads to escalation at that dose and sample size, and a de-escalation
boundary printed as `NA` that none leads to de-escalation.

## See also

[`iboin_boundary`](https://gosukehommaex.github.io/simFastBOIN/reference/iboin_boundary.md)

## Examples

``` r
bd <- iboin_boundary(
  target = 0.30, skeleton = c(0.10, 0.19, 0.30, 0.42, 0.54),
  prior_ess = 3, max_n = 30
)
print(bd, cohort_size = 3)
#> iBOIN decision boundaries
#>   target DLT rate : 0.3
#>   p_saf / p_tox   : 0.18 / 0.42
#>   cutoff_eli      : 0.95
#>   robust prior    : no
#> 
#> Prior information:
#>     Skeleton Prior ESS Pr(target) Pr(p_saf) Pr(p_tox)
#> DL1     0.10         3      0.326     0.445     0.230
#> DL2     0.19         3      0.332     0.392     0.276
#> DL3     0.30         3      0.335     0.332     0.333
#> DL4     0.42         3      0.333     0.271     0.396
#> DL5     0.54         3      0.326     0.215     0.459
#> 
#> Number of patients treated at the current dose (columns):
#>                                3 6 9 12 15 18 21 24 27 30
#> DL1 escalate if # of DLT <=    1 1 2  3  4  4  5  6  6  7
#> DL1 de-escalate if # of DLT >= 2 3 4  5  7  8  9 10 11 12
#> DL2 escalate if # of DLT <=    0 1 2  3  3  4  5  5  6  7
#> DL2 de-escalate if # of DLT >= 2 3 4  5  6  7  8  9 11 12
#> DL3 escalate if # of DLT <=    0 1 2  2  3  4  4  5  6  7
#> DL3 de-escalate if # of DLT >= 2 3 4  5  6  7  8  9 10 11
#> DL4 escalate if # of DLT <=    0 1 1  2  3  3  4  5  6  6
#> DL4 de-escalate if # of DLT >= 1 2 3  4  6  7  8  9 10 11
#> DL5 escalate if # of DLT <=    0 0 1  2  2  3  4  5  5  6
#> DL5 de-escalate if # of DLT >= 1 2 3  4  5  6  7  8 10 11
#> Eliminate if # of DLT >=       3 4 5  7  8  9 10 11 12 14
```
