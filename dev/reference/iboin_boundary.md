# Decision Boundaries of the iBOIN Design

Tabulate the escalation and de-escalation boundaries of the iBOIN design
of Zhou et al. (2021), which brings prior information on the DLT
probability of every dose into the BOIN design. The prior information is
a skeleton, the prior DLT probability of each dose, and a prior
effective sample size (PESS), the number of patients the prior is worth.
The boundaries depend on the dose and on the number of patients treated
there.

## Usage

``` r
iboin_boundary(
  target,
  skeleton,
  prior_ess,
  max_n,
  p_saf = NULL,
  p_tox = NULL,
  cutoff_eli = 0.95,
  extrasafe = FALSE,
  offset = 0.05,
  robust = FALSE
)
```

## Arguments

- target:

  Numeric scalar. Target DLT probability, for example 0.30.

- skeleton:

  Numeric vector. Prior DLT probability at each dose, strictly between 0
  and 1.

- prior_ess:

  Whole number or vector of whole numbers, one per dose. Prior effective
  sample size. A value of 0 gives the boundaries of the BOIN design at
  that dose.

- max_n:

  Integer scalar. Largest number of patients that can be treated at a
  single dose, which is the number of columns of the resulting tables.

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

- robust:

  Logical scalar. Use the robust prior of the article: when the prior
  MTD, the dose whose skeleton value is closest to `target`, has a
  number of at least half the number of doses, the doses above it get a
  prior effective sample size of 0. Defaults to `FALSE`.

## Value

An object of class `iboin_boundary`, which is a list with components

- prior:

  Matrix with one row per dose: the prior probabilities that the DLT
  probability equals `target`, `p_saf` and `p_tox`.

- lambda_e, lambda_d:

  Matrices, doses by `1:max_n` patients, of the escalation and
  de-escalation boundaries of the DLT rate.

- n:

  Number of patients, 1 to `max_n`.

- b_esc:

  Integer matrix. Escalate when the number of DLTs is at most this
  value; -1 when no number of DLTs escalates.

- b_deesc:

  Integer matrix. De-escalate when the number of DLTs reaches this
  value; a value above the number of patients means that no number of
  DLTs de-escalates.

- b_elim:

  Eliminate the dose when the number of DLTs reaches this value, `NA`
  when no elimination is possible at that sample size.

- b_stop:

  Safety stopping boundary at the lowest dose, all `NA` unless
  `extrasafe` is `TRUE`.

- searched:

  Logical matrix marking the boundaries found by search rather than by
  the closed form, see the details.

- prior_ess:

  The prior effective sample sizes used, after the robust prior has been
  applied.

- prior_mtd:

  The dose whose skeleton value is closest to `target`.

together with the design parameters used.

## Details

The BOIN design chooses its boundaries to minimize the probability of an
incorrect decision under three point hypotheses: the DLT probability of
the current dose equals `target` (stay), `p_saf` (escalate) or `p_tox`
(de-escalate). With equal prior probabilities of the three hypotheses
the boundaries are those of
[`boin_lambda`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/boin_lambda.md).
The iBOIN design derives the prior probabilities of the hypotheses at
each dose from the skeleton and the PESS (equation 5 of the article):
the posterior probabilities of the hypotheses given `x` DLTs among
`prior_ess` patients, averaged over the binomial distribution of `x`
with the skeleton value as DLT probability. The boundaries then follow
from equation 3 of the article, in which the prior probabilities shift
the BOIN boundaries by an amount that decreases as one over the number
of patients treated at the dose, so that the data override the prior as
they accumulate.

The design escalates when the observed DLT rate is at most `lambda_e`
and de-escalates when it exceeds `lambda_d`; the integer boundaries are
cut as in
[`boin_boundary`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/boin_boundary.md).
With the skeleton (0.10, 0.19, 0.30, 0.42, 0.54), a PESS of 3 at every
dose and a target of 0.30 they reproduce Table 1 of the article.

Where the closed form gives `lambda_e > lambda_d`, the article
determines the boundaries by a numerical search. Here the search is over
the integer boundaries `e < d` themselves: the pair that minimizes the
probability of an incorrect decision is used, and among pairs with the
same probability the smallest `e` and then the smallest `d`. This
happens mostly for a single patient at a dose.

The elimination boundary and the safety stopping boundary are those of
the BOIN design: they use a uniform prior, not the skeleton, so that a
misspecified prior cannot weaken the safety rules. As in the BOIN design
the de-escalation boundary is capped at the elimination boundary. The
escalation boundary is kept below the de-escalation boundary, so that a
dose that has just been eliminated is never followed by an escalation.

The PESS must be a whole number, because the prior probabilities sum
over the possible numbers of DLTs among `prior_ess` patients. The
article recommends a PESS between one third and one half of the number
of patients per dose, for example 2 or 3 for 30 patients and 5 doses.

## References

Zhou, Y., Lee, J. J., Wang, S., Bailey, S. and Yuan, Y. (2021).
Incorporating historical information to improve phase I clinical trials.
Pharmaceutical Statistics, 20, 1017-1034.
[doi:10.1002/pst.2121](https://doi.org/10.1002/pst.2121)

Liu S. and Yuan, Y. (2015). Bayesian Optimal Interval Designs for Phase
I Clinical Trials. Journal of the Royal Statistical Society: Series C,
64, 507-523.

## See also

[`sim_iboin`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_iboin.md),
[`boin_boundary`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/boin_boundary.md)

## Examples

``` r
# Table 1 of Zhou et al. (2021)
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

# Without prior information the boundaries are those of the BOIN design
iboin_boundary(0.30, c(0.10, 0.19, 0.30, 0.42, 0.54), 0, 9)$b_esc
#>     1 2 3 4 5 6 7 8 9
#> DL1 0 0 0 0 1 1 1 1 2
#> DL2 0 0 0 0 1 1 1 1 2
#> DL3 0 0 0 0 1 1 1 1 2
#> DL4 0 0 0 0 1 1 1 1 2
#> DL5 0 0 0 0 1 1 1 1 2
boin_boundary(0.30, 9)$b_esc
#> [1] 0 0 0 0 1 1 1 1 2
```
