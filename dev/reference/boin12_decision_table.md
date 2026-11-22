# Desirability Score Table of the BOIN12 Design

Tabulate the desirability of a dose for every possible outcome of the
patients treated there, as in Table 3 of Lin et al. (2020): the
posterior probability that the utility of the dose exceeds the
benchmark, its rank-based desirability score (RDS) and whether the dose
is admissible.

## Usage

``` r
boin12_decision_table(
  target,
  target_resp,
  n = c(0, 3, 6, 9),
  p_saf = NULL,
  p_tox = NULL,
  cutoff_eli = 0.95,
  cutoff_eli_resp = 0.9,
  psi00 = 40,
  psi11 = 60,
  benchmark_weight = 0.5,
  utility_prior = c(1, 1)
)
```

## Arguments

- target:

  Numeric scalar. Highest acceptable DLT probability, called phi_T by
  Lin et al. (2020). It is the target of the BOIN boundaries `lambda_e`
  and `lambda_d`, of the admissibility of a dose and of the MTD that
  bounds the OBD. Lin et al. recommend a value slightly above the target
  of a design that only considers toxicity, for example 0.35 instead of
  0.30.

- target_resp:

  Numeric scalar. Lowest acceptable response probability, phi_E. Note
  that this is a lower limit, as in
  [`sim_bf_boin_et`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_bf_boin_et.md),
  whereas the `target_resp` of
  [`sim_boin_et`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_boin_et.md)
  is a target.

- n:

  Integer vector. Numbers of patients at a dose to tabulate. Defaults to
  `c(0, 3, 6, 9)`.

- p_saf, p_tox:

  Numeric scalars or `NULL`. Highest DLT probability deemed
  subtherapeutic and lowest deemed overly toxic, which give the BOIN
  boundaries. Default to `0.6 * target` and `1.4 * target`.

- cutoff_eli:

  Numeric scalar. A dose is inadmissible, together with every higher
  dose, when the posterior probability that its DLT rate exceeds
  `target` is at least this cutoff. Defaults to 0.95.

- cutoff_eli_resp:

  Numeric scalar. A dose is inadmissible when the posterior probability
  that its response rate is below `target_resp` is at least this cutoff.
  Defaults to 0.90.

- psi00, psi11:

  Numeric scalars from 0 to 100. Utility of a patient without DLT and
  without response, and with both. The utility of a response without DLT
  is 100 and that of a DLT without response is 0. Default to 40 and 60,
  the example of Lin et al. (2020).

- benchmark_weight:

  Numeric scalar above 0 and at most 1. The utility benchmark is
  `benchmark_weight * u_low + (1 - benchmark_weight) * 100`, where
  `u_low` is the expected utility of a dose whose DLT and response
  probabilities are `target` and `target_resp`, the two outcomes being
  independent. Defaults to 0.5, the benchmark of Lin et al. (2020); Lu
  et al. (2025) use 0.6.

- utility_prior:

  Numeric vector of two positive numbers. Beta prior of the utility of a
  dose divided by 100. Defaults to `c(1, 1)`.

## Value

A data frame of class `boin12_decision_table` with one row per outcome
and columns `n`, `n_tox`, `n_resp`, `n_tox_resp` (patients with both a
DLT and a response; only when `psi00 + psi11` is not 100, since the
utility otherwise depends on the margins alone), `admissible`,
`desirability` and `score`. The design parameters are kept in the
attribute `"design"`, a list with `lambda_e`, `lambda_d`, `u_low` and
`u_benchmark` among others.

## Details

The utility of a patient is 100 for a response without DLT, `psi00` for
neither, `psi11` for both and 0 for a DLT without response. The utility
of a dose, divided by 100, has a beta posterior under the quasi-binomial
model with the prior `utility_prior`, and the desirability is its
posterior probability of exceeding `u_benchmark / 100`. The score ranks
the desirability among the admissible outcomes of the table, from 1 for
the lowest, equal desirabilities sharing the lowest rank. A dose is
inadmissible (`"E"` in the printed table) when the posterior
probability, under a uniform prior, that its DLT rate exceeds `target`
is at least `cutoff_eli`, or that its response rate is below
`target_resp` is at least `cutoff_eli_resp`. During a trial only the
order of the desirabilities of the doses compared matters, so the score
gives the same decisions.

## References

Lin, R., Zhou, Y., Yan, F., Li, D. and Yuan, Y. (2020). BOIN12: Bayesian
optimal interval phase I/II trial design for utility-based dose finding
in immunotherapy and targeted therapies. JCO Precision Oncology, 4,
1393-1402.
[doi:10.1200/PO.20.00257](https://doi.org/10.1200/PO.20.00257)

## See also

[`boin12_simulate`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/boin12_simulate.md),
[`sim_boin12`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_boin12.md)

## Examples

``` r
# Table 3 of Lin et al. (2020)
tab <- boin12_decision_table(target = 0.35, target_resp = 0.25)
tab
#> BOIN12 desirability score table
#>   DLT limit / response limit : 0.35 / 0.25
#>   lambda_e / lambda_d        : 0.276 / 0.419
#>   utility scores             : 100 (response only), 40 (neither), 60 (both), 0 (DLT only)
#>   utility benchmark          : 70.5
#>   E: inadmissible (DLT cutoff 0.95, response cutoff 0.9)
#> 
#>  n n_tox n_resp score
#>  0     0      0    60
#>  3     0      0    35
#>  3     0      1    55
#>  3     0      2    76
#>  3     0      3    91
#>  3     1      0    24
#>  3     1      1    44
#>  3     1      2    63
#>  3     1      3    80
#>  3     2      0    13
#>  3     2      1    31
#>  3     2      2    48
#>  3     2      3    69
#>  3     3      0     E
#>  3     3      1     E
#>  3     3      2     E
#>  3     3      3     E
#>  6     0      0    22
#>  6     0      1    38
#>  6     0      2    51
#>  6     0      3    67
#>  6     0      4    81
#>  6     0      5    93
#>  6     0      6   100
#>  6     1      0    15
#>  6     1      1    27
#>  6     1      2    42
#>  6     1      3    56
#>  6     1      4    72
#>  6     1      5    87
#>  6     1      6    96
#>  6     2      0     8
#>  6     2      1    19
#>  6     2      2    34
#>  6     2      3    47
#>  6     2      4    64
#>  6     2      5    77
#>  6     2      6    90
#>  6     3      0     4
#>  6     3      1    12
#>  6     3      2    22
#>  6     3      3    38
#>  6     3      4    51
#>  6     3      5    67
#>  6     3      6    81
#>  6     4      0     1
#>  6     4      1     6
#>  6     4      2    15
#>  6     4      3    27
#>  6     4      4    42
#>  6     4      5    56
#>  6     4      6    72
#>  6     5      0     E
#>  6     5      1     E
#>  6     5      2     E
#>  6     5      3     E
#>  6     5      4     E
#>  6     5      5     E
#>  6     5      6     E
#>  6     6      0     E
#>  6     6      1     E
#>  6     6      2     E
#>  6     6      3     E
#>  6     6      4     E
#>  6     6      5     E
#>  6     6      6     E
#>  9     0      0     E
#>  9     0      1    25
#>  9     0      2    36
#>  9     0      3    49
#>  9     0      4    61
#>  9     0      5    74
#>  9     0      6    85
#>  9     0      7    94
#>  9     0      8    99
#>  9     0      9   102
#>  9     1      0     E
#>  9     1      1    17
#>  9     1      2    29
#>  9     1      3    40
#>  9     1      4    53
#>  9     1      5    65
#>  9     1      6    78
#>  9     1      7    88
#>  9     1      8    97
#>  9     1      9   101
#>  9     2      0     E
#>  9     2      1    10
#>  9     2      2    20
#>  9     2      3    32
#>  9     2      4    45
#>  9     2      5    58
#>  9     2      6    70
#>  9     2      7    83
#>  9     2      8    92
#>  9     2      9    98
#>  9     3      0     E
#>  9     3      1     7
#>  9     3      2    14
#>  9     3      3    25
#>  9     3      4    36
#>  9     3      5    49
#>  9     3      6    61
#>  9     3      7    74
#>  9     3      8    85
#>  9     3      9    94
#>  9     4      0     E
#>  9     4      1     3
#>  9     4      2     9
#>  9     4      3    17
#>  9     4      4    29
#>  9     4      5    40
#>  9     4      6    53
#>  9     4      7    65
#>  9     4      8    78
#>  9     4      9    88
#>  9     5      0     E
#>  9     5      1     2
#>  9     5      2     5
#>  9     5      3    10
#>  9     5      4    20
#>  9     5      5    32
#>  9     5      6    45
#>  9     5      7    58
#>  9     5      8    70
#>  9     5      9    83
#>  9     6      0     E
#>  9     6      1     E
#>  9     6      2     E
#>  9     6      3     E
#>  9     6      4     E
#>  9     6      5     E
#>  9     6      6     E
#>  9     6      7     E
#>  9     6      8     E
#>  9     6      9     E
#>  9     7      0     E
#>  9     7      1     E
#>  9     7      2     E
#>  9     7      3     E
#>  9     7      4     E
#>  9     7      5     E
#>  9     7      6     E
#>  9     7      7     E
#>  9     7      8     E
#>  9     7      9     E
#>  9     8      0     E
#>  9     8      1     E
#>  9     8      2     E
#>  9     8      3     E
#>  9     8      4     E
#>  9     8      5     E
#>  9     8      6     E
#>  9     8      7     E
#>  9     8      8     E
#>  9     8      9     E
#>  9     9      0     E
#>  9     9      1     E
#>  9     9      2     E
#>  9     9      3     E
#>  9     9      4     E
#>  9     9      5     E
#>  9     9      6     E
#>  9     9      7     E
#>  9     9      8     E
#>  9     9      9     E

# A utility whose scores psi00 and psi11 do not add up to 100
boin12_decision_table(target = 0.35, target_resp = 0.25, n = c(0, 3),
                      psi00 = 30, psi11 = 50)
#> BOIN12 desirability score table
#>   DLT limit / response limit : 0.35 / 0.25
#>   lambda_e / lambda_d        : 0.276 / 0.419
#>   utility scores             : 100 (response only), 30 (neither), 50 (both), 0 (DLT only)
#>   utility benchmark          : 67.62
#>   E: inadmissible (DLT cutoff 0.95, response cutoff 0.9)
#> 
#>  n n_tox n_resp n_tox_resp score
#>  0     0      0          0    12
#>  3     0      0          0     4
#>  3     0      1          0    10
#>  3     0      2          0    15
#>  3     0      3          0    17
#>  3     1      0          0     2
#>  3     1      1          0     7
#>  3     1      1          1     6
#>  3     1      2          0    13
#>  3     1      2          1    11
#>  3     1      3          1    16
#>  3     2      0          0     1
#>  3     2      1          0     5
#>  3     2      1          1     3
#>  3     2      2          1     9
#>  3     2      2          2     7
#>  3     2      3          2    13
#>  3     3      0          0     E
#>  3     3      1          1     E
#>  3     3      2          2     E
#>  3     3      3          3     E
```
