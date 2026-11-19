# Select the Optimal Dose Combination of BOIN-ETC Trials

Select the optimal dose combination (ODC) of one or more trials of the
BOIN-ETC design of Kakizume et al. (2024) from their final data. The DLT
rates are made monotone along both drugs by bivariate isotonic
regression, the maximum tolerated dose combination (MTDC) of every row
bounds the candidates of that row, and the ODC is the candidate with the
highest estimated utility.

## Usage

``` r
boin_etc_select_obd(
  n_pts,
  n_tox,
  n_resp,
  n_tox_resp = NULL,
  target,
  psi00 = 40,
  psi11 = 60,
  utility_prior = c(1, 1),
  odc_tie = c("lowest", "all"),
  cutoff_eli = 0.95,
  p_eli = NULL,
  excluded = NULL,
  early_stop = NULL
)
```

## Arguments

- n_pts, n_tox, n_resp:

  Numeric matrices for a single trial, with one row per level of drug A
  and one column per level of drug B, or arrays with dimensions trial by
  level of drug A by level of drug B, of the patients treated, the DLTs
  and the responses.

- n_tox_resp:

  Array or matrix of the same size, or `NULL`. Patients with both a DLT
  and a response. Needed only when `psi00 + psi11` is not 100. Defaults
  to `NULL`.

- target:

  Numeric scalar. Target DLT probability, phi.

- psi00, psi11:

  Numeric scalars from 0 to 100. Utilities of a patient with neither a
  DLT nor a response and with both; a response without DLT scores 100
  and a DLT without response 0. Default to 40 and 60, the utilities
  (100, 40, 60, 0) of Kakizume et al. (2024).

- utility_prior:

  Numeric vector of two positive numbers. Beta prior of the utility of a
  combination divided by 100. Defaults to `c(1, 1)`.

- odc_tie:

  Character scalar. How tied candidates are reported: `"lowest"` (the
  default) selects the one with the lowest sum of the two levels;
  `"all"` also marks every tied candidate in `tied`, which
  [`sim_boin_etc`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin_etc.md)
  then counts as selected, as Kakizume et al. (2024) do. The ODC
  returned in `obd_a` and `obd_b` is the same in both cases.

- cutoff_eli, p_eli:

  As in
  [`boin_etc_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_etc_simulate.md).
  Used only when `excluded` is `NULL`.

- excluded:

  Logical array or matrix of the same size, or `NULL`. Combinations
  excluded by the DLT rule during the trial (component `excluded` of
  [`boin_etc_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_etc_simulate.md)),
  which are never selected. Defaults to `NULL`: a combination is
  excluded when its own final DLT data fail the rule of `cutoff_eli`.

- early_stop:

  Logical vector or `NULL`. Trials stopped early, which select no ODC.
  Defaults to `NULL`, no trial.

## Value

A list with components

- obd_a, obd_b:

  Integer vectors, the levels of drug A and drug B of the ODC of every
  trial (`NA` when none is selected).

- n_tied:

  Integer vector, the number of tied candidates.

- tied:

  Logical array, the tied candidates of every trial (with
  `odc_tie = "lowest"`, only the selected ODC).

- mtd_col:

  Integer matrix, trial by level of drug A, the level of drug B of the
  MTDC of every row.

- reason:

  Character vector: `"selected"`, `"early_stop"` or
  `"no_admissible_dose"`.

- tox_est:

  Array of the isotonic estimates of the DLT rates.

- utility:

  Array of the estimated utilities of the tried combinations, on the
  scale 0 to 100.

## Details

The DLT rate of every combination, tried or not, is estimated by
`(y + 0.05) / (n + 0.1)` (the beta(0.05, 0.05) posterior mean of
Kakizume et al., 2024) and made non-decreasing in both levels by the
unweighted least-squares fit of bivariate isotonic regression, as in the
code of the authors
([`Iso::biviso()`](https://rdrr.io/pkg/Iso/man/biviso.html) without
weights), computed by the algorithm of Dykstra and Robertson (1982)
until no estimate changes by `1e-13` in a cycle. In every row (level of
drug A) the MTDC is the combination whose estimate is closest to
`target`; among combinations within `1e-8` of the smallest distance,
which are tied in exact arithmetic, the highest level of drug B. The
candidates are the tried combinations that are not excluded and not
above the MTDC of their row.

The utility of a candidate is estimated by
`100 * (q + a) / (n + a + b)`, with `(a, b)` the `utility_prior` and the
quasi-binomial number of events
`q = (100 x + psi00 (n - y - x) + (psi00 + psi11 - 100) z) / 100` of `n`
patients with `y` DLTs, `x` responses and `z` both, the estimate of the
code of the authors (Kakizume et al., 2024, do not state it). The ODC is
the candidate with the highest estimate; the candidates within `1e-6` of
it are tied, and the one with the lowest sum of the levels, then the
first in column-major order, is selected.

The code of the authors also lets untried combinations (utility 50) and
excluded ones be selected, and counts every tied combination; on the
published trials, leaving out the untried ones changes the selection of
0.0 to 0.3 0.14 errors of about `1e-6`, so that which of tied MTDCs it
takes depends on them; the published percentages of trials selecting a
correct ODC are on average about 1 point lower than with exact ties (see
`inst/reproduce/reproduce_Kakizume2024.R`).

## References

Kakizume, T., Takeda, K., Taguri, M. and Morita, S. (2024). BOIN-ETC: A
Bayesian optimal interval design considering efficacy and toxicity to
identify the optimal dose combinations. Statistical Methods in Medical
Research, 33(4), 716-727.
[doi:10.1177/09622802241236936](https://doi.org/10.1177/09622802241236936)

Dykstra, R. L. and Robertson, T. (1982). An algorithm for isotonic
regression for two or more independent variables. The Annals of
Statistics, 10(3), 708-716.
[doi:10.1214/aos/1176345866](https://doi.org/10.1214/aos/1176345866)

## See also

[`sim_boin_etc`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin_etc.md),
[`boin_etc_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_etc_simulate.md)

## Examples

``` r
n_pts <- matrix(c(3, 6, 0, 9, 6, 3), nrow = 2)
n_tox <- matrix(c(0, 1, 0, 2, 2, 2), nrow = 2)
n_resp <- matrix(c(0, 2, 0, 5, 4, 2), nrow = 2)
boin_etc_select_obd(n_pts, n_tox, n_resp, target = 0.30)
#> $obd_a
#> [1] 1
#> 
#> $obd_b
#> [1] 3
#> 
#> $n_tied
#> [1] 1
#> 
#> $tied
#> , , B1
#> 
#>         A1    A2
#> [1,] FALSE FALSE
#> 
#> , , B2
#> 
#>         A1    A2
#> [1,] FALSE FALSE
#> 
#> , , B3
#> 
#>        A1    A2
#> [1,] TRUE FALSE
#> 
#> 
#> $mtd_col
#>      A1 A2
#> [1,]  3  2
#> 
#> $reason
#> [1] "selected"
#> 
#> $tox_est
#> , , B1
#> 
#>              A1        A2
#> [1,] 0.01612903 0.1721311
#> 
#> , , B2
#> 
#>             A1        A2
#> [1,] 0.3537801 0.3537801
#> 
#> , , B3
#> 
#>             A1        A2
#> [1,] 0.3537801 0.6612903
#> 
#> 
#> $utility
#> , , B1
#> 
#>      A1   A2
#> [1,] 44 52.5
#> 
#> , , B2
#> 
#>      A1       A2
#> [1,] NA 61.81818
#> 
#> , , B3
#> 
#>        A1 A2
#> [1,] 62.5 52
#> 
#> 
```
