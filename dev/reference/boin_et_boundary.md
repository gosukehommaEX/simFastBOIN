# Boundaries of the BOIN-ET Design

Find the lower and upper DLT boundaries and the lower response boundary
of the BOIN-ET design of Takeda, Taguri and Morita (2018). The three
boundaries minimize the probability of an incorrect dose assignment
under six hypotheses about the DLT and response probabilities at the
current dose, and are searched on a grid of step 0.01, as in the
'boinet' package.

## Usage

``` r
boin_et_boundary(
  target,
  target_resp,
  p_saf = NULL,
  p_tox = NULL,
  p_saf_resp = NULL,
  prior = rep(1/6, 6)
)
```

## Arguments

- target:

  Numeric scalar. Target DLT probability, for example 0.30.

- target_resp:

  Numeric scalar. Target response probability, for example 0.60.

- p_saf:

  Numeric scalar. Highest DLT probability deemed subtherapeutic.
  Defaults to `0.1 * target`, the value recommended by Takeda et al.
  (2018), which is far below the `0.6 * target` of the BOIN design.

- p_tox:

  Numeric scalar. Lowest DLT probability deemed overly toxic. Defaults
  to `1.4 * target`.

- p_saf_resp:

  Numeric scalar. Highest response probability deemed subtherapeutic.
  Defaults to `0.6 * target_resp`.

- prior:

  Numeric vector of length 6. Prior probabilities of the six hypotheses.
  Defaults to `rep(1 / 6, 6)`, the noninformative choice of the article.
  The third hypothesis, under which every decision is correct, does not
  enter the criterion.

## Value

An object of class `boin_et_boundary`, which is a list with components

- lambda1:

  Lower DLT boundary. The dose is escalated when the observed DLT rate
  is at most `lambda1` and the observed response rate is at most `eta1`.

- lambda2:

  Upper DLT boundary. The dose is de-escalated when the observed DLT
  rate is at least `lambda2`.

- eta1:

  Lower response boundary.

- prob_incorrect:

  The minimized probability of an incorrect decision.

- target, target_resp, p_saf, p_tox, p_saf_resp, prior:

  The design parameters.

## Details

The grids run from `p_saf` to `target` for `lambda1`, from `target` to
`p_tox` for `lambda2` and from `p_saf_resp` to `target_resp` for `eta1`,
all by 0.01. At each grid point the probability of an incorrect decision
is evaluated with binomial probabilities for 100 patients, and the first
minimum in the order of the grids is returned. This is the search of
[`boinet::gridoptim()`](https://rdrr.io/pkg/boinet/man/gridoptim.html),
including its expression for the fourth hypothesis, and with the default
design parameters and a target of 0.30 for the DLT and 0.60 for the
response it gives the boundaries 0.14, 0.35 and 0.48 reported by Takeda
et al. (2018).

## References

Takeda, K., Taguri, M. and Morita, S. (2018). BOIN-ET: Bayesian optimal
interval design for dose finding based on both efficacy and toxicity
outcomes. Pharmaceutical Statistics, 17(4), 383-395.
[doi:10.1002/pst.1864](https://doi.org/10.1002/pst.1864)

Yamaguchi, Y., Takeda, K. and Maruo, K. (2025). boinet: An R Package for
Bayesian Optimal Interval Design for Dose Finding Based on Both Efficacy
and Toxicity Outcomes. The R Journal, 17(2), 15-34.
[doi:10.32614/RJ-2025-012](https://doi.org/10.32614/RJ-2025-012)

## See also

[`boin_et_decision_table`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/boin_et_decision_table.md),
[`sim_boin_et`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_boin_et.md)

## Examples

``` r
bd <- boin_et_boundary(target = 0.30, target_resp = 0.60)
bd
#> BOIN-ET boundaries
#>   target DLT rate      : 0.3
#>   target response rate : 0.6
#>   p_saf / p_tox        : 0.03 / 0.42
#>   p_saf_resp           : 0.36
#> 
#>   lambda1 (escalate)   : 0.14
#>   lambda2 (de-escalate): 0.35
#>   eta1 (response)      : 0.48
#>   probability of an incorrect decision: 0.04703

c(bd$lambda1, bd$lambda2, bd$eta1)
#> [1] 0.14 0.35 0.48
```
