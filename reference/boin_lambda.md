# Escalation and De-escalation Interval Boundaries

Compute the two boundaries that define the BOIN decision interval. A
dose is escalated when the observed DLT rate falls at or below
`lambda_e` and de-escalated when it reaches `lambda_d`.

## Usage

``` r
boin_lambda(target, p_saf = NULL, p_tox = NULL)
```

## Arguments

- target:

  Numeric scalar. Target DLT probability, for example 0.30.

- p_saf:

  Numeric scalar. Highest DLT probability deemed subtherapeutic, so that
  dose escalation should be undertaken. Defaults to `0.6 * target`.

- p_tox:

  Numeric scalar. Lowest DLT probability deemed overly toxic, so that
  dose de-escalation is required. Defaults to `1.4 * target`.

## Value

A list with components `lambda_e` and `lambda_d`.

## Details

The boundaries minimize the probability of incorrect dose assignment
under a three-point hypothesis on the DLT probability at the current
dose. Values of `p_saf` and `p_tox` close to `target` should be avoided,
because the sample sizes of phase I trials cannot distinguish the target
rate from rates close to it.

## References

Liu S. and Yuan, Y. (2015). Bayesian Optimal Interval Designs for Phase
I Clinical Trials. Journal of the Royal Statistical Society: Series C,
64, 507-523.

## Examples

``` r
boin_lambda(target = 0.30)
#> $lambda_e
#> [1] 0.2364907
#> 
#> $lambda_d
#> [1] 0.3585195
#> 

boin_lambda(target = 0.25, p_saf = 0.12, p_tox = 0.40)
#> $lambda_e
#> [1] 0.1788381
#> 
#> $lambda_d
#> [1] 0.3219281
#> 
```
