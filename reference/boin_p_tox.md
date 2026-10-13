# Toxic Threshold Giving a Required De-escalation Boundary

Find the value of `p_tox` for which the BOIN de-escalation boundary
`lambda_d` equals a required value. This inverts
[`boin_lambda`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_lambda.md),
which maps `p_tox` to `lambda_d`.

## Usage

``` r
boin_p_tox(target, lambda_d)
```

## Arguments

- target:

  Numeric scalar. Target DLT probability.

- lambda_d:

  Numeric scalar. Required de-escalation boundary. Must lie strictly
  between `target` and 1.

## Value

A numeric scalar, the value of `p_tox` that produces `lambda_d`.

## Details

`lambda_d` increases strictly with `p_tox`, from `target` in the limit
as `p_tox` approaches `target` to 1 in the limit as `p_tox` approaches
1, so the solution exists and is unique for any `lambda_d` strictly
between those bounds. The root is found numerically.

A boundary close to the target needs a toxic threshold close to the
target, and the rest of the package refuses a `p_tox` within ten percent
of `target`, on the grounds that phase I sample sizes cannot tell the
two rates apart. The value is still returned, since the inversion is
well defined, but a warning names the smallest boundary that leaves a
usable threshold. For a target of 0.30 that smallest boundary is about
0.315.

A design team that wants to tighten the de-escalation boundary usually
states the boundary rather than the threshold behind it. This function
turns that statement into the `p_tox` to pass to the other functions, so
that the whole design follows from it.

## See also

[`boin_lambda`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_lambda.md),
[`boin_boundary`](https://gosukehommaex.github.io/simFastBOIN/reference/boin_boundary.md)

## Examples

``` r
# A target of 0.30 gives a de-escalation boundary of 0.359 by default
boin_lambda(target = 0.30)$lambda_d
#> [1] 0.3585195

# The toxic threshold that tightens it to 0.33
p_tox <- boin_p_tox(target = 0.30, lambda_d = 0.33)
p_tox
#> [1] 0.3609529

# which is what it claims to be
boin_lambda(target = 0.30, p_tox = p_tox)$lambda_d
#> [1] 0.33

# A boundary too close to the target leaves an unusable threshold, and says so
boin_p_tox(target = 0.30, lambda_d = 0.31)
#> Warning: a de-escalation boundary of 0.31 needs 'p_tox' = 0.32012, which is too close to 'target' for the other functions to accept. The smallest usable boundary for a target of 0.3 is about 0.3149.
#> [1] 0.3201199

# and shifts the de-escalation boundary down by one DLT at every cohort end
as.data.frame(boin_boundary(0.30, 18))[seq(3, 18, 3), ]
#>    n_pts escalate_if_dlt_leq deescalate_if_dlt_geq eliminate_if_dlt_geq
#> 3      3                   0                     2                    3
#> 6      6                   1                     3                    4
#> 9      9                   2                     4                    5
#> 12    12                   2                     5                    7
#> 15    15                   3                     6                    8
#> 18    18                   4                     7                    9
as.data.frame(boin_boundary(0.30, 18, p_tox = p_tox))[seq(3, 18, 3), ]
#>    n_pts escalate_if_dlt_leq deescalate_if_dlt_geq eliminate_if_dlt_geq
#> 3      3                   0                     1                    3
#> 6      6                   1                     2                    4
#> 9      9                   2                     3                    5
#> 12    12                   2                     4                    7
#> 15    15                   3                     5                    8
#> 18    18                   4                     6                    9
```
