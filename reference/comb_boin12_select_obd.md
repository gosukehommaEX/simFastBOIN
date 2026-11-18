# Select the OBD Combination of Comb-BOIN12 Trials

Select the optimal biological dose combination (OBDC) of one or more
trials of the Comb-BOIN12 design of Lu et al. (2025) from their final
data. The DLT rates are made monotone along both drugs by bivariate
isotonic regression, the maximum tolerated dose combination (MTDC)
bounds the candidates, and the OBDC is the admissible candidate with the
highest estimated utility.

## Usage

``` r
comb_boin12_select_obd(
  n_pts,
  n_tox,
  n_resp,
  target,
  target_resp,
  cutoff_eli = 0.95,
  cutoff_eli_resp = 0.9,
  extrasafe = FALSE,
  offset = 0.05,
  psi00 = 40,
  psi11 = 60,
  est_method = c("logistic", "observed"),
  obd_bound = c("rectangle", "contour"),
  dose_a = NULL,
  dose_b = NULL,
  logistic_scale = c(2.5, 2.5),
  n_draws = 4000,
  early_stop = NULL
)
```

## Arguments

- n_pts, n_tox, n_resp:

  Numeric matrices for a single trial, with one row per level of drug A
  and one column per level of drug B, or arrays with dimensions trial by
  level of drug A by level of drug B, of the patients treated, the DLTs
  and the responses.

- target, target_resp, cutoff_eli, cutoff_eli_resp, extrasafe, offset:

  As in
  [`comb_boin12_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/comb_boin12_simulate.md).
  The admissibility of the combinations is judged again from the final
  data with the same rules.

- psi00, psi11:

  Numeric scalars from 0 to 100. Utilities of a patient with neither a
  DLT nor a response and with both, as in
  [`sim_boin12`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin12.md).
  Default to 40 and 60.

- est_method:

  Character scalar. Estimate of the response probability: `"logistic"`
  (the default) the posterior mean under the logistic model of Lu et al.
  (2025), see Details; `"observed"` the observed response rate.

- obd_bound:

  Character scalar. Combinations from which the OBDC is chosen:
  `"rectangle"` (the default) those with both levels at most those of
  the MTDC, as Lu et al. (2025) write; `"contour"` those, in every
  column (level of drug B), at most the combination of that column whose
  isotonic DLT rate is closest to `target`, the rule that reproduces the
  results published by Lu et al. (2025). See Details.

- dose_a, dose_b:

  Numeric vectors or `NULL`. Dose values of the levels of drug A and of
  drug B in the logistic model. Default to `NULL`, the level numbers.
  Each is standardized to mean 0 and standard deviation 0.5.

- logistic_scale:

  Numeric vector of two positive numbers. Scales of the Cauchy priors of
  the intercept and of the other coefficients of the logistic model.
  Defaults to `c(2.5, 2.5)`, as in Lu et al. (2025).

- n_draws:

  Integer scalar. Number of points of the importance sampling of the
  logistic model. Defaults to 4000.

- early_stop:

  Logical vector or `NULL`. Trials stopped early, which select no OBDC.
  Defaults to `NULL`, no trial.

## Value

A list with components

- obd_a, obd_b:

  Integer vectors, the levels of drug A and drug B of the OBDC of every
  trial (`NA` when none is selected).

- mtd_a, mtd_b:

  Integer vectors, the levels of the MTDC.

- reason:

  Character vector: `"selected"`, `"early_stop"` or
  `"no_admissible_dose"`.

- tox_est:

  Array of the isotonic estimates of the DLT rates.

- resp_est:

  Array of the estimates of the response probabilities of the tried
  combinations (`NA` elsewhere and for the trials with no selection).

- utility:

  Array of the estimated utilities of the tried combinations, on the
  scale 0 to 100.

## Details

The DLT rates are estimated by `(y + 0.05) / (n + 0.1)` with weights
`n + 0.1`, as in
[`comb_boin_select_mtd`](https://gosukehommaex.github.io/simFastBOIN/reference/comb_boin_select_mtd.md),
and made non-decreasing in both levels by the weighted least-squares fit
of bivariate isotonic regression, computed by the algorithm of Dykstra
and Robertson (1982) until no estimate changes by `1e-13` in a cycle.
The MTDC is the tried combination whose estimate, plus `1e-5` times the
sum of its two levels, is closest to `target`: below the target the
small amount favors the higher combination, above it the lower one.
Among combinations whose distances differ by at most `1e-8`, such as two
combinations with the same sum of levels in one block of the isotonic
regression, the first in column-major order is the MTDC. Untried
combinations are never selected.

The utility of a combination is estimated by
`100 * pE * (1 - pT) + psi00 * (1 - pT) * (1 - pE) + psi11 * pT * pE`,
with `pT` the isotonic estimate and `pE` the estimate of `est_method`.
When `psi00 + psi11` is 100 it is a linear function of `pE - wT * pT`,
the estimated utility of Lu et al. (2025), with `wT = psi00 / psi11`.
The OBDC is the admissible tried candidate with the highest estimate;
among the candidates within `1e-6` of it, the one with the lowest sum of
the levels and then the first in column-major order.

Lu et al. (2025) state the candidates as the combinations (j, k) with j
and k at most the levels of the MTDC. Applied to the trials that
reproduce their allocation percentages, that rule gives selection
percentages far from those published; `obd_bound = "contour"`, which
bounds every level of drug B by its own MTD along drug A, agrees with
their Figure 2 and Supplementary Figure S2 (see
`inst/reproduce/reproduce_Lu2025.R`). It treats the two drugs
differently.

The logistic model is `logit(pE) = g0 + g1 a + g2 b + g3 a^2 + g4 b^2`,
with `a` and `b` the standardized doses (mean 0 and standard deviation
0.5 over the levels of each drug, the scaling of Gelman et al., 2008,
for these priors; Lu et al., 2025, do not state one) and independent
Cauchy priors centered at 0. Lu et al. draw the posterior with Gibbs
sampling; here its mean is computed by importance sampling with
`n_draws` Halton points: the proposal is a multivariate t distribution
with 4 degrees of freedom centered at the posterior mode, with the
inverse of the negative Hessian as scale matrix, and it is fitted twice
more to the weighted mean and covariance of the weighted points. No
random number is used, so that the selection does not change the random
number stream. In checks against long Markov chain Monte Carlo runs the
estimates were within 0.01.

## References

Lu, M., Zhang, J., Yuan, Y. and Lin, R. (2025). Comb-BOIN12: A
utility-based Bayesian optimal interval design for dose optimization in
cancer drug-combination trials. Statistics in Biopharmaceutical
Research, 17(2), 266-276.
[doi:10.1080/19466315.2024.2370403](https://doi.org/10.1080/19466315.2024.2370403)

Gelman, A., Jakulin, A., Pittau, M. G. and Su, Y.-S. (2008). A weakly
informative default prior distribution for logistic and other regression
models. The Annals of Applied Statistics, 2(4), 1360-1383.
[doi:10.1214/08-AOAS191](https://doi.org/10.1214/08-AOAS191)

Dykstra, R. L. and Robertson, T. (1982). An algorithm for isotonic
regression for two or more independent variables. The Annals of
Statistics, 10(3), 708-716.
[doi:10.1214/aos/1176345866](https://doi.org/10.1214/aos/1176345866)

## See also

[`sim_comb_boin12`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_comb_boin12.md),
[`comb_boin12_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/comb_boin12_simulate.md)

## Examples

``` r
n_pts <- matrix(c(3, 6, 0, 9, 12, 3), nrow = 2)
n_tox <- matrix(c(0, 1, 0, 1, 2, 1), nrow = 2)
n_resp <- matrix(c(0, 1, 0, 3, 6, 2), nrow = 2)
comb_boin12_select_obd(n_pts, n_tox, n_resp, target = 0.33,
                       target_resp = 0.2)
#> $obd_a
#> [1] 1
#> 
#> $obd_b
#> [1] 3
#> 
#> $mtd_a
#> [1] 2
#> 
#> $mtd_b
#> [1] 3
#> 
#> $reason
#> [1] "selected"
#> 
#> $tox_est
#> , , 1
#> 
#>            [,1]      [,2]
#> [1,] 0.01612903 0.1405229
#> 
#> , , 2
#> 
#>           [,1]      [,2]
#> [1,] 0.1405229 0.1405229
#> 
#> , , 3
#> 
#>           [,1]      [,2]
#> [1,] 0.1694215 0.3387097
#> 
#> 
#> $resp_est
#> , , 1
#> 
#>           [,1]      [,2]
#> [1,] 0.1298339 0.1631064
#> 
#> , , 2
#> 
#>      [,1]      [,2]
#> [1,]   NA 0.3492024
#> 
#> , , 3
#> 
#>          [,1]      [,2]
#> [1,] 0.490028 0.5914899
#> 
#> 
#> $utility
#> , , 1
#> 
#>          [,1]     [,2]
#> [1,] 47.14487 44.16547
#> 
#> , , 2
#> 
#>      [,1]     [,2]
#> [1,]   NA 55.33123
#> 
#> , , 3
#> 
#>          [,1]     [,2]
#> [1,] 62.62482 61.94101
#> 
#> 
```
