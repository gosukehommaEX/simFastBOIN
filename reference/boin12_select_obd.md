# Select the Optimal Biological Dose of BOIN12 Trials

Select the optimal biological dose (OBD) from the final data of one or
more trials of the BOIN12 design with the two-step procedure of Lin et
al. (2020): the MTD from the isotonic DLT rates, then the admissible
dose up to the MTD with the highest estimated utility.

## Usage

``` r
boin12_select_obd(
  n_pts,
  n_tox,
  n_resp,
  n_tox_resp = NULL,
  target,
  target_resp,
  cutoff_eli = 0.95,
  cutoff_eli_resp = 0.9,
  psi00 = 40,
  psi11 = 60,
  benchmark_weight = 0.5,
  utility_prior = c(1, 1),
  est_method = c("posterior_mean", "observed", "tail_probability"),
  early_stop = NULL
)
```

## Arguments

- n_pts:

  Integer vector of patients treated at each dose, or a matrix with one
  row per trial.

- n_tox, n_resp:

  DLT and response counts of the same shape as `n_pts`.

- n_tox_resp:

  Counts of the patients with both a DLT and a response, of the same
  shape as `n_pts`, or `NULL`. Needed unless `psi00 + psi11` is 100, in
  which case the utility depends on the margins alone. Defaults to
  `NULL`.

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
  [`sim_bf_boin_et`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_bf_boin_et.md),
  whereas the `target_resp` of
  [`sim_boin_et`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin_et.md)
  is a target.

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

- est_method:

  Character scalar. Estimate of the utility of a dose:
  `"posterior_mean"` (the default), the posterior mean of the
  quasi-binomial model given in Section 4 of the Supplement of Lin et
  al. (2020); `"observed"`, the observed mean utility, which the example
  of Section 6 of that Supplement reports; or `"tail_probability"`, the
  desirability used during the trial, which the 'escalation' package
  uses.

- early_stop:

  Logical vector or `NULL`. `TRUE` for the trials stopped early, which
  select no OBD. Defaults to `NULL`, no early stop.

## Value

A list with components

- obd:

  Integer vector of the selected dose levels, `NA` when no dose is
  selected.

- reason:

  Character vector: `"selected"`, `"early_stop"` or
  `"no_admissible_dose"`.

- utility:

  Numeric matrix of the utility estimates (on the scale 0 to 100; for
  `"tail_probability"`, 100 times the probability), `NA` for the doses
  not tried.

## Details

Only the doses with at least one patient are considered. The DLT rates
are made monotone by isotonic regression with equal weights, and the MTD
is the dose whose estimate is closest to `target`, ties to the higher
dose as in the 'escalation' package. The OBD is the admissible dose, in
the sense of
[`boin12_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/boin12_simulate.md),
at or below the MTD with the highest estimate, ties to the lower dose.
When no dose up to the MTD is admissible, no OBD is selected.

## References

Lin, R., Zhou, Y., Yan, F., Li, D. and Yuan, Y. (2020). BOIN12: Bayesian
optimal interval phase I/II trial design for utility-based dose finding
in immunotherapy and targeted therapies. JCO Precision Oncology, 4,
1393-1402.
[doi:10.1200/PO.20.00257](https://doi.org/10.1200/PO.20.00257)

## See also

[`sim_boin12`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin12.md),
[`boin12_simulate`](https://gosukehommaex.github.io/simFastBOIN/reference/boin12_simulate.md)

## Examples

``` r
# The example of Section 6 of the Supplement of Lin et al. (2020)
boin12_select_obd(
  n_pts = c(3, 3, 6, 12, 6), n_tox = c(0, 1, 1, 3, 3),
  n_resp = c(0, 1, 2, 6, 3), target = 0.35, target_resp = 0.25,
  est_method = "observed"
)
#> $obd
#> [1] 4
#> 
#> $reason
#> [1] "selected"
#> 
#> $utility
#>      DL1      DL2      DL3 DL4 DL5
#> [1,]  40 46.66667 53.33333  60  50
#> 
```
