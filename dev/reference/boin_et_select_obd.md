# Select the Optimal Biological Dose of BOIN-ET Trials

Select the optimal biological dose (OBD) from the final data of one or
more trials of the BOIN-ET design, with the estimators and selection
rules of the 'boinet' package.

## Usage

``` r
boin_et_select_obd(
  n_pts,
  n_tox,
  n_resp,
  target,
  target_resp,
  p_saf = NULL,
  p_tox = NULL,
  p_saf_resp = NULL,
  cutoff_eli = 0.95,
  cutoff_eli_resp = 0.99,
  p_eli = NULL,
  p_eli_resp = NULL,
  est_method = c("observed", "multi_iso", "fp_logistic"),
  obd_method = c("max_resp", "utility_weighted", "utility_truncated_linear",
    "utility_scoring"),
  w1 = 0.33,
  w2 = 1.09,
  t_low = NULL,
  t_upp = NULL,
  e_low = NULL,
  e_upp = NULL,
  psi00 = 40,
  psi11 = 60,
  early_stop = NULL,
  boinet_compatible = FALSE
)
```

## Arguments

- n_pts:

  Integer vector of patients treated at each dose, or a matrix with one
  row per trial.

- n_tox:

  DLT counts of the same shape as `n_pts`.

- n_resp:

  Response counts of the same shape as `n_pts`.

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

- cutoff_eli:

  Numeric scalar. A dose is no longer admissible when the posterior
  probability that its DLT rate is below `p_eli` falls below
  `1 - cutoff_eli`. Defaults to 0.95.

- cutoff_eli_resp:

  Numeric scalar. A dose is no longer admissible when the posterior
  probability that its response rate is above `p_eli_resp` falls below
  `1 - cutoff_eli_resp`. Defaults to 0.99.

- p_eli:

  Numeric scalar or `NULL`. DLT probability of the rule above. Defaults
  to `NULL`, which uses `target` as the 'boinet' package does. Takeda et
  al. (2018) used 0.40 with a cutoff of 0.90.

- p_eli_resp:

  Numeric scalar or `NULL`. Response probability of the rule above.
  Defaults to `NULL`, which uses `p_saf_resp` as the 'boinet' package
  does. Takeda et al. (2018) used 0.20 with a cutoff of 0.90.

- est_method:

  Character scalar. Estimator of the response rates: `"observed"` (the
  default) for the observed rates, `"multi_iso"` for the model-averaged
  unimodal isotonic regression of Lin and Yin (2017), or `"fp_logistic"`
  for the fractional polynomial logistic regression of Takeda et al.
  (2018), which needs the 'mfp' package.

- obd_method:

  Character scalar. Selection rule: `"max_resp"` (the default) for the
  admissible dose with the highest estimated response rate up to the
  estimated MTD, or the dose that maximizes the utility
  `"utility_weighted"`, `"utility_truncated_linear"` or
  `"utility_scoring"` among the admissible doses.

- w1, w2:

  Numeric scalars. Weights of the weighted utility
  `q - w1 * p - w2 * p * (p > p_tox)`, where `p` and `q` are the
  estimated DLT and response rates. Default to 0.33 and 1.09.

- t_low, t_upp:

  Numeric scalars or `NULL`. DLT rates between which the truncated
  linear utility falls from one to zero. Default to `p_saf` and `p_tox`.

- e_low, e_upp:

  Numeric scalars or `NULL`. Response rates between which the truncated
  linear utility rises from zero to one. Default to `p_saf_resp / 2` and
  `target_resp`.

- psi00, psi11:

  Numeric scalars. Scores of the scoring utility for a patient with
  neither DLT nor response and for a patient with both, the scores of a
  response without DLT and of a DLT without response being 100 and 0.
  Default to 40 and 60.

- early_stop:

  Logical vector or `NULL`. `TRUE` for the trials stopped early, which
  select no OBD. Defaults to `NULL`, no early stop.

- boinet_compatible:

  Logical scalar. Apply the rules of the 'boinet' package (version
  1.6.0) exactly, in which a dose excluded by its DLT data does not
  exclude the doses above it. Defaults to `FALSE`; see Details.

## Value

A list with components

- obd:

  Integer vector of the selected dose levels, `NA` when no dose is
  selected.

- reason:

  Character vector: `"selected"`, `"early_stop"` or
  `"no_admissible_dose"`.

- probe:

  Numeric matrix of the response estimates, `NA` for the doses not
  tried.

## Details

Only the doses with at least one patient are considered. A dose is
admissible when the posterior probability that its DLT rate is below
`p_eli` is at least `1 - cutoff_eli` and the posterior probability that
its response rate is above `p_eli_resp` is at least
`1 - cutoff_eli_resp`, under uniform priors, and when no lower dose is
excluded by its DLT data. When no dose is admissible, no OBD is
selected. Otherwise the DLT rates of the doses tried are estimated by
isotonic regression of the observed rates with equal weights. Under
`"max_resp"` the estimated MTD is the highest dose whose estimate is as
close to `target` as that of the closest admissible dose, and the OBD is
the admissible dose up to the MTD with the highest response estimate.
Ties go to the lowest dose. These are the rules of
[`boinet::obd.select()`](https://rdrr.io/pkg/boinet/man/obd.select.html),
except that in 'boinet' a dose excluded by its DLT data does not exclude
the doses above it, which can then be selected;
`boinet_compatible = TRUE` applies this rule, and
[`boin_et_simulate`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/boin_et_simulate.md)
explains the change.

The model-based estimates are computed in R once for every distinct set
of final data; the selection itself runs in C++.

## References

Takeda, K., Taguri, M. and Morita, S. (2018). BOIN-ET: Bayesian optimal
interval design for dose finding based on both efficacy and toxicity
outcomes. Pharmaceutical Statistics, 17(4), 383-395.
[doi:10.1002/pst.1864](https://doi.org/10.1002/pst.1864)

Yamaguchi, Y., Takeda, K., Yoshida, S. and Maruo, K. (2024). Optimal
biological dose selection in dose-finding trials with model-assisted
designs based on efficacy and toxicity: a simulation study. Journal of
Biopharmaceutical Statistics, 34(3), 379-393.
[doi:10.1080/10543406.2023.2202259](https://doi.org/10.1080/10543406.2023.2202259)

Lin, R. and Yin, G. (2017). STEIN: A simple toxicity and efficacy
interval design for seamless phase I/II clinical trials. Statistics in
Medicine, 36(26), 4106-4120.
[doi:10.1002/sim.7428](https://doi.org/10.1002/sim.7428)

## See also

[`sim_boin_et`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_boin_et.md),
[`boin_et_simulate`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/boin_et_simulate.md)

## Examples

``` r
# One trial: patients, DLTs and responses at four doses
boin_et_select_obd(
  n_pts = c(3, 6, 12, 3), n_tox = c(0, 1, 3, 2), n_resp = c(0, 2, 7, 2),
  target = 0.30, target_resp = 0.60
)
#> $obd
#> [1] 3
#> 
#> $reason
#> [1] "selected"
#> 
#> $probe
#>      DL1       DL2       DL3       DL4
#> [1,]   0 0.3333333 0.5833333 0.6666667
#> 

# The OBD under a utility
boin_et_select_obd(
  n_pts = c(3, 6, 12, 3), n_tox = c(0, 1, 3, 2), n_resp = c(0, 2, 7, 2),
  target = 0.30, target_resp = 0.60, obd_method = "utility_scoring"
)$obd
#> [1] 3
```
