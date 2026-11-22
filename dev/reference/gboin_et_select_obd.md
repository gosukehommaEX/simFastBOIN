# Select the Optimal Biological Dose of gBOIN-ET Trials

Select the optimal biological dose (OBD) from the final data of one or
more trials of the gBOIN-ET design: the patients and the total
normalized toxicity and efficacy scores at every dose.

## Usage

``` r
gboin_et_select_obd(
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
  est_method = c("observed", "fp_linear", "fp_logistic"),
  obd_method = c("max_resp", "utility_weighted", "utility_truncated_linear",
    "utility_scoring"),
  min_obd_sample = 1,
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

  Numeric vector or matrix of the same shape as `n_pts`. Total of the
  normalized toxicity scores of the patients at each dose.

- n_resp:

  Numeric vector or matrix of the same shape as `n_pts`. Total of the
  normalized efficacy scores of the patients at each dose.

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

  Character scalar. Estimator of the mean efficacy scores: `"observed"`
  (the default) for the observed means, `"fp_linear"` for the linear
  fractional polynomial of Takeda et al. (2022), or `"fp_logistic"` for
  the fractional polynomial logistic regression of the 'boinet' package,
  which needs the 'mfp' package.

- obd_method:

  Character scalar. Selection rule: `"max_resp"` (the default) for the
  admissible dose with the highest estimated response rate up to the
  estimated MTD, or the dose that maximizes the utility
  `"utility_weighted"`, `"utility_truncated_linear"` or
  `"utility_scoring"` among the admissible doses.

- min_obd_sample:

  Integer scalar. Only doses given to at least this many patients can be
  selected. Defaults to 1. The simulations of Takeda et al. (2022) use
  6.

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

  Numeric matrix of the estimates of the mean efficacy scores, `NA` for
  the doses not tried.

## Details

The rules are those of
[`boin_et_select_obd`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/boin_et_select_obd.md)
with the mean scores in the place of the observed rates. A dose is
admissible when the posterior probability that its mean toxicity score
is below `p_eli` is at least `1 - cutoff_eli`, the posterior probability
that its mean efficacy score is above `p_eli_resp` is at least
`1 - cutoff_eli_resp`, and no lower dose is excluded by its toxicity
data. The mean toxicity scores of the doses tried are estimated by
isotonic regression with equal weights. Under `"max_resp"` the estimated
MTD is the highest dose whose estimate is as close to `target` as that
of the closest admissible dose, and the OBD is the admissible dose up to
the MTD with the highest estimated efficacy score; ties go to the lowest
dose. A dose with fewer than `min_obd_sample` patients takes part in
finding the MTD but cannot be selected; when no dose can be selected the
reason is `"no_admissible_dose"`.

The estimated efficacy scores (under `"max_resp"`) and the utilities are
compared after rounding to ten decimals, and ties go to the lowest dose.
Doses whose scores or utilities are equal as fractions then tie whatever
the rounding error of the platform; this happens regularly with
`"utility_scoring"`, whose default scores give different doses the same
utility.
[`boin_et_select_obd`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/boin_et_select_obd.md)
follows 'boinet' and compares the unrounded values, so that for whole
numbers of DLTs and responses the two functions select the same dose
under `"max_resp"` and can differ in such ties under a utility. The
estimates in `probe` are rounded to ten decimals.

`est_method = "fp_linear"` fits, by least squares to the scores of the
patients, the 36 models `b0 + b1 * dose^p1 + b2 * dose^p2` with the
powers -2, -1, -0.5, 0 (the logarithm), 0.5, 1, 2 and 3, where a
repeated power `p` stands for `dose^p` and `dose^p * log(dose)` and the
doses are the dose levels 1, 2 and so on, and takes the fitted values of
the model with the smallest residual sum of squares. This is the
estimator of Section 2.3 of Takeda et al. (2022), as their program
computes it (a normal model with identity link). With three doses or
fewer every model passes through the observed means, which are then
used. The fit runs in C++.

The program of Takeda et al. (2022) differs from this function in two
further points that `min_obd_sample = 6` and `boinet_compatible = TRUE`
bring close: it selects only among doses with at least six patients, and
it excludes no dose at the end of the trial.

## References

Takeda, K., Morita, S. and Taguri, M. (2022). gBOIN-ET: The generalized
Bayesian optimal interval design for optimal dose-finding accounting for
ordinal graded efficacy and toxicity in early clinical trials.
Biometrical Journal, 64(7), 1178-1191.
[doi:10.1002/bimj.202100263](https://doi.org/10.1002/bimj.202100263)

## See also

[`sim_gboin_et`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_gboin_et.md),
[`gboin_et_simulate`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/gboin_et_simulate.md),
[`boin_et_select_obd`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/boin_et_select_obd.md)

## Examples

``` r
# One trial: patients and total normalized scores at four doses
gboin_et_select_obd(
  n_pts = c(3, 6, 12, 9), n_tox = c(0, 1, 3.5, 4),
  n_resp = c(0.5, 2, 8, 6.5),
  target = 0.47 / 1.5, target_resp = 0.95 / 1.5
)
#> $obd
#> [1] 3
#> 
#> $reason
#> [1] "selected"
#> 
#> $probe
#>            DL1       DL2       DL3       DL4
#> [1,] 0.1666667 0.3333333 0.6666667 0.7222222
#> 

# The efficacy scores estimated by the linear fractional polynomial
gboin_et_select_obd(
  n_pts = c(3, 6, 12, 9), n_tox = c(0, 1, 3.5, 4),
  n_resp = c(0.5, 2, 8, 6.5),
  target = 0.47 / 1.5, target_resp = 0.95 / 1.5, est_method = "fp_linear"
)$probe
#>            DL1       DL2       DL3       DL4
#> [1,] 0.1290505 0.3691855 0.6555476 0.7256849
```
