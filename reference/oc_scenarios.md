# True Dose-Toxicity and Dose-Response Scenarios

Tabulate a list of scenarios, as given to the `_multi` simulation
functions, with the true utility of every dose and the doses that are
the true maximum tolerated dose (MTD) and the true optimal biological
dose (OBD). The result is the truth against which
[`oc_compare`](https://gosukehommaex.github.io/simFastBOIN/reference/oc_compare.md)
judges the designs, and it has
[`print()`](https://rdrr.io/r/base/print.html) and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods.

## Usage

``` r
oc_scenarios(
  scenarios,
  target,
  target_resp = NULL,
  tox_limit = target,
  psi00 = 40,
  psi11 = 60,
  resp_cor = 0,
  resp_cor_type = c("copula", "binary"),
  true_mtd = NULL,
  true_obd = NULL
)
```

## Arguments

- scenarios:

  List of scenarios. For a single agent, each element is a numeric
  vector of DLT probabilities or a list with the component `p_true`, and
  optionally `p_resp` (response probabilities) and `name`. For two drugs
  given together, `p_true` (and `p_resp`) are matrices with one row per
  level of drug A and one column per level of drug B. Either every
  scenario or none has `p_resp`. Names of the list are used when `name`
  is absent, and "Scenario 1", "Scenario 2", ... otherwise.

- target:

  Numeric scalar. Target DLT probability; the true MTD is the dose whose
  DLT probability is closest to it.

- target_resp:

  Numeric scalar or `NULL`. Lowest acceptable response probability, as
  `target_resp` of
  [`sim_boin12`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin12.md);
  a dose with a lower response probability cannot be the true OBD.
  Defaults to `NULL`, no such limit. The `target_resp` of
  [`sim_boin_et`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin_et.md),
  the target response probability of BOIN-ET, is a different quantity.

- tox_limit:

  Numeric scalar. Highest acceptable DLT probability; a dose with a
  higher DLT probability cannot be the true OBD. Defaults to `target`.

- psi00, psi11:

  Numeric scalars. Utilities of a patient with neither a DLT nor a
  response and with both, on the scale where a response without a DLT
  scores 100 and a DLT without a response scores 0. Default to 40 and
  60, as in
  [`sim_boin12`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin12.md).

- resp_cor:

  Numeric scalar. Correlation between the DLT and the response, used for
  the probabilities of the four outcomes. Defaults to 0, independence.

- resp_cor_type:

  Character scalar. How `resp_cor` is read, as in
  [`sim_boin12`](https://gosukehommaex.github.io/simFastBOIN/reference/sim_boin12.md):
  `"copula"` (the default), the correlation of the latent normal
  variables, or `"binary"`, the correlation of the two binary outcomes.

- true_mtd, true_obd:

  `NULL` (the default) or a list with one element per scenario, in the
  order of `scenarios` or named by the scenarios, giving the true MTDs
  or OBDs directly when the definition of an article differs from the
  one used here: a vector of dose levels for a single agent, a
  two-column matrix of the levels of drug A and drug B for a
  combination, or `NULL` for the computed doses. An empty vector means
  no such dose.

## Value

A data frame of class `c("oc_scenarios", "data.frame")` with one row per
scenario and dose (or combination) and the columns

- scenario:

  Factor of the scenario names, in the order given.

- dose:

  Dose level, or `dose_a` and `dose_b` for a combination.

- p_true:

  True DLT probability.

- p_resp, utility:

  True response probability and utility, when the scenarios have
  response probabilities.

- true_mtd:

  Logical, the dose is a true MTD.

- true_obd:

  Logical, the dose is a true OBD, when the scenarios have response
  probabilities.

The arguments are kept as attributes, with `type` (`"single"` or
`"combination"`) and `dims` (the number of doses, or the numbers of
levels of the two drugs).

## Details

The utility of a dose is
`100 * Pr(response, no DLT) + psi00 * Pr(neither) + psi11 * Pr(both)`.
The true MTDs are all the doses whose DLT probability is closest to
`target` (several when tied). The true OBDs are all the doses with the
highest utility among the doses whose DLT probability is at most
`tox_limit` and, when `target_resp` is given, whose response probability
is at least `target_resp`; a scenario may have none. Ties are decided
with a tolerance of `1e-10`.

Some articles use other definitions, for example Lin et al. (2020), who
count every dose within 3 points of the highest utility as an OBD; give
such doses through `true_obd`.

## References

Lin, R., Zhou, Y., Yan, F., Li, D. and Yuan, Y. (2020). BOIN12: Bayesian
optimal interval phase I/II trial design for utility-based dose finding
in immunotherapy and targeted therapies. JCO Precision Oncology, 4,
1393-1402.
[doi:10.1200/PO.20.00257](https://doi.org/10.1200/PO.20.00257)

## See also

[`oc_results`](https://gosukehommaex.github.io/simFastBOIN/reference/oc_results.md),
[`oc_compare`](https://gosukehommaex.github.io/simFastBOIN/reference/oc_compare.md),
[`print.oc_scenarios`](https://gosukehommaex.github.io/simFastBOIN/reference/print.oc_scenarios.md),
[`plot.oc_scenarios`](https://gosukehommaex.github.io/simFastBOIN/reference/plot.oc_scenarios.md)

## Examples

``` r
truth <- oc_scenarios(
  list(
    "Scenario 1" = list(p_true = c(0.05, 0.12, 0.27, 0.35, 0.50),
                        p_resp = c(0.20, 0.35, 0.36, 0.37, 0.38)),
    "Scenario 2" = list(p_true = c(0.05, 0.12, 0.27, 0.35, 0.50),
                        p_resp = c(0.01, 0.05, 0.30, 0.60, 0.60))
  ),
  target = 0.35, target_resp = 0.25
)
truth
#> True scenarios: 2 scenario(s), 5 doses 
#> True MTD: DLT probability closest to 0.35 
#> True OBD: highest utility among the doses with DLT probability <= 0.35 and response probability >= 0.25 
#> Utilities: 100 (response, no DLT), 40 (neither), 60 (both), 0 (DLT, no response)
#> 
#> Scenario 1
#>             1    2    3    4    5
#> p_true   0.05 0.12 0.27 0.35 0.50
#> p_resp   0.20 0.35 0.36 0.37 0.38
#> utility  50.0 56.2 50.8 48.2 42.8
#> true MTD                   *     
#> true OBD         *               
#> 
#> Scenario 2
#>             1    2    3    4    5
#> p_true   0.05 0.12 0.27 0.35 0.50
#> p_resp   0.01 0.05 0.30 0.60 0.60
#> utility  38.6 38.2 47.2 62.0 56.0
#> true MTD                   *     
#> true OBD                   *     
```
