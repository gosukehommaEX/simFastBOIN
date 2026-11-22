# Compare Several Designs Against the True Doses

Judge the results of several designs against the true doses of the
scenarios: how often each design selects the true MTD (designs that
select an MTD) or the true OBD (designs that select an OBD), how many
patients it treats there, how often it selects an overdose and how many
patients it treats at overdoses, how often it selects no dose, and its
average sample size and duration. The result has
[`print()`](https://rdrr.io/r/base/print.html) and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods.

## Usage

``` r
oc_compare(..., truth, overdose_cutoff = NULL)
```

## Arguments

- ...:

  Design results, named by the design, as for
  [`oc_results`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/oc_results.md).

- truth:

  The result of
  [`oc_scenarios`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/oc_scenarios.md)
  for the same scenarios, with the same scenario names. The DLT
  probabilities of every scenario must agree with those stored in the
  results. Designs that select an OBD need scenarios with response
  probabilities.

- overdose_cutoff:

  Numeric scalar or `NULL`. A dose whose true DLT probability is above
  this value is an overdose. Defaults to `NULL`, the `tox_limit` of
  `truth`.

## Value

A data frame of class `c("oc_compare", "data.frame")` with one row per
design and scenario and the columns

- design, target_type, scenario:

  As in
  [`oc_results`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/oc_results.md).

- pct_correct:

  Percentage of trials selecting a true MTD or OBD; `NA` when the
  scenario has none.

- n_pts_correct, pct_patients_correct:

  Average number of patients treated at the true MTDs or OBDs, and its
  percentage of all the patients.

- pct_overdose_selected:

  Percentage of trials selecting an overdose.

- n_pts_overdose, pct_patients_overdose:

  Average number of patients treated at overdoses, and its percentage of
  all the patients.

- pct_no_selection:

  Percentage of trials selecting no dose.

- total_n_pts:

  Average number of patients per trial.

- duration:

  Average trial duration, `NA` for the designs without time.

with the attribute `overdose_cutoff`.

## Details

Every quantity is computed from the per-dose averages of the results
(`sel_percent` and `n_pts_dose`), so nothing is simulated again and the
results need not keep the trials. A dose is an overdose when its DLT
probability exceeds `overdose_cutoff` by more than `1e-12`.

## See also

[`oc_scenarios`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/oc_scenarios.md),
[`oc_results`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/oc_results.md),
[`print.oc_compare`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/print.oc_compare.md),
[`plot.oc_compare`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/plot.oc_compare.md)

## Examples

``` r
scen <- list(
  "Scenario 1" = list(p_true = c(0.05, 0.12, 0.27, 0.35, 0.50),
                      p_resp = c(0.20, 0.35, 0.36, 0.37, 0.38)),
  "Scenario 2" = list(p_true = c(0.05, 0.12, 0.27, 0.35, 0.50),
                      p_resp = c(0.01, 0.05, 0.30, 0.60, 0.60))
)
truth <- oc_scenarios(scen, target = 0.35, target_resp = 0.25)
boin12 <- sim_boin12_multi(target = 0.35, target_resp = 0.25,
                           scenarios = scen, n_cohort = 12,
                           cohort_size = 3, n_trials = 200, seed = 1)
boin <- sim_boin_multi(target = 0.35, scenarios = scen, n_cohort = 12,
                       cohort_size = 3, n_trials = 200, seed = 1)
oc_compare(BOIN12 = boin12, BOIN = boin, truth = truth)
#> Comparison of 2 design(s) with the true doses
#> Correct: the true MTD for designs selecting an MTD, the true OBD for designs selecting an OBD
#> Overdose: true DLT probability above 0.35 
#> 
#> Scenario 1
#>              Correct (%) Pts at correct (%) Overdose sel. (%) Pts at overdose
#> BOIN12 (OBD)        49.5               35.1               1.0             2.1
#> BOIN (MTD)          39.0               25.2              13.5             3.3
#>              No selection (%) Patients
#> BOIN12 (OBD)              0.0     36.0
#> BOIN (MTD)                0.0     33.6
#> 
#> Scenario 2
#>              Correct (%) Pts at correct (%) Overdose sel. (%) Pts at overdose
#> BOIN12 (OBD)        49.5               30.2               5.0             3.0
#> BOIN (MTD)          39.0               25.2              13.5             3.3
#>              No selection (%) Patients
#> BOIN12 (OBD)              8.0     35.2
#> BOIN (MTD)                0.0     33.6
```
