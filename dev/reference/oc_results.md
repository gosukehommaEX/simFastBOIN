# Operating Characteristics of Several Designs in One Table

Collect the results of the simulation functions of different designs
(and of
[`oc_3p3`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/oc_3p3.md))
into one long data frame, so that the designs can be tabulated and drawn
side by side. The result has
[`print()`](https://rdrr.io/r/base/print.html) and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods;
[`oc_compare`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/oc_compare.md)
judges the same results against the true doses.

## Usage

``` r
oc_results(...)
```

## Arguments

- ...:

  Design results, named by the design, for example
  `oc_results(BOIN12 = res_a, "BOIN-ET" = res_b)`. Each is the result of
  a simulation function for one scenario (such as
  [`sim_boin()`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_boin.md)
  or
  [`sim_boin12()`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/sim_boin12.md)),
  of its `_multi` version, of
  [`oc_3p3`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/oc_3p3.md),
  or a list of such single-scenario results named by the scenarios (the
  way to give
  [`oc_3p3()`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/oc_3p3.md)
  several scenarios). Unnamed results are named after their design. All
  results must be for a single agent or all for a drug combination, with
  the same number of doses.

## Value

A data frame of class `c("oc_results", "data.frame")` with the columns

- design:

  Factor of the design names, in the order given.

- target_type:

  `"MTD"` or `"OBD"`, what the design selects.

- scenario:

  Factor of the scenario names. A single-scenario result is named
  "Scenario 1".

- dose:

  Dose level, or `dose_a` and `dose_b` for a combination; `NA` for the
  quantities of the whole trial.

- metric:

  Factor: per dose, `"sel_percent"` (percentage of trials selecting the
  dose), `"n_pts"`, `"n_tox"` and `"n_resp"` (average numbers of
  patients, DLTs and responses); per trial, `"percent_no_selection"`
  (percentage of trials selecting no dose), `"total_n_pts"`,
  `"total_n_tox"`, `"total_n_resp"` and `"duration"` (average trial
  duration, for the designs with time). Quantities a design does not
  record are left out.

- value:

  The value.

with the attributes `type` (`"single"` or `"combination"`) and `dims`.

## Details

The values are those of the components `sel_percent`, `n_pts_dose`,
`n_tox_dose`, `n_resp_dose`, `percent_no_mtd` or `percent_no_obd`,
`total_n_pts`, `total_n_tox`, `total_n_resp` and `duration_mean` or
`duration` of the results; nothing is simulated again.

## See also

[`oc_scenarios`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/oc_scenarios.md),
[`oc_compare`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/oc_compare.md),
[`print.oc_results`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/print.oc_results.md),
[`plot.oc_results`](https://gosukehommaex.github.io/simFastBOIN/dev/reference/plot.oc_results.md)

## Examples

``` r
scen <- list(
  "Scenario 1" = list(p_true = c(0.05, 0.12, 0.27, 0.35, 0.50),
                      p_resp = c(0.20, 0.35, 0.36, 0.37, 0.38))
)
boin12 <- sim_boin12_multi(target = 0.35, target_resp = 0.25,
                           scenarios = scen, n_cohort = 12,
                           cohort_size = 3, n_trials = 200, seed = 1)
boin <- sim_boin_multi(target = 0.35, scenarios = scen, n_cohort = 12,
                       cohort_size = 3, n_trials = 200, seed = 1)
res <- oc_results(BOIN12 = boin12, BOIN = boin)
head(res)
#> Operating characteristics of 1 design(s) in 1 scenario(s)
#> 
#> Scenario 1
#> Selection (%) by dose
#>           1    2    3    4   5
#> BOIN12 16.5 49.5 22.5 10.5 1.0
#> Patients by dose
#>          1
#> BOIN12 6.8
res
#> Operating characteristics of 2 design(s) in 1 scenario(s)
#> 
#> Scenario 1
#> Selection (%) by dose
#>           1    2    3    4    5
#> BOIN12 16.5 49.5 22.5 10.5  1.0
#> BOIN    0.0  9.0 38.5 39.0 13.5
#> Patients by dose
#>          1    2    3   4   5
#> BOIN12 6.8 12.6  8.7 5.7 2.1
#> BOIN   4.0  7.1 10.7 8.5 3.3
#> Per trial (averages)
#>        No selection (%) Patients DLTs Responses
#> BOIN12              0.0     36.0  7.5      11.8
#> BOIN                0.0     33.6  8.7         -
```
