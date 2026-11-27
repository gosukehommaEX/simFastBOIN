# Computation time

The simulation engines of the package are written in C++. This article
gives the time that each design takes to simulate the trials of one
scenario and to summarize their operating characteristics, on one
machine. The times are read from `inst/benchmark/results/benchmark.csv`,
written by `inst/benchmark/benchmark.R`.

## What is measured

The 10 designs that have a reference implementation in another R package
are run by both with the same settings and the same number of trials:
those of the comparison of the design in `inst/validation/`, or of the
README for the BOIN design, with the number of trials used there for the
reference package. The 11 designs without a reference implementation are
run by the package alone, with the settings of the example in their
vignette and 10,000 trials. The settings are listed in the section
“Settings” below.

The time is the elapsed time of one call. A fast call is repeated until
the repetitions take at least one second, and the time per call is kept.
Every time was measured 3 times, alternating between the two
implementations, and the tables give the median. Printed output and
warnings of the calls are discarded. The package was installed before it
was timed: `devtools::load_all()` compiles the C++ code without
optimization.

## The machine

The times were measured on 2026-11-26 under Windows 11 x64 (build
26200), on an Intel(R) Core(TM) Ultra 7 155H with 16 physical cores,
with R 4.6.1. The simulations of the package run on one core. The
versions of the packages were

| Package            | Version    |
|:-------------------|:-----------|
| simFastBOIN        | 2.2.0.9000 |
| BOIN               | 2.7.2      |
| boinet             | 1.6.0      |
| bfboinet           | 1.1.0      |
| UnifiedDoseFinding | 0.1.10     |
| escalation         | 0.2.3      |

## Designs with a reference implementation

| Design | Trials | Reference | Reference (s) | simFastBOIN | simFastBOIN (s) | Times faster |
|:---|:---|:---|:---|:---|:---|:---|
| BOIN | 10,000 | BOIN::get.oc() | 12 | sim_boin() | 0.044 | 260 |
| BOIN combination | 1,000 | BOIN::get.oc.comb() | 0.8 | sim_comb_boin() | 0.13 | 5.9 |
| Waterfall | 500 | BOIN::get.oc.comb(mtd.contour = TRUE) | 2 | sim_waterfall() | 0.08 | 25 |
| BOIN-ET | 5,000 | boinet::boinet() | 25 | sim_boin_et() | 0.044 | 560 |
| TITE-BOIN-ET | 2,000 | boinet::tite.boinet() | 58 | sim_tite_boin_et() | 0.061 | 960 |
| gBOIN-ET | 5,000 | boinet::gboinet() | 25 | sim_gboin_et() | 0.14 | 170 |
| TITE-gBOIN-ET | 5,000 | boinet::tite.gboinet() | 150 | sim_tite_gboin_et() | 0.4 | 390 |
| BF-BOIN-ET | 1,000 | bfboinet::get.oc.backboinet() | 340 | sim_bf_boin_et() | 0.032 | 11,000 |
| gBOIN | 4,000 | UnifiedDoseFinding::get_oc_QuasiBOIN() | 6.8 | sim_gboin() | 0.021 | 320 |
| BOIN12 | 2,000 | escalation::simulate_trials() | 140 | sim_boin12() | 0.042 | 3,300 |

The package was between 5.9 and 11,000 times faster than the reference
implementation. The ratio is smallest for the two drug-combination
designs, 5.9 for the BOIN drug-combination design and 25 for the
waterfall design: both packages select the MTD of every trial with the
bivariate isotonic regression of
[`Iso::biviso()`](https://rdrr.io/pkg/Iso/man/biviso.html), called from
R once per trial. It is largest for BF-BOIN-ET, whose reference
implementation takes 340 seconds for 1,000 trials.

## Designs without a reference implementation

| Design                | Trials | simFastBOIN            | Seconds |
|:----------------------|:-------|:-----------------------|:--------|
| iBOIN                 | 10,000 | sim_iboin()            | 0.053   |
| Accelerated BOIN      | 10,000 | sim_boin()             | 0.039   |
| TITE-BOIN, imputation | 10,000 | sim_tite_boin()        | 0.074   |
| TITE-BOIN, ESS        | 10,000 | sim_tite_boin()        | 0.093   |
| BF-BOIN               | 10,000 | sim_bf_boin()          | 0.16    |
| BE-BOIN               | 10,000 | sim_be_boin()          | 0.17    |
| TITE-BOIN12           | 10,000 | sim_tite_boin12()      | 0.84    |
| U-BOIN                | 10,000 | sim_uboin()            | 0.19    |
| Comb-BOIN12           | 10,000 | sim_comb_boin12()      | 96      |
| TITE-Comb-BOIN12      | 10,000 | sim_tite_comb_boin12() | 94      |
| BOIN-ETC              | 10,000 | sim_boin_etc()         | 0.12    |

9 of the 11 designs take less than one second for 10,000 trials.
Comb-BOIN12 and TITE-Comb-BOIN12 take about 96 and 94 seconds: by
default (`est_method = "logistic"`) they estimate the response
probabilities at the end of every trial by the posterior mean under the
logistic model of Lu et al. (2025), computed by importance sampling with
4,000 points (`n_draws = 4000`). With `est_method = "observed"`, the
observed response rates, the selection is much faster.

## Settings

| Design | Setting |
|:---|:---|
| BOIN | Five doses, 20 cohorts of three, n_earlystop = 18 (README) |
| BOIN combination | Lin and Yin 1 (3 x 5), 20 cohorts of three, n_earlystop = 12 (compare-comb-with-BOIN.R) |
| Waterfall | Scenario 9 of Zhang and Yuan (2016) (3 x 5), cohorts of three, boin_compatible = TRUE (compare-waterfall-with-BOIN.R) |
| BOIN-ET | Scenario 1 of Takeda et al. (2018), 12 cohorts of three (compare-with-boinet.R) |
| TITE-BOIN-ET | Scenario 1 of Takeda et al. (2018), windows of 30 and 45 days, a patient every 10 days (compare-with-boinet.R) |
| gBOIN-ET | Scenario 1 of Takeda et al. (2022), probabilities of Table S1, 10 cohorts of three (compare-with-gboinet.R) |
| TITE-gBOIN-ET | Scenario 1 of Takeda et al. (2023), windows of 30 and 60 days, 12 cohorts of three (compare-with-tite-gboinet.R) |
| BF-BOIN-ET | Example of ‘bfboinet’, 10 cohorts of three, windows of one month, three patients a month (compare-with-bfboinet.R) |
| gBOIN | Scenario 1 of Mu et al. (2019), 10 cohorts of three, elim_rule = “count” (compare-with-UnifiedDoseFinding.R) |
| BOIN12 | Scenario 1 of Lin et al. (2020), 12 cohorts of three (compare-with-escalation.R) |
| iBOIN | Scenario 3 of Zhou et al. (2021), prior_ess = 3, 10 cohorts of three (vignette) |
| Accelerated BOIN | titration_doses = 2, n_earlystop = 6, eight cohorts of three (vignette) |
| TITE-BOIN, imputation | Window of three months, two patients a month, 10 cohorts of three (vignette) |
| TITE-BOIN, ESS | Window of three months, two patients a month, 10 cohorts of three (vignette) |
| BF-BOIN | Setting of Zhao et al. (2024), window of one month, three patients a month (vignette) |
| BE-BOIN | Scenario of Zhao et al. (2024), window of three months, two patients a month (vignette) |
| TITE-BOIN12 | Scenario 1 of Lin et al. (2020), windows of 45 and 60 days, three patients every 30 days, 12 cohorts of three (vignette) |
| U-BOIN | Scenario 1 of Simulation A of Zhou, Lee and Yuan (2019), 18 cohorts of three (vignette) |
| Comb-BOIN12 | Scenario 16 of Lu et al. (2025) (2 x 3), 14 cohorts of three (vignette) |
| TITE-Comb-BOIN12 | Scenario 16 of Lu et al. (2025), windows of one and two months, two patients a month (vignette) |
| BOIN-ETC | Scenario 1 of Kakizume et al. (2024) (2 x 3), eight cohorts of three (vignette) |

The scripts named in the settings are in `inst/validation/`, and the
scenarios of the articles are taken from the scripts of
`inst/reproduce/`.

## Measuring on your own machine

From the root of the source package, install the package, restart R and
run

``` r

library(simFastBOIN)
source("inst/benchmark/benchmark.R")
res <- run_benchmark()
```

with the reference packages ‘BOIN’, ‘boinet’, ‘bfboinet’,
‘UnifiedDoseFinding’ and ‘escalation’ installed; a design whose package
is missing is recorded as such. The run writes
`inst/benchmark/results/benchmark.csv`. Most of the time goes to the
reference packages.

## References

Lu, M., Zhang, J., Yuan, Y. and Lin, R. (2025). Comb-BOIN12: A
utility-based Bayesian optimal interval design for dose optimization in
cancer drug-combination trials. *Statistics in Biopharmaceutical
Research*, 17(2), 266-276.
