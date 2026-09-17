# simFastBOIN <img src="man/figures/logo.png" align="right" height="139" alt="" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/gosukehommaEX/simFastBOIN/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/gosukehommaEX/simFastBOIN/actions/workflows/R-CMD-check.yaml)
[![CRAN status](https://www.r-pkg.org/badges/version/simFastBOIN)](https://CRAN.R-project.org/package=simFastBOIN)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![CRAN downloads](https://cranlogs.r-pkg.org/badges/grand-total/simFastBOIN)](https://CRAN.R-project.org/package=simFastBOIN)
[![Monthly downloads](https://cranlogs.r-pkg.org/badges/simFastBOIN)](https://CRAN.R-project.org/package=simFastBOIN)
<!-- badges: end -->

Simulation tools for Bayesian optimal interval (BOIN) designs in phase I
dose-finding trials.

## Overview

simFastBOIN tabulates the BOIN decision boundaries, simulates the design to
obtain operating characteristics, and applies the MTD selection rule used at the
end of a trial. The simulation engine is written in C++.

The engine consumes exactly one uniform random variate per patient, drawn in
enrollment order, and applies the decision rules in the same order as
`BOIN::get.oc()`. With the same seed and matching arguments the two
implementations agree trial by trial, not merely on average. The test suite
checks this directly against the BOIN package rather than relying on a tolerance,
so the claim is verifiable rather than asserted.

## Installation

```r
install.packages("simFastBOIN")
```

Or the development version:

```r
# install.packages("devtools")
devtools::install_github("gosukehommaEX/simFastBOIN")
```

Version 1.4.0 and later contain compiled code, so a working C++ toolchain is
needed to install from source (Rtools on Windows, Xcode command line tools on
macOS).

## Quick start

```r
library(simFastBOIN)

oc <- sim_boin(
  target = 0.30,
  p_true = c(0.05, 0.15, 0.30, 0.45, 0.60),
  n_cohort = 10,
  cohort_size = 3,
  n_trials = 10000,
  seed = 123
)

oc
```

Compare several dose-toxicity scenarios under the same design:

```r
sim_boin_multi(
  target = 0.30,
  scenarios = list(
    "MTD at dose 2"   = c(0.15, 0.30, 0.45, 0.60, 0.75),
    "MTD at dose 4"   = c(0.02, 0.06, 0.15, 0.30, 0.50),
    "All doses toxic" = c(0.40, 0.55, 0.65, 0.75, 0.85)
  ),
  n_cohort = 10,
  cohort_size = 3,
  n_trials = 10000,
  seed = 123
)
```

Look at the boundaries the design will actually use:

```r
print(boin_boundary(target = 0.30, max_n = 18, extrasafe = TRUE), cohort_size = 3)
```

## Functions

**Design**

| Function | Purpose |
|---|---|
| `boin_lambda()` | Escalation and de-escalation interval boundaries |
| `boin_boundary()` | Integer decision boundaries by sample size |
| `boin_decision_table()` | Decision table indexed by DLTs and patients |

**Simulation**

| Function | Purpose |
|---|---|
| `sim_boin()` | Operating characteristics for one scenario |
| `sim_boin_multi()` | Operating characteristics across scenarios |
| `boin_simulate()` | Raw trial data from the engine |
| `boin_isotonic()` | Isotonic estimates of the dose-toxicity curve |
| `boin_select_mtd()` | MTD selection from completed trials |

## Design options

| Argument | Effect |
|---|---|
| `titration` | Treat one patient per dose until the first DLT, then switch to full cohorts |
| `extrasafe` | Add a stopping rule at the lowest dose that triggers before elimination |
| `bound_mtd` | Refuse to select a dose whose estimate exceeds the de-escalation boundary |
| `n_earlystop` | Stop once the current dose has accrued this many patients and the design would stay |
| `start_dose` | Dose level for the first cohort |
| `min_mtd_sample` | Smallest number of patients for a dose to be eligible as the MTD |

Note that `n_earlystop` defaults to 18 here, whereas `BOIN::get.oc()` defaults to
100, which in practice switches the rule off. Set it explicitly when comparing
the two.

## Why the trials stop

`boin_simulate()` records a reason for every trial, and `sim_boin()` reports the
distribution in `stop_reason_percent`.

| Reason | Meaning |
|---|---|
| `lowest_dose_eliminated` | The lowest dose was eliminated for toxicity |
| `lowest_dose_too_toxic` | The `extrasafe` rule stopped the trial at the lowest dose |
| `n_earlystop` | Enough patients had accrued at the current dose |
| `max_sample_size` | The maximum number of patients was reached |
| `max_cohorts` | All planned cohorts were completed |

## Checking the agreement with the BOIN package yourself

```r
reference <- BOIN::get.oc(
  target = 0.30, p.true = c(0.05, 0.15, 0.25, 0.45, 0.60),
  ncohort = 20, cohortsize = 3, n.earlystop = 18, ntrial = 10000, seed = 6
)

ours <- sim_boin(
  target = 0.30, p_true = c(0.05, 0.15, 0.25, 0.45, 0.60),
  n_cohort = 20, cohort_size = 3, n_earlystop = 18, n_trials = 10000, seed = 6
)

all.equal(unname(ours$sel_percent), reference$selpercent)
all.equal(unname(ours$n_pts_dose), reference$npatients)
all.equal(unname(ours$n_tox_dose), reference$ntox)
all.equal(ours$percent_no_mtd, reference$percentstop)
```

## Speed

Five doses, 20 cohorts of three, `n_earlystop = 18` and 10,000 simulated trials,
measured on one Windows machine:

| | Elapsed |
|---|---|
| `BOIN::get.oc()` | 9.98 s |
| `sim_boin()` | 0.07 s |

That is about two orders of magnitude for this design. The ratio depends on the
design and on the machine, so the code below re-measures it rather than asking
you to take the table on trust.

```r
scenario <- c(0.05, 0.15, 0.25, 0.45, 0.60)

system.time(
  BOIN::get.oc(target = 0.30, p.true = scenario, ncohort = 20, cohortsize = 3,
               n.earlystop = 18, ntrial = 10000, seed = 6)
)

system.time(
  sim_boin(target = 0.30, p_true = scenario, n_cohort = 20, cohort_size = 3,
           n_earlystop = 18, n_trials = 10000, seed = 6)
)
```

## Upgrading from version 1.3.2

Version 1.4.0 renames most functions, changes the argument order of `sim_boin()`
and corrects several defects that affected numerical results. Simulations run
with the same seed will not reproduce the version 1.3.2 output. See
[NEWS.md](NEWS.md) for the full list and the mapping of old names to new ones.

## References

Liu, S. and Yuan, Y. (2015). Bayesian Optimal Interval Designs for Phase I
Clinical Trials. *Journal of the Royal Statistical Society: Series C*, 64,
507-523.

Yan, F., Zhang, L., Zhou, Y., Pan, H., Liu, S. and Yuan, Y. (2020). BOIN: An R
Package for Designing Single-Agent and Drug-Combination Dose-Finding Trials Using
Bayesian Optimal Interval Designs. *Journal of Statistical Software*, 94(13),
1-32.

## License

MIT (c) Gosuke Homma
