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

The engine draws one uniform variate per patient, in enrollment order, and
applies the decision rules in the same order as
`BOIN::get.oc()`. With the same seed and matching arguments the two
implementations agree trial by trial, not merely on average. The test suite
checks this directly against the BOIN package rather than relying on a tolerance,
so the claim is verifiable rather than asserted. A wider comparison ships with the
package as `inst/validation/compare-with-BOIN.R`, which runs 200 configurations
covering every combination of the options the two packages share, including
`extrasafe` and `titration`.

Version 2.1.0 adds the time-to-event BOIN (TITE-BOIN) design, for toxicity that
appears late in the assessment window or accrual that outpaces it. New patients
are treated while earlier ones are still being followed, with the pending
outcomes handled either by the single mean imputation of Yuan et al. (2018) or
by the effective sample size of Lin and Yuan (2020). The decision table
reproduces the published tables of Yuan et al. (2018) entry by entry, and
whenever no patient is pending the simulated trials are identical to those of
`sim_boin()` with the same seed.

## Installation

```r
install.packages("simFastBOIN")
```

Or the development version:

```r
# install.packages("devtools")
devtools::install_github("gosukehommaEX/simFastBOIN")
```

Version 2.0.0 contains compiled code, so a working C++ toolchain is
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

Read the dose decisions off a table, or off a figure for a protocol:

```r
decisions <- boin_decision_table(target = 0.30, max_n = 18)

print(decisions, cohort_size = 3)
plot(decisions)                      # needs ggplot2
```

When toxicity can appear late, or patients arrive faster than they can be
assessed, use the time-to-event design. Its decisions also depend on how long the
pending patients have been followed:

```r
tite <- tite_boin_decision_table(target = 0.30, max_n = 15)
print(tite, cohort_size = 3)

sim_tite_boin(
  target = 0.30,
  p_true = c(0.05, 0.15, 0.30, 0.45, 0.60),
  n_cohort = 10,
  cohort_size = 3,
  window = 3,          # a three month assessment window
  accrual_rate = 2,    # two patients a month
  n_trials = 10000,
  seed = 123
)
```

## Functions

**Design**

| Function | Purpose |
|---|---|
| `boin_lambda()` | Escalation and de-escalation interval boundaries |
| `boin_boundary()` | Integer decision boundaries by sample size |
| `boin_decision_table()` | Decision table indexed by DLTs and patients |
| `boin_stopping_table()` | Safety stopping boundary as a two-row table |

**Simulation**

| Function | Purpose |
|---|---|
| `sim_boin()` | Operating characteristics for one scenario |
| `sim_boin_multi()` | Operating characteristics across scenarios |
| `boin_simulate()` | Raw trial data from the engine |
| `boin_isotonic()` | Isotonic estimates of the dose-toxicity curve |
| `boin_select_mtd()` | MTD selection from completed trials |

**Time-to-event BOIN (TITE-BOIN)**

| Function | Purpose |
|---|---|
| `tite_boin_decision_table()` | Decision table with pending patients, for either method |
| `sim_tite_boin()` | Operating characteristics, trial duration and suspensions of accrual |
| `sim_tite_boin_multi()` | The same across scenarios |
| `tite_boin_simulate()` | Raw trial data from the time-to-event engine |

**3+3 design**

| Function | Purpose |
|---|---|
| `oc_3p3()` | Operating characteristics in closed form |
| `sim_3p3()` | The same by simulation, as a check on the closed form |
| `decision_table_3p3()` | Decision table of the design |

## Design options

| Argument | Effect |
|---|---|
| `titration` | Treat one patient per dose until the first DLT, then switch to full cohorts |
| `extrasafe` | Add a stopping rule at the lowest dose that triggers before elimination |
| `bound_mtd` | Refuse to select a dose whose estimate exceeds the de-escalation boundary |
| `n_earlystop` | Stop once the current dose has accrued this many patients and the design would stay |
| `start_dose` | Dose level for the first cohort, ignored under `titration` |
| `stay_on_1_of_3` | Make one DLT out of three a stay rather than a de-escalation |
| `min_mtd_sample` | Smallest number of patients for a dose to be eligible as the MTD |

Note that `n_earlystop` defaults to 18 here, whereas `BOIN::get.oc()` defaults to
100, which in practice switches the rule off. Set it explicitly when comparing
the two.

## Late-onset toxicity: the TITE-BOIN design

| `method` | Pending outcomes handled by | Boundaries refer to | Default suspension rule |
|---|---|---|---|
| `"imputation"` | Single mean imputation (Yuan et al., 2018) | STFT | Suspend accrual when more than half of the patients at the dose are pending |
| `"ess"` | Effective sample size (Lin and Yuan, 2020) | ESS | Escalate only once two patients at the dose have completed the assessment |

STFT is the total follow-up time of the pending patients divided by the length of
the assessment window, and ESS adds to it the number of patients who have
completed the assessment. Neither table depends on the length of the window.
Both suspension rules can be set for either method with `max_pending_ratio` and
`min_completed`.

A simulation needs the length of the window and the accrual rate in the same
unit of time. Arrivals can be exponential (the default), uniform or fixed, and
the time to DLT Weibull, with a share `late_fraction` of the DLTs in the second
half of the window, or uniform. The follow-up of the pending patients can be
weighted by a piecewise uniform prior for the time to DLT through
`prior_weights`; both articles evaluate their designs with equal weights, the
default.

How the implementation is checked:

* `tite_boin_decision_table()` reproduces all 640 entries of Table 1 (target 0.2)
  of Yuan et al. (2018) and of Table S1 (target 0.3) of its supplementary
  appendix, and the dose decisions of the trial example in Figure 1 of that
  article. The published tables allow de-escalation when the observed DLT rate
  equals the target, which equation (5) of the appendix, read literally, does
  not; the package follows the tables.
* With no patient pending at any decision, the trials of `tite_boin_simulate()`
  and the summaries of `sim_tite_boin()` are identical to those of
  `boin_simulate()` and `sim_boin()` under the same seed.
* The simulated trials agree exactly with a separate implementation written in
  Python and run on a replica of R's random number stream.

Neither article reports operating characteristics that can be reproduced
number for number: Yuan et al. (2018) give them as figures relative to the 3+3
design and as averages over randomly generated scenarios, and the tables of Lin
and Yuan (2020) are for the keyboard and mTPI designs.

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

The time-to-event engine was timed on a similar design, five doses and 20
cohorts of three with the true DLT rates `c(0.05, 0.15, 0.30, 0.45, 0.60)`, a
three month window and two patients a month, 10,000 trials, all three runs
together on one Windows machine:

| | Elapsed |
|---|---|
| `sim_boin()` | 0.05 s |
| `sim_tite_boin()`, `method = "imputation"` | 0.09 s |
| `sim_tite_boin()`, `method = "ess"` | 0.11 s |

Following patients over time costs about twice as much as the BOIN engine. The
BOIN package has no time-to-event design, so there is no reference
implementation of TITE-BOIN to time it against.

## Comparing with the 3+3 design

The traditional 3+3 design is included as a comparator. Its operating
characteristics are obtained in closed form rather than by simulation, so they
carry no Monte Carlo error.

```r
oc_3p3(p_true = c(0.30, 0.48, 0.67))

# The other definition of the MTD in common use
oc_3p3(p_true = c(0.30, 0.48, 0.67), mtd_rule = "expand")
```

`sim_3p3()` simulates the same design and exists to confirm the closed form. The
result of either has the same component names as the result of `sim_boin()`, so
the two designs can be tabulated side by side.

`decision_table_3p3()` shows the rule itself, with `print()` and `plot()`
methods like those of the BOIN decision table.

```r
decision_table_3p3(mtd_rule = "expand")
```

## Upgrading from version 1.3.2

Version 2.0.0 renames most functions, changes the argument order of `sim_boin()`
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

Yuan, Y., Lin, R., Li, D., Nie, L. and Warren, K. E. (2018). Time-to-Event
Bayesian Optimal Interval Design to Accelerate Phase I Trials. *Clinical Cancer
Research*, 24(20), 4921-4930.

Lin, R. and Yuan, Y. (2020). Time-to-Event Model-Assisted Designs for
Dose-Finding Trials with Delayed Toxicity. *Biostatistics*, 21(4), 807-824.

## License

MIT (c) Gosuke Homma
