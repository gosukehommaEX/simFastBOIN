## Submission

This is version 2.1.0 of simFastBOIN, an update of version 2.0.0 currently on
CRAN. It adds the time-to-event BOIN (TITE-BOIN) design of Yuan et al. (2018)
and Lin and Yuan (2020), for toxicity that is observed late in the assessment
window: a decision table, a simulation engine and summaries of the operating
characteristics. It also adds the BF-BOIN design of Zhao et al. (2024) and the
BE-BOIN design of Chen et al. (2026), which backfill patients to lower doses
while the dose escalation waits, and a decision table for the 3+3 design. No
existing function or result changes.

The new engines are compiled code in src/tite.cpp, src/tite_core.h,
src/backfill.cpp and src/backfill_core.h, written in standard C++ with Rcpp as
before; no dependency was added. The DESCRIPTION cites the TITE-BOIN and
backfilling articles by DOI.

The TITE-BOIN decision table reproduces every entry of the published tables of
Yuan et al. (2018) and Chen et al. (2025, 2026), and the BF-BOIN simulations
reproduce the operating characteristics of Table 4 of Zhao et al. (2024)
within simulation error; the tests check both. Whenever no patient is pending
at a decision, or no dose is opened for backfilling, the simulated trials are
identical to those of the existing engines under the same seed, which the
tests also check.

## Test environments

TODO

## R CMD check results

TODO

## Downstream dependencies

TODO
