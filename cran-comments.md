## Submission

This is version 2.1.0 of simFastBOIN, an update of version 2.0.0 currently on
CRAN. It adds the time-to-event BOIN (TITE-BOIN) design of Yuan et al. (2018)
and Lin and Yuan (2020), for toxicity that is observed late in the assessment
window: a decision table, a simulation engine and summaries of the operating
characteristics. It also adds the BF-BOIN design of Zhao et al. (2024) and the
BE-BOIN design of Chen et al. (2026), which backfill patients to lower doses
while the dose escalation waits, the drug-combination BOIN design of Lin and
Yin (2017), with a function that gives the combination for the next cohort of
a trial in progress, and a decision table for the 3+3 design. No existing
function or result changes.

The new engines are compiled code in src/tite.cpp, src/tite_core.h,
src/backfill.cpp, src/backfill_core.h, src/comb.cpp and src/comb_core.h,
written in standard C++ with Rcpp as before. The package 'Iso', previously
suggested, is now imported: the drug-combination design selects the MTD with
its bivariate isotonic regression, as the reference implementation in 'BOIN'
does. The DESCRIPTION cites the TITE-BOIN, backfilling and drug-combination
articles by DOI.

The TITE-BOIN decision table reproduces every entry of the published tables of
Yuan et al. (2018) and Chen et al. (2025, 2026), and the BF-BOIN simulations
reproduce the operating characteristics of Table 4 of Zhao et al. (2024)
within simulation error; the tests check both. Whenever no patient is pending
at a decision, or no dose is opened for backfilling, the simulated trials are
identical to those of the existing engines under the same seed, which the
tests also check. The drug-combination simulations reproduce the published
example of Yan et al. (2020) and agree exactly with BOIN::get.oc.comb() under
the same seed; the tests check both. The function for a trial in progress,
comb_boin_next(), takes its decision with the same compiled code as these
simulations, and the tests check that it reproduces the decision of the
simulated trials after every cohort.

## Test environments

* Local Windows, R 4.6.0 (2026-04-24 ucrt): `devtools::check(remote = TRUE)`
* GitHub Actions: macOS (R release), Windows (R release), Ubuntu (R devel,
  release and oldrel-1)
* win-builder: R Under development (unstable) (2026-10-10 r90659 ucrt)
* win-builder: R version 4.6.1 (2026-06-24 ucrt)

## R CMD check results

0 errors | 0 warnings | 0 notes locally. The checks on GitHub Actions also
pass.

On win-builder (R-devel and R release) there is 1 NOTE:

    Possibly misspelled words in DESCRIPTION:
      Zhao (26:8)
      backfilled (25:38)

Both are spelled correctly: Zhao is the surname of the first author of a cited
article, and "backfilled" is the term the cited articles use for patients
treated at a lower dose while the dose escalation waits.

## Downstream dependencies

There are currently no downstream dependencies for this package.
