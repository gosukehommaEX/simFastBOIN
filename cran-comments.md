## Submission

This is version 2.1.0 of simFastBOIN, an update of version 2.0.0 currently on
CRAN. It adds the time-to-event BOIN (TITE-BOIN) design of Yuan et al. (2018)
and Lin and Yuan (2020), for toxicity that is observed late in the assessment
window: a decision table, a simulation engine and summaries of the operating
characteristics. It also adds a decision table for the 3+3 design. No existing
function or result changes.

The new engine is compiled code in src/tite.cpp and src/tite_core.h, written in
standard C++ with Rcpp as before; no dependency was added. The DESCRIPTION cites
the two TITE-BOIN articles by DOI.

The decision table reproduces every entry of the two published tables of Yuan
et al. (2018), 640 states in all, and the tests check this. Whenever no patient
is pending at a decision, the simulated trials are identical to those of the
existing BOIN engine under the same seed, which the tests also check.

## Test environments

TODO

## R CMD check results

TODO

## Downstream dependencies

TODO
