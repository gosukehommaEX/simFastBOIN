## Submission

This is version 2.2.0 of simFastBOIN, an update of version 2.1.0 currently on
CRAN. It adds the BOIN-ET design of Takeda et al. (2018), which decides from
both the toxicity and the response and selects the optimal biological dose,
with its time-to-event version (TITE-BOIN-ET, Takeda et al., 2020) and its
backfill version (BF-BOIN-ET, Takeda et al., 2025): the boundaries, a decision
table, simulation engines, the selection of the optimal biological dose and a
new vignette. No existing function or result changes.

The new engines are compiled code in src/boinet.cpp, src/boinet_core.h,
src/tite_boinet_core.h and src/bf_boinet_core.h, written in standard C++ with
Rcpp as before. 'boinet', 'bfboinet' and 'mfp' are new suggested packages: the
tests compare the boundaries and the selected doses with 'boinet', the scripts
in inst/validation compare the simulations with 'boinet' and 'bfboinet', and
'mfp' fits the optional fractional polynomial estimator after a check that it
is installed. The DESCRIPTION cites the three articles by DOI.

The boundaries and the decision table reproduce those published by Takeda et
al. (2018) exactly, and the boundaries and the doses selected from simulated
trials agree exactly with boinet::gridoptim() and boinet::obd.select(); the
tests check both. The scripts in inst/reproduce set the simulations beside the
published tables of the three articles.

## Test environments

[[to be filled after the checks]]

## R CMD check results

[[to be filled after the checks]]

## Downstream dependencies

There are currently no downstream dependencies for this package.
