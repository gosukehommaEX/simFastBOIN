## Submission

This is version 2.2.0 of simFastBOIN, an update of version 2.1.0 currently on
CRAN. It adds designs that decide from both the toxicity and the response and
select the optimal biological dose:

* the BOIN-ET design (Takeda et al., 2018) with its time-to-event version
  (TITE-BOIN-ET, Takeda et al., 2020) and its backfill version (BF-BOIN-ET,
  Takeda et al., 2025);
* the utility-based BOIN12 design (Lin et al., 2020) with its time-to-event
  version (TITE-BOIN12, Zhou et al., 2022);
* for two drugs given together, the Comb-BOIN12 design with its
  time-to-event version (Lu et al., 2025) and the BOIN-ETC design (Kakizume
  et al., 2024).

It also adds the accelerated BOIN design (Kojima et al., 2024) through a new
argument of the existing BOIN functions, functions that tabulate and plot the
results of several designs side by side, and two vignettes. Existing functions
only gain arguments whose defaults keep their results, so no existing result
changes.

The new engines are compiled code in src/ (boinet.cpp, boin12.cpp,
tite_boin12.cpp, comb_boin12.cpp, tite_comb_boin12.cpp and boin_etc.cpp, with
their headers), written in standard C++ with Rcpp as before. 'boinet',
'bfboinet' and 'mfp' are new suggested packages: the tests compare the
boundaries and the selected doses with 'boinet', the scripts in
inst/validation compare the simulations with 'boinet' and 'bfboinet', and
'mfp' fits an optional estimator after a check that it is installed. The
DESCRIPTION cites the articles by DOI.

The tests check the published decision tables and compare the engines trial
by trial with fixtures computed by separate implementations. The scripts in
inst/reproduce set the simulations beside the published tables of the
articles; their tables of judgments are installed in inst/reproduce/results
and inst/validation/results and summarized in the vignettes.

## Test environments

[[to be filled after the checks]]

## R CMD check results

[[to be filled after the checks]]

## Downstream dependencies

There are currently no downstream dependencies for this package.
