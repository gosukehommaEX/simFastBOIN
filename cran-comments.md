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

* local: Windows 11 x64, R 4.6.1
* win-builder: R-devel (2026-11-17 r90926 ucrt) and R 4.6.2 (2026-10-29 ucrt)
* GitHub Actions: macOS-latest (release), windows-latest (release),
  ubuntu-latest (devel, release, oldrel-1)

## R CMD check results

0 errors | 0 warnings | 1 note

The note is raised by the incoming feasibility check on win-builder, under both
R-devel and R-release:

    Possibly misspelled words in DESCRIPTION:
      Kakizume (42:55)
      Kojima (44:73)
      Takeda (33:8, 35:6, 36:14)
      Zhou (38:68)
      backfill (35:66)
      titration (43:49)

The words are spelled correctly. Kakizume, Kojima, Takeda and Zhou are family
names of the authors of the cited articles. To backfill is to treat the
patients who arrive while the dose escalation waits at lower doses, and
titration is the escalation by single patients at the lowest doses before
full cohorts start, as in the 'BOIN' package; both are terms of dose-finding
trials.

The local check with remote = TRUE gave no errors, warnings or notes, the URLs
were checked with urlchecker, and the R CMD check workflow, which fails on
warnings, passed on all five GitHub Actions platforms.

## Downstream dependencies

There are currently no downstream dependencies for this package.
