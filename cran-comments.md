## Submission

This is version 2.3.0 of simFastBOIN, an update of version 2.2.0 currently on
CRAN. It adds six designs:

* the iBOIN design (Zhou et al., 2021), which brings prior information on
  the DLT probabilities into the decision boundaries;
* the gBOIN design (Mu et al., 2019), which decides from toxicity grades;
* the gBOIN-ET design (Takeda et al., 2022) and its time-to-event version
  TITE-gBOIN-ET (Takeda et al., 2023), which select the optimal biological
  dose from graded toxicity and graded efficacy;
* for two drugs given together, the waterfall design (Zhang and Yuan, 2016),
  which finds the contour of maximum tolerated doses;
* the U-BOIN design (Zhou, Lee and Yuan, 2019), which selects the optimal
  biological dose by the utility of the joint outcome of a patient.

It also adds a vignette that sets designs and scenarios side by side.
Existing functions keep their results, except that BOIN-ET, TITE-BOIN-ET
and BF-BOIN-ET now compare the utilities of the doses after rounding, so
that a tie between two doses no longer depends on rounding error and the
selection is the same on every platform (NEWS.md).

The new engines are compiled code in src/ (iboin.cpp, gboin.cpp,
gboinet.cpp, tite_gboinet.cpp, waterfall.cpp and uboin.cpp, with their
headers), written in standard C++ with Rcpp as before. There are no new
dependencies. The DESCRIPTION cites the new articles by DOI.

As in the previous versions, the tests check the published tables and
compare the engines trial by trial with fixtures computed by separate
implementations, and the scripts in inst/reproduce and inst/validation set
the simulations beside the published tables and the reference
implementations; their tables of judgments are installed in
inst/reproduce/results and inst/validation/results and summarized in the
vignettes. A script in inst/benchmark times the package against the
reference implementations, with its results in inst/benchmark/results.

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
      Zhang (36:66)
      gBOIN (24:17, 43:47)
      iBOIN (22:57)

The words are spelled correctly. Zhang is the family name of the first author
of a cited article (Zhang and Yuan, 2016). iBOIN and gBOIN are the names that
the cited articles give to their designs, the BOIN design with prior
information (Zhou et al., 2021) and the generalized BOIN design (Mu et al.,
2019), as BOIN is the name of the Bayesian optimal interval design.

The same check listed three URLs in README.md, of pages of the pkgdown site of
the package, as not found. The pages are new in this version and were published
when the site of version 2.3.0 was built, after the win-builder checks; they
now resolve. The R CMD check workflow, which fails on warnings, passed on all
five GitHub Actions platforms.

## Downstream dependencies

There are currently no downstream dependencies for this package.
