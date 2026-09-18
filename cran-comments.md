## Submission

This is version 2.0.0 of simFastBOIN, a major revision of the version currently
on CRAN, 1.3.2. Versions 1.4.0 and 1.5.0 were development versions on GitHub and
were never submitted.

The simulation engine has been rewritten in C++, so the package now has compiled
code and depends on Rcpp. Several defects that affected numerical results were
corrected, and the user-facing functions were renamed. The previous names remain
available and issue a deprecation warning, and NEWS.md documents every change and
the mapping between old and new names. The major version number reflects these
breaking changes.

Dependencies were reduced: Iso and kableExtra are no longer imported, and knitr
moved from Imports to Suggests.

The package reproduces the BOIN package trial by trial for a given seed, and this
is checked rather than asserted. The test suite compares against BOIN directly
where that package is installed, and inst/validation/compare-with-BOIN.R runs the
same comparison over 200 configurations covering every option combination,
including extrasafe and titration. All 200 agree to within floating point.

## Test environments

* local: Windows 11 x64 (build 26200), R 4.6.0 (2026-04-24 ucrt)
* GitHub Actions: macOS-latest (release), windows-latest (release),
  ubuntu-latest (release, devel, oldrel-1)
* win-builder: R-devel (2026-09-16 r90549 ucrt) and R 4.6.1
* macOS builder: R 4.6.1 Patched, aarch64-apple-darwin23

## R CMD check results

0 errors | 0 warnings | 1 note

The note is raised by the incoming feasibility check on win-builder, under both
R-devel and R-release:

    Possibly misspelled words in DESCRIPTION:
      comparator (20:29)

The word is spelled correctly. A comparator is the design or treatment that
another one is measured against, which is the sense used here: the traditional
3+3 design is included so that the BOIN results can be set beside it.

The local check, all five GitHub Actions platforms and the macOS builder were
clean, with no errors, warnings or notes.

## Downstream dependencies

There are currently no reverse dependencies on CRAN.
