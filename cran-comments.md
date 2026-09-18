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

* local: Windows 11, R <<FILL IN: R.version.string>>
* GitHub Actions: windows-latest (release), macOS-latest (release),
  ubuntu-latest (release, devel, oldrel-1)
* win-builder: devel and release <<FILL IN: confirm after devtools::check_win_*>>

## R CMD check results

0 errors | 0 warnings | 0 notes

## Downstream dependencies

There are currently no reverse dependencies on CRAN.
