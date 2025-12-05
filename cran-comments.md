## R CMD check results

0 errors | 0 warnings | 0 notes

## Test environments

* local: macOS (R 4.4.0)
* GitHub Actions (via usethis):
  - macOS-latest (release)
  - windows-latest (release)
  - ubuntu-latest (devel)
  - ubuntu-latest (release)
  - ubuntu-latest (oldrel-1)

## Submission notes

This is a new release (version 1.3.0) with the following updates:

* Added `verbose` parameter to control progress message output
* Enhanced documentation with standardized roxygen2 format
* Improved vignette output for R Markdown integration
* All tests pass successfully
* No breaking changes from previous version

## Reverse dependencies

There are currently no reverse dependencies for this package.

## Additional notes

* Package provides fast simulation tools for Bayesian Optimal Interval (BOIN) designs
* Results are validated against the BOIN package (within 0.5% accuracy)
* All examples run successfully without errors
* Vignettes build without warnings
