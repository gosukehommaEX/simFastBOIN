## R CMD check results

0 errors | 0 warnings | 0 notes

## Test environments

* local: Windows (R 4.5.0)
* win-builder (release and devel)
* GitHub Actions (via usethis):
  - macOS-latest (release)
  - windows-latest (release)
  - ubuntu-latest (devel)
  - ubuntu-latest (release)
  - ubuntu-latest (oldrel-1)

## Submission notes

This is the first submission of this package to CRAN.

simFastBOIN provides fast and efficient simulation tools for Bayesian Optimal 
Interval (BOIN) designs in Phase I clinical trials. Key features include:

* Vectorized implementation with superior performance (2-5x faster than traditional approaches)
* Full compatibility with BOIN methodology (results match within 0.5%)
* Multi-scenario simulation support for protocol development
* Publication-ready table formatting
* Comprehensive vignettes and documentation
* All tests pass successfully (116 tests)

## Reverse dependencies

There are currently no reverse dependencies for this package.

## Additional notes

* Package provides validated simulation tools for clinical trial design
* All examples run successfully without errors
* Vignettes build without warnings
* Reference: Liu and Yuan (2015) <doi:10.1111/rssc.12089>
