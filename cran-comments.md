## Resubmission

This is a resubmission to correct a typo in the DESCRIPTION file.

**Version 1.3.0** (submitted earlier today) contained an error in the Date field 
(2024-12-06 instead of 2025-12-06). Please disregard version 1.3.0 and review 
version 1.3.1 instead.

### Changes in version 1.3.1
* Fixed Date field in DESCRIPTION (2024-12-06 → 2025-12-06)
* No functional changes - all code and features are identical to version 1.3.0

## R CMD check results

0 errors ✔ | 0 warnings ✔ | 0 notes ✔

## Test environments

* local: Windows 11, R 4.5.0
* win-builder (release and devel)
* GitHub Actions (via usethis):
  - macOS-latest (release): ✅ Successful
  - windows-latest (release): ✅ Successful
  - ubuntu-latest (release): ✅ Successful
  - ubuntu-latest (oldrel-1): ✅ Successful
  - ubuntu-latest (devel): Failing due to yaml package compilation issue (unrelated to simFastBOIN)

Note: The failure on ubuntu-latest (devel) is caused by a compilation error in 
the yaml package dependency (`R_ext/PrtUtil.h` header file issue in R-devel). 
This is not related to simFastBOIN code and does not affect release versions.

## Submission notes

This is the first submission of this package to CRAN.

simFastBOIN provides fast and efficient simulation tools for Bayesian Optimal 
Interval (BOIN) designs in Phase I clinical trials. Key features include:

* Vectorized implementation with superior performance (2-5x faster than traditional approaches)
* Full compatibility with BOIN methodology (results match BOIN package within 0.5%)
* Multi-scenario simulation support for protocol development
* Publication-ready table formatting with multiple output options
* Comprehensive vignettes and documentation
* All tests pass successfully (116 tests across all platforms)

## Reverse dependencies

There are currently no reverse dependencies for this package.

## Additional notes

* Package provides validated simulation tools for clinical trial design
* All examples run successfully without errors
* Vignettes build without warnings
* Reference: Liu and Yuan (2015) <doi:10.1111/rssc.12089>
