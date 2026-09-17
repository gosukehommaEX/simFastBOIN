## Submission

This is a major revision of simFastBOIN, currently on CRAN at version 1.3.2.

The simulation engine has been rewritten in C++, so the package now has compiled
code and depends on Rcpp. Several defects that affected numerical results were
corrected, and the user-facing functions were renamed. The previous names remain
available and issue a deprecation warning, and NEWS.md documents every change and
the mapping between old and new names.

Dependencies were reduced: Iso and kableExtra are no longer imported, and knitr
moved from Imports to Suggests.

## Test environments

* (to be completed before submission)

## R CMD check results

* (to be completed before submission)

## Downstream dependencies

There are currently no reverse dependencies on CRAN.
