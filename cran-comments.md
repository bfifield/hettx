## R CMD check results

0 errors | 0 warnings | 0 notes

## Test environments

* local macOS Sonoma 14.7.1, R 4.4.2
* R-hub: Ubuntu Linux (R-devel)
* R-hub: Windows (R-devel)

## Downstream dependencies

There are no reverse dependencies for this package on CRAN.

## Notes

This is a major version update (0.1.3 -> 1.0.0) that modernizes the
codebase. All user-facing dot-case functions have been renamed to
snake_case with backward-compatible deprecated aliases. Dependencies
on formula.tools and plyr have been removed. New tidy() and glance()
methods have been added. All changes are verified to be behavioral
no-ops via cross-branch equivalence tests.
