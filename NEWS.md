# hettx 1.0.1

* Fixed CRAN submission NOTEs (excluded `cran-comments.md` and
  `CRAN-SUBMISSION` from package build).

# hettx 1.0.0

## New features

* Added `tidy()` and `glance()` methods for all result objects
  (`FRTCI.test`, `RI.regression.result`, `RI.R2.result`) via the
  `generics` package, enabling tidyverse-style workflows.

* Added `test_stat_info()` with structured output and dynamic formatting,
  replacing the old hard-coded text display.

## Function renaming

* All user-facing functions have been renamed from dot-case to snake_case
  for consistency with modern R conventions:
    - `detect_idiosyncratic()` and `estimate_systematic()` (unchanged)
    - `get.p.value()` -> `get_p_value()`
    - `variance.ratio.test()` -> `variance_ratio_test()`
    - `make.randomized.dat()` -> `make_randomized_dat()`
    - `make.randomized.compliance.dat()` -> `make_randomized_compliance_dat()`
    - `test.stat.info()` -> `test_stat_info()`
    - All test statistic functions: `KS.stat` -> `KS_stat`,
      `SKS.stat` -> `SKS_stat`, `SKS.stat.cov` -> `SKS_stat_cov`, etc.
    - `make.linear.data()` -> `make_linear_data()`,
      `make.quadradic.data()` -> `make_quadradic_data()`,
      `make.skew.data()` -> `make_skew_data()`

* The old dot-case names continue to work but emit deprecation warnings
  via `.Deprecated()`. They will be removed in a future release.

## Dependency changes

* Removed `formula.tools` dependency; replaced with base R formula parsing
  (`all.vars()`, manual LHS/RHS extraction).

* Removed `plyr` dependency; replaced `plyr::ddply()`/`plyr::summarize()`
  with base R `aggregate()` in test statistic functions.

* Removed unused `dplyr`, `tidyr`, and `purrr` imports.

* Added `generics` to Imports (for `tidy()` and `glance()` methods).

## Code quality improvements

* Standardized assignment operators to `<-` throughout the codebase.

* Replaced `class(x) %in%` with `inherits()` for proper S3 class checking.

* Replaced `sapply()` with `vapply()` (explicit `FUN.VALUE`) across the
  codebase for type-safe apply operations.

* Renamed internal dot-case helper functions to snake_case.

* Added roxygen documentation to internal helper functions in
  `helper_frt.R`.

## Infrastructure

* Bumped minimum R version from 2.14.0 to 4.0.0.

* Removed unnecessary `stringsAsFactors = FALSE` arguments (default
  since R 4.0.0).

* Removed redundant manual `@usage` roxygen tags (17 instances).

* Updated `expect_is()` to `expect_s3_class()` in tests.

* Fixed CRAN NOTE: arXiv references now use DOI format
  (`doi:10.48550/arXiv.XXXX.XXXXX`).

* Added `importFrom(stats, terms)` to NAMESPACE.

* Added R-hub GitHub Actions workflow for cross-platform CRAN check
  testing.

* Added cross-branch equivalence tests verifying numeric results match
  the previous release.
