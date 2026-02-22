# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

hettx is an R package implementing Fisherian and Neymanian methods for detecting and measuring treatment effect variation in randomized experiments. Based on Ding, Feller, and Miratrix (2016, JRSS-B and 2018, JASA).

## Build & Development Commands

```bash
# Generate roxygen documentation
R -e 'devtools::document()'

# Build package
R CMD BUILD hettx --resave-data

# Run CRAN checks
R CMD CHECK hettx_*.tar.gz --as-cran

# Install locally
R CMD INSTALL hettx_*.tar.gz

# Run all tests
R -e 'devtools::test()'

# Run a single test file
R -e 'testthat::test_file("tests/testthat/test_detect_idiosyncratic.R")'

# Full build/check cycle (uses builder.sh)
bash builder.sh
```

## Architecture

### Two main entry points

- **`detect_idiosyncratic()`** (`R/detect_idiosyncratic.R`) — Tests whether there is unexplained treatment effect variation via Fisherian Randomization Tests. Dispatches to three internal engines: `FRTplug()` (plug-in test), `FRTCI()` (confidence interval test), `FRTCI_interact()` (interaction test), all in `R/helper_frt.R` and `R/helper_FRTCI.R`.

- **`estimate_systematic()`** (`R/estimate_systematic.R`) — Estimates systematic treatment effect heterogeneity explained by covariates. Supports ITT, LATE, and Oracle modes. Core estimation logic lives in `R/helper_est_beta.R`.

### Helper modules

- `R/helper_stat_calc.R` — All test statistic functions (KS, SKS, rq-based variants, with/without covariates/interactions). These are passed by name or function reference to `detect_idiosyncratic()`.
- `R/helper_R2.R` — R-squared estimation for treatment effect variation.
- `R/helper_est_beta.R` — Effect estimation internals for ITT/LATE/Oracle (largest file, ~590 lines).
- `R/make_fake_data.R` — Data generation functions for simulations and testing.

### S3 object system

Results are returned as S3 objects (`FRTCI.test`, `RI.regression.result`, `RI.R2.result`) with print/summary/plot methods defined across `R/S3_print.R`, `R/S3_summary.R`, `R/S3_plot.R`.

### Key design patterns

- Formula-based interface (R formula objects) for specifying outcomes, treatment, covariates, and interactions.
- Parallelization via `foreach`/`doParallel` controlled by `n.cores` parameter in main functions.
- Monte Carlo permutation-based inference: `B` parameter controls number of permutations.

## Testing

Uses `testthat`. Three test files cover:
- `test_detect_idiosyncratic.R` — Hypothesis testing with various test statistics
- `test_estimate_systematic.R` — Effect estimation (Oracle, OLS, LATE)
- `test_data_generators.R` — Data generation validation

Tests involve randomization, so some use `set.seed()` for reproducibility. Many tests are slow due to Monte Carlo permutations.

## Datasets

- `ToyData` — Small simulated dataset (500 obs) used in examples/tests
- `Penn46_ascii` — Pennsylvania reemployment bonus experiment (6384 obs)

## Modernization Checklist

### High Priority
- [x] Replace `formula.tools` dependency with base R (`all.vars()`, manual LHS/RHS parsing)
- [x] Replace `plyr::ddply()`/`plyr::summarize()` with base R in `R/helper_stat_calc.R`
- [x] Audit `tidyr::gather()` import — removed (unused); also removed unused `dplyr` and `purrr` imports

### Medium Priority
- [x] Standardize assignment operators to `<-` throughout (currently mixes `=` and `<-`)
- [x] Replace `class(x) %in%` with `inherits()` in `detect_idiosyncratic.R` and `helper_est_beta.R`
- [x] Replace `sapply()` with `vapply()` (explicit `FUN.VALUE`) across codebase
- [x] Rename dot-case internal functions to snake_case (e.g., `get.tau.vector` -> `get_tau_vector`)

### Low Priority
- [ ] Bump minimum R version from 2.14.0 to 4.0.0+; drop `stringsAsFactors = FALSE`
- [ ] Remove manual `@usage` roxygen tags (17 instances)
- [ ] Replace `expect_is()` with `expect_s3_class()` in tests
- [ ] Consider `future`/`furrr` to replace `foreach`/`doParallel`
- [ ] Add roxygen docs to internal helpers in `helper_frt.R`
