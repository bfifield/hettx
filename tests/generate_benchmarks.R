#!/usr/bin/env Rscript
##
## Generate benchmark results from the current hettx installation.
##
## Usage:
##   1. Install the master branch:  R CMD INSTALL .
##   2. Run this script:            Rscript tests/generate_benchmarks.R
##   3. The file tests/master_benchmarks.rds will be created.
##
## These benchmarks are used by tests/testthat/test_cross_branch_equivalence.R
## to verify that code changes are behavioral no-ops.

library(hettx)

results <- list()

# ===== estimate_systematic vignette =====
set.seed(42)
df <- make_randomized_dat(10000, gamma.vec = c(1, 1, 1, 2), beta.vec = c(-1, -1, 1, 0))

# Basic RI estimation
rs <- estimate_systematic(Yobs ~ Z, interaction.formula = ~ A + B, data = df)
results$est_basic_coef <- coef(rs)
results$est_basic_vcov <- vcov(rs)
results$est_basic_SE <- SE(rs)
results$est_basic_confint <- confint(rs)

# OLS
M.ols.ours <- estimate_systematic(Yobs ~ Z, ~ A + B, data = df, method = "OLS")
results$est_ols_coef <- coef(M.ols.ours)
results$est_ols_beta_hat <- M.ols.ours$beta.hat

# OLS vs lm comparison
M0 <- lm(Yobs ~ (A + B) * Z, data = df)
results$est_ols_vs_lm_diff <- M.ols.ours$beta - coef(M0)[4:6]

# Control formula
rs_ctrl <- estimate_systematic(Yobs ~ Z, interaction.formula = ~ A + B,
                                control.formula = ~ C, data = df)
results$est_ctrl_coef <- coef(rs_ctrl)

# Both formulas
rsA2 <- estimate_systematic(Yobs ~ Z, ~ A + B, ~ A + B + C, data = df)
results$est_both_coef <- coef(rsA2)

# OLS with control
rsB <- estimate_systematic(Yobs ~ Z, ~ A + B, ~ C, data = df, method = "OLS")
results$est_ols_ctrl_coef <- coef(rsB)
rsB2 <- estimate_systematic(Yobs ~ Z, ~ A + B, ~ A + B + C, data = df, method = "OLS")
results$est_ols_both_coef <- coef(rsB2)

# Oracle
Moracle <- estimate_systematic(Y.1 + Y.0 ~ Z, ~ A + B, data = df)
results$est_oracle_coef <- coef(Moracle)
results$est_oracle_SE <- SE(Moracle)

# R2
set.seed(42)
df2 <- make_randomized_dat(1000, beta.vec = c(-1, 1, 1))
rs2 <- estimate_systematic(Yobs ~ Z, ~ A + B, data = df2, method = "OLS")
r2 <- R2(rs2)
results$r2_vals <- list(R2.lower = r2$R2.lower, R2.lower.sharp = r2$R2.lower.sharp,
                        R2.upper = r2$R2.upper, Sdd = r2$Sdd)

# R2 with idiosyncratic variation
set.seed(42)
df3 <- make_randomized_dat(1000, beta.vec = c(-1, 1, 1), ideo.sd = 3)
rs3 <- estimate_systematic(Yobs ~ Z, ~ A + B, data = df3, method = "OLS")
r2b <- R2(rs3)
results$r2b_vals <- list(R2.lower = r2b$R2.lower, R2.lower.sharp = r2b$R2.lower.sharp,
                         R2.upper = r2b$R2.upper, Sdd = r2b$Sdd)

# Non-compliance / LATE
set.seed(42)
beta_nc <- c(-1, 6, 0)
data_nc <- make_randomized_compliance_dat(10000, beta.vec = beta_nc)
rs_nc <- estimate_systematic(Yobs ~ D | Z, ~ A + B, data = data_nc)
results$late_coef <- coef(rs_nc)
results$late_beta_hat <- rs_nc$beta.hat
results$late_SE <- SE(rs_nc)

# R2 for LATE
r2_nc <- R2(rs_nc)
results$r2_late <- list(R2.lower = r2_nc$R2.lower, R2.lower.sharp = r2_nc$R2.lower.sharp,
                        R2.upper = r2_nc$R2.upper, Sdd = r2_nc$Sdd)

# 2SLS
rs2SLS <- estimate_systematic(Yobs ~ Z | D, ~ A + B, data = data_nc, method = "2SLS")
results$tsls_coef <- coef(rs2SLS)
results$tsls_beta_hat <- rs2SLS$beta.hat
results$tsls_SE <- SE(rs2SLS)

# ===== detect_idiosyncratic vignette =====
data(ToyData)

# Variance ratio test
results$vrt <- variance_ratio_test(ToyData$Y, ToyData$Z)

# Basic detect_idiosyncratic
set.seed(42)
tst1 <- detect_idiosyncratic(Y ~ Z, data = ToyData, B = 20, grid.size = 11, verbose = FALSE)
results$detect_basic <- get_p_value(tst1)

# With covariates
set.seed(42)
tst2 <- detect_idiosyncratic(Y ~ Z, data = ToyData,
                              control.formula = ~ x1 + x2 + x3 + x4,
                              B = 20, test.stat = "SKS_stat_cov", verbose = FALSE)
results$detect_cov <- get_p_value(tst2)

# With interaction formula
set.seed(42)
tst3 <- detect_idiosyncratic(Y ~ Z, data = ToyData,
                              interaction.formula = ~ x1 + x2,
                              control.formula = ~ x3 + x4,
                              B = 20, test.stat = "SKS_stat_int_cov", verbose = FALSE)
results$detect_int_cov <- get_p_value(tst3)

output_path <- file.path(dirname(sys.frame(1)$ofile), "master_benchmarks.rds")
cat("Saving benchmarks to:", output_path, "\n")
saveRDS(results, output_path)
cat("Done. Captured", length(results), "result sets.\n")
