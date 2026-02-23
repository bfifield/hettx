## Cross-branch equivalence tests
##
## These tests verify that numeric results from the modernize-codebase branch
## match the master branch exactly. Benchmarks were captured from master using
## the vignette computations with fixed seeds.
##
## To regenerate benchmarks, run:
##   Rscript tests/generate_benchmarks.R
##
## The benchmark file is: tests/master_benchmarks.rds

library(testthat)
library(hettx)

context("Cross-branch equivalence with master")

benchmark_path <- file.path(
  testthat::test_path(), "..", "master_benchmarks.rds"
)

skip_if_not(file.exists(benchmark_path),
            "Master branch benchmarks not found. Run tests/generate_benchmarks.R first.")

master <- readRDS(benchmark_path)

# ===== estimate_systematic: Basic estimation =====

test_that("Basic RI estimation matches master", {
  set.seed(42)
  df <- make_randomized_dat(10000, gamma.vec = c(1, 1, 1, 2), beta.vec = c(-1, -1, 1, 0))
  rs <- estimate_systematic(Yobs ~ Z, interaction.formula = ~ A + B, data = df)

  expect_equal(coef(rs), master$est_basic_coef, tolerance = 1e-10)
  expect_equal(vcov(rs), master$est_basic_vcov, tolerance = 1e-10)
  expect_equal(SE(rs), master$est_basic_SE, tolerance = 1e-10)
  expect_equal(confint(rs), master$est_basic_confint, tolerance = 1e-10)
})

# ===== estimate_systematic: OLS =====

test_that("OLS estimation matches master", {
  set.seed(42)
  df <- make_randomized_dat(10000, gamma.vec = c(1, 1, 1, 2), beta.vec = c(-1, -1, 1, 0))
  M.ols.ours <- estimate_systematic(Yobs ~ Z, ~ A + B, data = df, method = "OLS")

  expect_equal(coef(M.ols.ours), master$est_ols_coef, tolerance = 1e-10)
  expect_equal(M.ols.ours$beta.hat, master$est_ols_beta_hat, tolerance = 1e-10)

  # OLS vs lm comparison
  M0 <- lm(Yobs ~ (A + B) * Z, data = df)
  diff <- M.ols.ours$beta - coef(M0)[4:6]
  expect_equal(diff, master$est_ols_vs_lm_diff, tolerance = 1e-10)
})

# ===== estimate_systematic: Control formula =====

test_that("Control formula estimation matches master", {
  set.seed(42)
  df <- make_randomized_dat(10000, gamma.vec = c(1, 1, 1, 2), beta.vec = c(-1, -1, 1, 0))

  rs_ctrl <- estimate_systematic(Yobs ~ Z, interaction.formula = ~ A + B,
                                  control.formula = ~ C, data = df)
  expect_equal(coef(rs_ctrl), master$est_ctrl_coef, tolerance = 1e-10)

  rsA2 <- estimate_systematic(Yobs ~ Z, ~ A + B, ~ A + B + C, data = df)
  expect_equal(coef(rsA2), master$est_both_coef, tolerance = 1e-10)
})

# ===== estimate_systematic: OLS with control =====

test_that("OLS with control formula matches master", {
  set.seed(42)
  df <- make_randomized_dat(10000, gamma.vec = c(1, 1, 1, 2), beta.vec = c(-1, -1, 1, 0))

  rsB <- estimate_systematic(Yobs ~ Z, ~ A + B, ~ C, data = df, method = "OLS")
  expect_equal(coef(rsB), master$est_ols_ctrl_coef, tolerance = 1e-10)

  rsB2 <- estimate_systematic(Yobs ~ Z, ~ A + B, ~ A + B + C, data = df, method = "OLS")
  expect_equal(coef(rsB2), master$est_ols_both_coef, tolerance = 1e-10)
})

# ===== estimate_systematic: Oracle =====

test_that("Oracle estimation matches master", {
  set.seed(42)
  df <- make_randomized_dat(10000, gamma.vec = c(1, 1, 1, 2), beta.vec = c(-1, -1, 1, 0))

  Moracle <- estimate_systematic(Y.1 + Y.0 ~ Z, ~ A + B, data = df)
  expect_equal(coef(Moracle), master$est_oracle_coef, tolerance = 1e-10)
  expect_equal(SE(Moracle), master$est_oracle_SE, tolerance = 1e-10)
})

# ===== R2 =====

test_that("R2 estimation matches master", {
  set.seed(42)
  df2 <- make_randomized_dat(1000, beta.vec = c(-1, 1, 1))
  rs2 <- estimate_systematic(Yobs ~ Z, ~ A + B, data = df2, method = "OLS")
  r2 <- R2(rs2)

  expect_equal(r2$R2.lower, master$r2_vals$R2.lower, tolerance = 1e-10)
  expect_equal(r2$R2.lower.sharp, master$r2_vals$R2.lower.sharp, tolerance = 1e-10)
  expect_equal(r2$R2.upper, master$r2_vals$R2.upper, tolerance = 1e-10)
  expect_equal(r2$Sdd, master$r2_vals$Sdd, tolerance = 1e-10)
})

test_that("R2 with idiosyncratic variation matches master", {
  set.seed(42)
  df3 <- make_randomized_dat(1000, beta.vec = c(-1, 1, 1), ideo.sd = 3)
  rs3 <- estimate_systematic(Yobs ~ Z, ~ A + B, data = df3, method = "OLS")
  r2b <- R2(rs3)

  expect_equal(r2b$R2.lower, master$r2b_vals$R2.lower, tolerance = 1e-10)
  expect_equal(r2b$R2.lower.sharp, master$r2b_vals$R2.lower.sharp, tolerance = 1e-10)
  expect_equal(r2b$R2.upper, master$r2b_vals$R2.upper, tolerance = 1e-10)
  expect_equal(r2b$Sdd, master$r2b_vals$Sdd, tolerance = 1e-10)
})

# ===== Non-compliance / LATE =====

test_that("LATE estimation matches master", {
  set.seed(42)
  beta_nc <- c(-1, 6, 0)
  data_nc <- make_randomized_compliance_dat(10000, beta.vec = beta_nc)

  rs_nc <- estimate_systematic(Yobs ~ D | Z, ~ A + B, data = data_nc)
  expect_equal(coef(rs_nc), master$late_coef, tolerance = 1e-10)
  expect_equal(rs_nc$beta.hat, master$late_beta_hat, tolerance = 1e-10)
  expect_equal(SE(rs_nc), master$late_SE, tolerance = 1e-10)

  # R2 for LATE
  r2_nc <- R2(rs_nc)
  expect_equal(r2_nc$R2.lower, master$r2_late$R2.lower, tolerance = 1e-10)
  expect_equal(r2_nc$R2.lower.sharp, master$r2_late$R2.lower.sharp, tolerance = 1e-10)
  expect_equal(r2_nc$R2.upper, master$r2_late$R2.upper, tolerance = 1e-10)
  expect_equal(r2_nc$Sdd, master$r2_late$Sdd, tolerance = 1e-10)
})

test_that("2SLS estimation matches master", {
  set.seed(42)
  data_nc <- make_randomized_compliance_dat(10000, beta.vec = c(-1, 6, 0))

  rs2SLS <- estimate_systematic(Yobs ~ Z | D, ~ A + B, data = data_nc, method = "2SLS")
  expect_equal(coef(rs2SLS), master$tsls_coef, tolerance = 1e-10)
  expect_equal(rs2SLS$beta.hat, master$tsls_beta_hat, tolerance = 1e-10)
  expect_equal(SE(rs2SLS), master$tsls_SE, tolerance = 1e-10)
})

# ===== detect_idiosyncratic =====

test_that("Variance ratio test matches master", {
  data(ToyData)
  vrt <- variance_ratio_test(ToyData$Y, ToyData$Z)

  # Compare all numeric fields
  expect_equal(vrt$p.value, master$vrt$p.value, tolerance = 1e-10)
  expect_equal(vrt$F, master$vrt$F, tolerance = 1e-10)
  expect_equal(vrt$ratio, master$vrt$ratio, tolerance = 1e-10)
})

test_that("Basic detect_idiosyncratic matches master", {
  data(ToyData)
  set.seed(42)
  tst1 <- detect_idiosyncratic(Y ~ Z, data = ToyData, B = 20, grid.size = 11, verbose = FALSE)
  pv1 <- get_p_value(tst1)

  expect_equal(pv1, master$detect_basic, tolerance = 1e-10)
})

test_that("Covariate-adjusted detect_idiosyncratic matches master", {
  data(ToyData)
  set.seed(42)
  tst2 <- detect_idiosyncratic(Y ~ Z, data = ToyData,
                                control.formula = ~ x1 + x2 + x3 + x4,
                                B = 20, test.stat = "SKS_stat_cov", verbose = FALSE)
  pv2 <- get_p_value(tst2)

  expect_equal(pv2, master$detect_cov, tolerance = 1e-10)
})

test_that("Interaction detect_idiosyncratic matches master", {
  data(ToyData)
  set.seed(42)
  tst3 <- detect_idiosyncratic(Y ~ Z, data = ToyData,
                                interaction.formula = ~ x1 + x2,
                                control.formula = ~ x3 + x4,
                                B = 20, test.stat = "SKS_stat_int_cov", verbose = FALSE)
  pv3 <- get_p_value(tst3)

  expect_equal(pv3, master$detect_int_cov, tolerance = 1e-10)
})
