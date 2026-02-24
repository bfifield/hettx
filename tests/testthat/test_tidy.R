library(testthat)
library(generics)
library(hettx)

context("tidy() and glance() methods")

test_that("tidy.FRTCI.test returns expected structure", {
    data(ToyData)
    tst <- detect_idiosyncratic(Y ~ Z, data = ToyData, B = 20,
                                grid.size = 11, verbose = FALSE)

    td <- tidy(tst)
    expect_s3_class(td, "data.frame")
    expect_true(all(c("statistic", "p.value", "p.value.plug",
                       "method", "test.stat", "estimate", "std.error") %in% names(td)))
    expect_equal(td$statistic, tst$statistic)
    expect_equal(td$p.value, tst$p.value)
    expect_equal(td$estimate, tst$te.hat)
})

test_that("glance.FRTCI.test returns one-row data.frame", {
    data(ToyData)
    tst <- detect_idiosyncratic(Y ~ Z, data = ToyData, B = 20,
                                grid.size = 11, verbose = FALSE)

    gl <- glance(tst)
    expect_s3_class(gl, "data.frame")
    expect_equal(nrow(gl), 1)
    expect_true(all(c("statistic", "p.value", "B", "gamma", "n") %in% names(gl)))
    expect_equal(gl$B, 20)
    expect_equal(gl$n, tst$n)
})

test_that("tidy.RI.regression.result returns one row per coefficient", {
    set.seed(42)
    df <- make_randomized_dat(100, beta.vec = c(-1, -1, 1))
    es <- estimate_systematic(Yobs ~ Z, data = df,
                              interaction.formula = ~ A + B)

    td <- tidy(es)
    expect_s3_class(td, "data.frame")
    expect_true(all(c("term", "estimate", "std.error") %in% names(td)))
    expect_equal(nrow(td), length(coef(es)))
    expect_equal(td$estimate, unname(coef(es)))
    expect_equal(td$std.error, unname(SE(es)))
})

test_that("glance.RI.regression.result returns one-row data.frame", {
    set.seed(42)
    df <- make_randomized_dat(100, beta.vec = c(-1, -1, 1))
    es <- estimate_systematic(Yobs ~ Z, data = df,
                              interaction.formula = ~ A + B)

    gl <- glance(es)
    expect_s3_class(gl, "data.frame")
    expect_equal(nrow(gl), 1)
    expect_true(all(c("method", "ATE", "SE.ATE", "chisq.stat",
                       "p.value", "SD.Y0", "SD.Y1") %in% names(gl)))
    expect_equal(gl$ATE, es$ATE)
    expect_equal(gl$chisq.stat, as.numeric(es$chisq.stat))
})

test_that("tidy.RI.R2.result works for ITT", {
    set.seed(42)
    df <- make_randomized_dat(100, beta.vec = c(-1, -1, 1))
    es <- estimate_systematic(Yobs ~ Z, data = df,
                              interaction.formula = ~ A + B)
    r2 <- R2(es)

    td <- tidy(r2)
    expect_s3_class(td, "data.frame")
    expect_true(all(c("term", "R2.lower", "R2.lower.sharp", "R2.upper") %in% names(td)))
    expect_equal(nrow(td), 1)
})

test_that("glance.RI.R2.result works for ITT", {
    set.seed(42)
    df <- make_randomized_dat(100, beta.vec = c(-1, -1, 1))
    es <- estimate_systematic(Yobs ~ Z, data = df,
                              interaction.formula = ~ A + B)
    r2 <- R2(es)

    gl <- glance(r2)
    expect_s3_class(gl, "data.frame")
    expect_equal(nrow(gl), 1)
    expect_equal(gl$type, "ITT")
    expect_true("Sdd" %in% names(gl))
})

test_that("tidy.RI.R2.result works for LATE", {
    set.seed(42)
    df <- make_randomized_compliance_dat(100)
    es <- estimate_systematic(Yobs ~ D | Z, data = df,
                              interaction.formula = ~ A + B)
    r2 <- R2(es)

    td <- tidy(r2)
    expect_s3_class(td, "data.frame")
    expect_equal(nrow(td), 3)
    expect_true(all(c("term", "R2.lower", "R2.lower.sharp", "R2.upper") %in% names(td)))
})

test_that("glance.RI.R2.result works for LATE", {
    set.seed(42)
    df <- make_randomized_compliance_dat(100)
    es <- estimate_systematic(Yobs ~ D | Z, data = df,
                              interaction.formula = ~ A + B)
    r2 <- R2(es)

    gl <- glance(r2)
    expect_s3_class(gl, "data.frame")
    expect_equal(nrow(gl), 1)
    expect_equal(gl$type, "LATE")
    expect_true(all(c("LATE", "ITT", "prop_compliers") %in% names(gl)))
})
