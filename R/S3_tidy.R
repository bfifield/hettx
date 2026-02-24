#' Tidy a FRTCI.test result
#'
#' @param x A \code{FRTCI.test} object from \code{detect_idiosyncratic()}.
#' @param ... Additional arguments (ignored).
#'
#' @return A data.frame with columns: statistic, p.value, p.value.plug,
#'   method, test.stat, estimate, std.error.
#'
#' @importFrom generics tidy
#' @export
tidy.FRTCI.test <- function(x, ...) {
    data.frame(
        statistic = x$statistic,
        p.value = x$p.value,
        p.value.plug = x$p.value.plug,
        method = x$method,
        test.stat = x$test.stat,
        estimate = x$te.hat,
        std.error = x$te.SE
    )
}

#' Glance at a FRTCI.test result
#'
#' @param x A \code{FRTCI.test} object from \code{detect_idiosyncratic()}.
#' @param ... Additional arguments (ignored).
#'
#' @return A one-row data.frame with model-level summary statistics.
#'
#' @importFrom generics glance
#' @export
glance.FRTCI.test <- function(x, ...) {
    data.frame(
        statistic = x$statistic,
        p.value = x$p.value,
        p.value.plug = x$p.value.plug,
        estimate = x$te.hat,
        std.error = x$te.SE,
        B = x$B,
        gamma = x$gamma,
        n = x$n
    )
}

#' Tidy an RI.regression.result
#'
#' @param x An \code{RI.regression.result} object from \code{estimate_systematic()}.
#' @param ... Additional arguments (ignored).
#'
#' @return A data.frame with one row per coefficient, containing columns:
#'   term, estimate, std.error.
#'
#' @importFrom generics tidy
#' @export
tidy.RI.regression.result <- function(x, ...) {
    co <- coef(x)
    se <- SE(x)
    data.frame(
        term = names(co),
        estimate = unname(co),
        std.error = unname(se)
    )
}

#' Glance at an RI.regression.result
#'
#' @param x An \code{RI.regression.result} object from \code{estimate_systematic()}.
#' @param ... Additional arguments (ignored).
#'
#' @return A one-row data.frame with model-level summary statistics.
#'
#' @importFrom generics glance
#' @export
glance.RI.regression.result <- function(x, ...) {
    data.frame(
        method = x$method,
        ATE = x$ATE,
        SE.ATE = x$SE.ATE,
        chisq.stat = as.numeric(x$chisq.stat),
        p.value = as.numeric(x$p.value),
        SD.Y0 = x$SD.Y0,
        SD.Y1 = x$SD.Y1
    )
}

#' Tidy an RI.R2.result
#'
#' @param x An \code{RI.R2.result} object from \code{R2()}.
#' @param ... Additional arguments (ignored).
#'
#' @return A data.frame with R-squared bound estimates. For ITT results,
#'   one row. For LATE results, three rows (Compliers, Noncompliers,
#'   Covariates and compliers).
#'
#' @importFrom generics tidy
#' @export
tidy.RI.R2.result <- function(x, ...) {
    s <- summary(x)
    r2 <- s$hte_r2
    data.frame(
        term = rownames(r2),
        R2.lower = r2[[1]],
        R2.lower.sharp = r2[[2]],
        R2.upper = r2[[3]],
        row.names = NULL
    )
}

#' Glance at an RI.R2.result
#'
#' @param x An \code{RI.R2.result} object from \code{R2()}.
#' @param ... Additional arguments (ignored).
#'
#' @return A one-row data.frame with model-level summary statistics.
#'
#' @importFrom generics glance
#' @export
glance.RI.R2.result <- function(x, ...) {
    s <- summary(x)
    if (s$method == "ITT") {
        data.frame(
            type = s$method,
            Sdd = s$hte_variance_systematic
        )
    } else {
        data.frame(
            type = s$method,
            Sdd = s$hte_variance_systematic_compliers,
            LATE = s$LATE,
            ITT = s$ITT,
            prop_compliers = s$prop_compliers
        )
    }
}
