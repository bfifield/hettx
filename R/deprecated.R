## Deprecated function aliases
##
## These provide backward compatibility for the old dot-case function names.
## Each wrapper calls the new snake_case version and emits a deprecation warning.

#' Deprecated functions in hettx
#'
#' These functions have been renamed to use snake_case. The old dot-case
#' names still work but will emit a deprecation warning. Please update
#' your code to use the new names.
#'
#' @param ... Arguments passed to the replacement function.
#' @name hettx-deprecated
#' @keywords internal
NULL

# --- Test statistic functions (helper_stat_calc.R) ---

#' @rdname hettx-deprecated
#' @export
test.stat.info <- function() {
    .Deprecated("test_stat_info")
    test_stat_info()
}

#' @rdname hettx-deprecated
#' @export
KS.stat <- function(...) {
    .Deprecated("KS_stat")
    KS_stat(...)
}

#' @rdname hettx-deprecated
#' @export
SKS.stat <- function(...) {
    .Deprecated("SKS_stat")
    SKS_stat(...)
}

#' @rdname hettx-deprecated
#' @export
SKS.stat.cov.pool <- function(...) {
    .Deprecated("SKS_stat_cov_pool")
    SKS_stat_cov_pool(...)
}

#' @rdname hettx-deprecated
#' @export
SKS.stat.cov <- function(...) {
    .Deprecated("SKS_stat_cov")
    SKS_stat_cov(...)
}

#' @rdname hettx-deprecated
#' @export
SKS.stat.int.cov.pool <- function(...) {
    .Deprecated("SKS_stat_int_cov_pool")
    SKS_stat_int_cov_pool(...)
}

#' @rdname hettx-deprecated
#' @export
SKS.stat.int.cov <- function(...) {
    .Deprecated("SKS_stat_int_cov")
    SKS_stat_int_cov(...)
}

#' @rdname hettx-deprecated
#' @export
SKS.stat.cov.rq <- function(...) {
    .Deprecated("SKS_stat_cov_rq")
    SKS_stat_cov_rq(...)
}

#' @rdname hettx-deprecated
#' @export
rq.stat <- function(...) {
    .Deprecated("rq_stat")
    rq_stat(...)
}

#' @rdname hettx-deprecated
#' @export
rq.stat.cond.cov <- function(...) {
    .Deprecated("rq_stat_cond_cov")
    rq_stat_cond_cov(...)
}

#' @rdname hettx-deprecated
#' @export
rq.stat.uncond.cov <- function(...) {
    .Deprecated("rq_stat_uncond_cov")
    rq_stat_uncond_cov(...)
}

#' @rdname hettx-deprecated
#' @export
WSKS.t <- function(...) {
    .Deprecated("WSKS_t")
    WSKS_t(...)
}

#' @rdname hettx-deprecated
#' @export
SKS.pool.t <- function(...) {
    .Deprecated("SKS_pool_t")
    SKS_pool_t(...)
}

# --- Data generation functions (make_fake_data.R) ---

#' @rdname hettx-deprecated
#' @export
make.linear.data <- function(...) {
    .Deprecated("make_linear_data")
    make_linear_data(...)
}

#' @rdname hettx-deprecated
#' @export
make.quadradic.data <- function(...) {
    .Deprecated("make_quadradic_data")
    make_quadradic_data(...)
}

#' @rdname hettx-deprecated
#' @export
make.skew.data <- function(...) {
    .Deprecated("make_skew_data")
    make_skew_data(...)
}

#' @rdname hettx-deprecated
#' @export
make.randomized.dat <- function(...) {
    .Deprecated("make_randomized_dat")
    make_randomized_dat(...)
}

#' @rdname hettx-deprecated
#' @export
make.randomized.compliance.dat <- function(...) {
    .Deprecated("make_randomized_compliance_dat")
    make_randomized_compliance_dat(...)
}

# --- Detection/estimation functions ---

#' @rdname hettx-deprecated
#' @export
get.p.value <- function(...) {
    .Deprecated("get_p_value")
    get_p_value(...)
}

#' @rdname hettx-deprecated
#' @export
variance.ratio.test <- function(...) {
    .Deprecated("variance_ratio_test")
    variance_ratio_test(...)
}
