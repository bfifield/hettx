## Deprecated function aliases
##
## These provide backward compatibility for the old dot-case function names.
## Each wrapper calls the new snake_case version and emits a deprecation warning.

# --- Test statistic functions (helper_stat_calc.R) ---

#' @rdname test_stat_info
#' @export
test.stat.info <- function() {
    .Deprecated("test_stat_info")
    test_stat_info()
}

#' @rdname KS_stat
#' @export
KS.stat <- function(...) {
    .Deprecated("KS_stat")
    KS_stat(...)
}

#' @rdname SKS_stat
#' @export
SKS.stat <- function(...) {
    .Deprecated("SKS_stat")
    SKS_stat(...)
}

#' @rdname SKS_stat_cov
#' @export
SKS.stat.cov.pool <- function(...) {
    .Deprecated("SKS_stat_cov_pool")
    SKS_stat_cov_pool(...)
}

#' @rdname SKS_stat_cov
#' @export
SKS.stat.cov <- function(...) {
    .Deprecated("SKS_stat_cov")
    SKS_stat_cov(...)
}

#' @rdname SKS_stat_int_cov
#' @export
SKS.stat.int.cov.pool <- function(...) {
    .Deprecated("SKS_stat_int_cov_pool")
    SKS_stat_int_cov_pool(...)
}

#' @rdname SKS_stat_int_cov
#' @export
SKS.stat.int.cov <- function(...) {
    .Deprecated("SKS_stat_int_cov")
    SKS_stat_int_cov(...)
}

#' @rdname SKS_stat_cov_rq
#' @export
SKS.stat.cov.rq <- function(...) {
    .Deprecated("SKS_stat_cov_rq")
    SKS_stat_cov_rq(...)
}

#' @rdname rq_stat
#' @export
rq.stat <- function(...) {
    .Deprecated("rq_stat")
    rq_stat(...)
}

#' @rdname rq_stat
#' @export
rq.stat.cond.cov <- function(...) {
    .Deprecated("rq_stat_cond_cov")
    rq_stat_cond_cov(...)
}

#' @rdname rq_stat
#' @export
rq.stat.uncond.cov <- function(...) {
    .Deprecated("rq_stat_uncond_cov")
    rq_stat_uncond_cov(...)
}

#' @rdname WSKS_t
#' @export
WSKS.t <- function(...) {
    .Deprecated("WSKS_t")
    WSKS_t(...)
}

#' @rdname SKS_pool_t
#' @export
SKS.pool.t <- function(...) {
    .Deprecated("SKS_pool_t")
    SKS_pool_t(...)
}

# --- Data generation functions (make_fake_data.R) ---

#' @rdname make_linear_data
#' @export
make.linear.data <- function(...) {
    .Deprecated("make_linear_data")
    make_linear_data(...)
}

#' @rdname make_linear_data
#' @export
make.quadradic.data <- function(...) {
    .Deprecated("make_quadradic_data")
    make_quadradic_data(...)
}

#' @rdname make_linear_data
#' @export
make.skew.data <- function(...) {
    .Deprecated("make_skew_data")
    make_skew_data(...)
}

#' @rdname make_randomized_dat
#' @export
make.randomized.dat <- function(...) {
    .Deprecated("make_randomized_dat")
    make_randomized_dat(...)
}

#' @rdname make_randomized_dat
#' @export
make.randomized.compliance.dat <- function(...) {
    .Deprecated("make_randomized_compliance_dat")
    make_randomized_compliance_dat(...)
}

# --- Detection/estimation functions ---

#' @rdname get_p_value
#' @export
get.p.value <- function(...) {
    .Deprecated("get_p_value")
    get_p_value(...)
}

#' @rdname variance_ratio_test
#' @export
variance.ratio.test <- function(...) {
    .Deprecated("variance_ratio_test")
    variance_ratio_test(...)
}
