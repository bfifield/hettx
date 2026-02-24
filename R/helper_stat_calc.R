##
## This file has all the various test statistic functions provided in the package.
##



#' test_stat_info
#'
#' A list of test statistics for detect.idiosyncratic(), and information on use cases when each is appropriate.
#'
#' @examples
#' test_stat_info()
#'
#' @export
test_stat_info <- function() {

    stats <- list(
        "No covariate adjustment or interactions" = list(
            KS_stat = "Calculate classic (not shifted) KS statistic. If tau passed, Y1 will be shifted by tau.",
            SKS_stat = "Shifted KS statistic. Calculate KS distance between Y0 and Y1 shifted by tau.",
            rq_stat = "Classic KS statistic via quantile regression with covariates."
        ),
        "Covariate adjustment, no interactions" = list(
            SKS_stat_cov_pool = "Shifted KS statistic with covariates to increase precision.",
            SKS_stat_cov = paste("Shifted KS statistic with covariates with model for outcomes calculated on",
                                 "control group only. This avoids 'splitting' the treatment variation between",
                                 "treatment and control groups. We recommend this method over the 'pool' method",
                                 "in SKS_stat_cov_pool."),
            SKS_stat_cov_rq = "Shifted KS statistic via quantile regression with covariates.",
            rq_stat_cond_cov = paste("KS statistic via quantile regression with covariates. Conditional approach;",
                                     "see Koenker and Xiao (2002) for more information."),
            rq_stat_uncond_cov = paste("KS statistic via quantile regression with covariates. Unconditional",
                                       "approach; see Firpo (2007) for more information.")
        ),
        "Interactions, with or without covariate adjustment" = list(
            SKS_stat_int_cov_pool = paste("Shifted KS statistic with a linear treatment effect model and optional",
                                          "covariate adjustment. This will attempt to remove any systematic variation",
                                          "corresponding to the specified interaction model and then return an SKS",
                                          "statistic on the residuals to measure any variation 'left over'. This is",
                                          "the test statistic used in Ding, Feller, and Miratrix (2016), JRSS-B."),
            SKS_stat_int_cov = paste("Similar to SKS_stat_int_cov_pool, but this method first adjusts for baseline",
                                     "and then models treatment effects on the residuals to not split treatment",
                                     "effects. We recommend this method over the 'pool' method in",
                                     "SKS_stat_int_cov_pool.")
        ),
        "Interactions, no covariate adjustment" = list(
            WSKS_t = paste("Calculates the shifted KS statistic within each group specified in the interaction",
                           "model, and then aggregates together as a weighted average. Should be used when the",
                           "interaction model is a single categorical or factor covariate."),
            SKS_pool_t = paste("Subtract off group level treatment effect estimates and then look at KS statistic",
                               "on residuals. Should be used when the interaction model is a single categorical",
                               "or factor covariate.")
        )
    )

    width <- getOption("width", 80)
    name_width <- max(nchar(unlist(lapply(stats, names))))
    indent <- paste(rep(" ", name_width + 5), collapse = "")

    cat("\nAvailable test statistics for detect_idiosyncratic():\n")

    for (section in names(stats)) {
        cat("\n", section, ":\n", sep = "")
        entries <- stats[[section]]
        for (i in seq_along(entries)) {
            label <- sprintf("  %-*s - ", name_width, names(entries)[i])
            desc_lines <- strwrap(entries[[i]], width = width - nchar(indent),
                                  initial = "", prefix = "")
            cat(label, desc_lines[1], "\n", sep = "")
            if (length(desc_lines) > 1) {
                for (line in desc_lines[-1]) {
                    cat(indent, line, "\n", sep = "")
                }
            }
        }
    }
    cat("\n")

    invisible(stats)
}

#' KS_stat
#'
#' Calculate classic (not shifted) KS statistic; code is a modified version of R's ks.test().
#'
#' If tau passed, Y1 will be shifted by tau.
#'
#' @param Y Observed outcome vector
#' @param Z Treatment assigment vector
#' @param tau Value of treatment effect for shifting Y1. Default is NULL (Y1 not shifted).
#' @param alternative Direction of test ("two.sided", "less", "greater")
#'
#' @examples
#' df <- make_randomized_dat( 1000, gamma.vec=c(1,1,1,2), beta.vec=c(-1,-1,1,0) )
#' KS_stat(df$Yobs, df$Z)
#'
#' @return The value of the test.
#' @export
#'
#' @seealso detect_idiosyncratic
KS_stat <- function( Y, Z, tau = NULL, alternative = c("two.sided", "less", "greater") ) {
    x <- Y[Z==1]
    y <- Y[Z==0]
    if ( !is.null( tau ) ) {
        x <- x - tau
    }
    alternative <- match.arg(alternative)

    n.x <- length(x)
    n.y <- length(y)

    w <- c(x, y)
    z <- cumsum(ifelse(order(w) <= n.x, 1/n.x, -1/n.y))
    if (length(unique(w)) < (n.x + n.y)) {
        z <- z[c(which(diff(sort(w)) != 0), n.x + n.y)]
    }
    STATISTIC <- switch(alternative, two.sided = max(abs(z)),
                        greater = max(z), less = -min(z))

    STATISTIC
}

#' SKS_stat
#'
#' Shifted kolmogorov-smirnov statistic. Calculate KS distance between Y0 and Y1
#' shifted by sample tau.
#'
#' @param Y Observed outcome vector
#' @param Z Treatment assigment vector
#'
#' @examples
#' df <- make_randomized_dat( 1000, gamma.vec=c(1,1,1,2), beta.vec=c(-1,-1,1,0) )
#' SKS_stat(df$Yobs, df$Z)
#'
#' @return The value of the test.
#'
#' @seealso KS_stat, SKS_stat_cov
#' @seealso detect_idiosyncratic
#'
#' @export
SKS_stat <- function(Y, Z)
{
    Y1 <- Y[Z==1]
    Y0 <- Y[Z==0]

    Y1.star   <- Y1 - mean(Y1)
    Y0.star   <- Y0 - mean(Y0)

    unique.points <- c(Y1.star, Y0.star)

    Fn1 <- ecdf(Y1.star)
    Fn0 <- ecdf(Y0.star)

    difference <- Fn1(unique.points) - Fn0(unique.points)

    return(max(abs(difference)))

}

#' SKS_stat_cov_pool
#'
#' SKS_stat_cov_pool is the shifted kolmogorov-smirnov statistic with covariates
#' to increase precision.  This is the test statistic used Ding, Feller, and
#' Miratrix (2016), JRSS-B.
#'
#' @examples
#' df <- make_randomized_dat( 1000, gamma.vec=c(1,1,1,2), beta.vec=c(-1,-1,1,0) )
#' SKS_stat_cov_pool(df$Yobs, df$Z, df$A)
#'
#' @rdname SKS_stat_cov
#'
#' @export
SKS_stat_cov_pool <- function(Y, Z, X)
{
    this.lm <- lm(Y ~ Z + X)
    Y1.star <- this.lm$res[Z == 1]
    Y0.star <- this.lm$res[Z == 0]

    unique.points <- c(Y1.star, Y0.star)

    Fn1 <- ecdf(Y1.star)
    Fn0 <- ecdf(Y0.star)

    difference <- Fn1(unique.points) - Fn0(unique.points)

    return(max(abs(difference)))
}

#' SKS_stat_cov
#'
#' SKS_stat_cov is the shifted kolmogorov-smirnov statistic with covariates
#' with model for outcomes calculated on control group only.
#' This avoids "splitting" the treatment variation between tx
#' and co groups.
#' We recommend this method over the "pool" method.
#'
#' @param Y  Observed outcome vector
#' @param Z  Treatment assigment vector
#' @param X Additional pre-treatment covariates to adjust for in estimation, but not to interact with treatment.
#'
#' @examples
#' df <- make_randomized_dat( 1000, gamma.vec=c(1,1,1,2), beta.vec=c(-1,-1,1,0) )
#' SKS_stat_cov(df$Yobs, df$Z, df$A)
#'
#' @return The value of the test.
#'
#' @export
SKS_stat_cov <- function(Y, Z, X)
{
    this.lm <- lm(Y ~ X, subset = Z == 0)

    Y0.hat <- predict(this.lm, newdata = as.data.frame(X))
    Y0.res <- Y - Y0.hat

    Y1.star <- Y0.res[Z == 1] - mean(Y0.res[Z == 1])
    Y0.star <- Y0.res[Z == 0] - mean(Y0.res[Z == 0])

    unique.points <- c(Y1.star, Y0.star)

    Fn1 <- ecdf(Y1.star)
    Fn0 <- ecdf(Y0.star)

    difference <- Fn1(unique.points) - Fn0(unique.points)

    return(max(abs(difference)))
}

#' SKS_stat_int_cov_pool
#'
#' SKS_stat_int_cov_pool is a shifted kolmogorov-smirnov statistic with a linear
#' treatment effect model defined by W. It will attempt to remove any systematic
#' variation corresponding to W and then return a SKS statistic on the residuals
#' to measure any variation "left over".
#'
#'
#' @examples
#' df <- make_randomized_dat( 1000, gamma.vec=c(1,1,1,2), beta.vec=c(-1,-1,1,0) )
#' SKS_stat_int_cov_pool(Y = df$Yobs, Z = df$Z, W = df$A, X = df$B)
#'
#'
#' @rdname SKS_stat_int_cov
#'
#' @export
SKS_stat_int_cov_pool <- function( Y, Z, W, X=NULL )
{
    if ( !is.null( X ) ) {
        this.lm <- lm( Y ~ Z + X + W + Z:W )
    } else {
        this.lm <- lm( Y ~ Z + W + Z:W )
    }
    Y1.star <- this.lm$res[Z == 1]
    Y0.star <- this.lm$res[Z == 0]

    unique.points <- c(Y1.star, Y0.star)

    Fn1 <- ecdf(Y1.star)
    Fn0 <- ecdf(Y0.star)

    difference <- Fn1(unique.points) - Fn0(unique.points)

    return(max(abs(difference)))
}

#' SKS_stat_int_cov
#'
#' SKS_stat_int_cov() is a Shifted kolmogorov-smirnov statistic with a linear
#' treatment effect model defined by W. It will attempt to remove any systematic
#' variation corresponding to W and then return a SKS statistic on the residuals
#' to measure any variation "left over".
#'
#' X are _additional_ covariates to adjust for beyond those involved in
#' treatment effect model.  It will automatically ajust for W as well.  Do not
#' put a covariate in for both X and W.
#'
#' This is the test statistic used in Ding, Feller, and Miratrix (2016), JRSS-B.
#'
#' SKS_stat_int_cov first adjusts for baseline and then models treatment effect
#' on the residuals to not split treatment effects (see the vignette for more
#' information on this).
#'
#' We recommend SKS_stat_int_cov over the "pool" method.
#'
#' @inheritParams SKS_stat_cov
#' @param W Additional pre-treatment covariates to interact with T to define
#'   linear model of treatment effects.
#'
#' @examples
#' df <- make_randomized_dat( 1000, gamma.vec=c(1,1,1,2), beta.vec=c(-1,-1,1,0) )
#' SKS_stat_int_cov(Y = df$Yobs, Z = df$Z, W = df$A, X = df$B)
#'
#' @export
SKS_stat_int_cov <- function( Y, Z, W, X=NULL )
{
    ## First wipe out Y0 predicted by X via linear model
    if ( !is.null( X ) ) {
        this.lm <- lm(Y ~ X, subset = Z == 0)

        Y0.hat <- predict(this.lm, newdata = as.data.frame(X))
        Y <- Y - Y0.hat
    }

    ## now model treatment effect
    this.lm <- lm( Y ~ Z + W + Z:W )

    Y1.star <- this.lm$res[Z == 1]
    Y0.star <- this.lm$res[Z == 0]

    unique.points <- c(Y1.star, Y0.star)

    Fn1 <- ecdf(Y1.star)
    Fn0 <- ecdf(Y0.star)

    difference <- Fn1(unique.points) - Fn0(unique.points)

    return(max(abs(difference)))
}



##
## Other possible test statistics
##

#' SKS_stat_cov_rq
#'
#' Shifted kolmogorov-smirnov statistic with covariates and quantile regression.
#'
#' @inheritParams SKS_stat_cov
#'
#' @examples
#' df <- make_randomized_dat( 1000, gamma.vec=c(1,1,1,2), beta.vec=c(-1,-1,1,0) )
#' SKS_stat_cov_rq(df$Yobs, df$Z, df$A)
#'
#' @return The test statistic value.
#'
#' @export
SKS_stat_cov_rq <- function(Y, Z, X)
{

    this.lm <- lm(Y ~ Z + X)
    this.rq <- suppressWarnings(rq(Y ~ Z + X, tau = seq(0.05, 0.95, by = 0.05)))

    z.rq <- coef(this.rq)["Z",]
    return( max(z.rq - coef(this.lm)["Z"]) )

}

#' rq_stat
#'
#' rq_stat is the Kolmogorov-smirnov statistic via quantile regression with covariates without further adjustment.
#'
#' Warning: This function supresses all warnings of the `rq()` method call.
#'
#' @param Y  Observed outcome vector
#' @param Z  Treatment assigment vector
#' @param rq.pts Sequence of quantile points at which to evaluate the test. Default is seq(.1, .9, by = .1). Should not go beyond 0 and 1.
#'
#' @examples
#' df <- make_randomized_dat( 1000, gamma.vec=c(1,1,1,2), beta.vec=c(-1,-1,1,0) )
#' rq_stat(df$Yobs, df$Z)
#'
#' @return The value of the test.
#' @export
rq_stat <- function(Y, Z, rq.pts = seq(0.1, 0.9, by = 0.1))
{

    if(min(rq.pts) <= 0 | max(rq.pts) >= 1){
        stop("All values of rq.pts must be strictly greater than 0 and strictly less than 1.")
    }

    this.lm <- lm(Y ~ Z)
    this.rq <- suppressWarnings(rq(Y ~ Z, tau = rq.pts))

    z.rq <- coef(this.rq)["Z",]
    return( max(z.rq - coef(this.lm)["Z"]) )

}

#' rq_stat_cond_cov
#'
#' rq_stat_cond_cov does Kolmogorov-smirnov statistic via quantile regression
#' with covariates, with a conditional approach; see Koenker and Xiao (2002).
#'
#' Warning: This function supresses all warnings of the `rq()` method call.
#'
#' @param X Additional pre-treatment covariates to adjust for in estimation, but
#'   not to interact with treatment.
#'
#' @examples
#' df <- make_randomized_dat( 1000, gamma.vec=c(1,1,1,2), beta.vec=c(-1,-1,1,0) )
#' rq_stat_cond_cov(df$Yobs, df$Z, df$A)
#'
#' @rdname rq_stat
#' @export
rq_stat_cond_cov <- function(Y, Z, X, rq.pts = seq(0.1, 0.9, by = 0.1))
{

    if(min(rq.pts) <= 0 | max(rq.pts) >= 1){
        stop("All values of rq.pts must be strictly greater than 0 and strictly less than 1.")
    }

    this.lm <- lm(Y ~ Z + X)
    this.rq <- suppressWarnings(rq(Y ~ Z + X, tau = rq.pts))

    z.rq <- coef(this.rq)["Z",]
    return( max(z.rq - coef(this.lm)["Z"]) )

}


#' rq_stat_uncond_cov
#'
#' rq_stat_uncond_cov implements a Kolmogorov-smirnov statistic via quantile regression with covariates,
#' unconditional approach; see Firpo (2007).
#'
#' @examples
#' df <- make_randomized_dat( 1000, gamma.vec=c(1,1,1,2), beta.vec=c(-1,-1,1,0) )
#' rq_stat_uncond_cov(df$Yobs, df$Z, df$A)
#'
#' @rdname rq_stat
#' @export
rq_stat_uncond_cov <- function(Y, Z, X, rq.pts = seq(0.1, 0.9, by = 0.1))
{

    if(min(rq.pts) <= 0 | max(rq.pts) >= 1){
        stop("All values of rq.pts must be strictly greater than 0 and strictly less than 1.")
    }

    ## propensity score model
    this.glm <- glm(Z ~ X, family = binomial(link = logit))

    pscore <- predict(this.glm, type = "response")
    ipw.weights <- ifelse(Z, 1/pscore, 1/(1 - pscore))


    this.lm <- lm(Y ~ Z, weights = ipw.weights)
    this.rq <- suppressWarnings(rq(Y ~ Z, tau = rq.pts, weights = ipw.weights))

    z.rq <- coef(this.rq)["Z",]
    return( max(z.rq - coef(this.lm)["Z"]) )

}


##
## Test statistics for categorical covariates as discussed in the
## JRSS B Paper.
##

#' WSKS_t
#'
#' Weighted average of the group-level SKS statistics.  This is useful for a
#' blocked experiment.
#'
#' @param Y Observed outcome vector
#' @param Z Treatment assigment vector
#' @param W A a factor or categorical covariate.
#'
#' @examples
#' df <- make_randomized_dat( 1000, gamma.vec=c(1,1,1,2), beta.vec=c(-1,-1,1,0) )
#' df$W <- sample(c("A", "B", "C"), nrow(df), replace = TRUE)
#' WSKS_t(df$Yobs, df$Z, df$W)
#'
#' @return The value of the test.
#' @export
WSKS_t <- function( Y, Z, W ) {

    dat <- data.frame(Y=Y, Z=Z, W=W)
    dd <- do.call(rbind, lapply(split(dat, dat$W), function(d) {
        data.frame(t.sks = SKS_stat(d$Y, d$Z), n.k = nrow(d))
    }))
    n <- length(Y)
    return( sum( dd$t.sks * dd$n.k / n ) )
}

#' SKS_pool_t
#'
#' Subtract off group level treatment effect estimates and then look
#' at KS statistic on residuals.
#'
#' Distinct from the interacted lm in that the control units are not
#' shifted and centered with respect to eachother.
#'
#' @examples
#' df <- make_randomized_dat( 1000, gamma.vec=c(1,1,1,2), beta.vec=c(-1,-1,1,0) )
#' df$W <- sample(c("A", "B", "C"), nrow(df), replace = TRUE)
#' SKS_pool_t(df$Yobs, df$Z, df$W)
#'
#' @inheritParams WSKS_t
#'
#' @export
SKS_pool_t <- function( Y, Z, W ) {

    dat <- data.frame( Y=Y, Z=Z, W=W )
    mean1 <- tapply(dat$Y[dat$Z == 1], dat$W[dat$Z == 1], mean)
    mean0 <- tapply(dat$Y[dat$Z == 0], dat$W[dat$Z == 0], mean)
    taus <- mean1 - mean0
    return( KS_stat( dat$Y - dat$Z*taus[dat$W], dat$Z ) )
}
