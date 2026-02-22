#' Generate truncated multivariate normal draws
#'
#' Draws from a multivariate normal distribution, rejecting any samples
#' outside the confidence region defined by alpha.
#'
#' @param n Number of draws to generate.
#' @param mu Mean vector.
#' @param sigma Covariance matrix.
#' @param alpha Tail probability to exclude (defines the chi-squared confidence region).
#'
#' @return An n x K matrix of draws from the truncated distribution.
#'
#' @keywords internal
rcrmtvnorm <- function(n, mu, sigma, alpha)  {
    ## dimension
    K <- length(mu)

    ## return a matrix
    Mat <- matrix(0, n, K)

    ## critical value of chisq distribution of dof = K
    ctv.chisq <- qchisq(p = 1 - alpha, df = K)

    count <- 0
    repeat{

        draw.try <- rmvnorm(n = 1, mean = mu, sigma = sigma)
        draw.try <- as.vector(draw.try)

        ##within confidence region or not?
        chisq.try <- t(draw.try - mu)%*%solve(sigma, draw.try - mu)
        if(chisq.try <= ctv.chisq)
        {
            count <- count + 1
            Mat[count, ] <- draw.try

            if(count == n) break
        }## if

    }## repeat

    return(Mat)
}

#' Generate a grid of treatment effect values to test
#'
#' Creates a sequence of tau values spanning a confidence interval around
#' the estimated treatment effect, with oversampling near the point estimate.
#'
#' @param Y Outcome vector.
#' @param Z Treatment assignment vector (0/1).
#' @param X Optional covariate matrix for regression-adjusted estimation.
#' @param gamma Width of the confidence interval. Default is 0.001.
#' @param grid.size Number of grid points. Default is 21.
#' @param grid.gamma Controls oversampling density near the estimate. Default is 100*gamma.
#'
#' @return A numeric vector of tau values with attributes \code{te.hat},
#'   \code{te.se}, and \code{te.MOE}.
#'
#' @keywords internal
get_tau_vector <- function( Y, Z, X = NULL, gamma=0.001, grid.size=21, grid.gamma = 100*gamma ) {

    if ( is.null( X ) ) {

        te.hat <- mean(Y[Z == 1]) - mean(Y[Z == 0])
        te.se <- sqrt( var(Y[Z == 1])/sum(Z) + var(Y[Z == 0])/sum(1 - Z))

    } else {
        lm.tau <- lm( Y ~ Z + X )

        te.hat <- as.numeric(coef(lm.tau)["Z"])
        te.se <- as.numeric(sqrt( diag(vcov(lm.tau))["Z"]))
    }

    ## oversample points near the estimated tau-hat
    te.MOE <- qnorm(1 - gamma/2)*te.se
    te.vec <- te.hat + (te.MOE/qnorm(1-grid.gamma)) * qnorm( seq( grid.gamma, 1-grid.gamma, length.out=grid.size ) )
    te.vec

    attr( te.vec, "te.hat" ) <- te.hat
    attr( te.vec, "te.se" ) <- te.se
    attr( te.vec, "te.MOE" ) <- te.MOE

    te.vec
}

#' Build a grid of treatment effect models for interaction testing
#'
#' Estimates a linear model \code{Y ~ Z + W + Z:W (+X)} and samples from
#' the confidence region of the Z-related coefficients. For each sampled
#' model, computes individual imputed treatment effects and science tables
#' of imputed potential outcomes.
#'
#' @param Y Outcome vector.
#' @param Z Treatment assignment vector (0/1).
#' @param W Covariate matrix for treatment effect interactions.
#' @param X Optional covariate matrix for additional adjustment.
#' @param gamma Tail probability for the confidence region. Default is 0.0001.
#' @param grid.size Number of models to sample. Default is 150.
#'
#' @return A list with components:
#'   \describe{
#'     \item{te.grid}{grid.size x p matrix, each row a different effect model.}
#'     \item{te.mat}{N x grid.size matrix of individual treatment effects.}
#'     \item{Y0.mat}{N x grid.size matrix of imputed control potential outcomes.}
#'     \item{Y1.mat}{N x grid.size matrix of imputed treated potential outcomes.}
#'   }
#'
#' @keywords internal
get_testing_grid <- function( Y, Z, W, X=NULL, gamma=0.0001, grid.size=150 ) {

    ## get sample of treatment effect models to calculate p-values for
    if ( !is.null(X) ) {
        lm.tau <- lm(Y ~ Z + W + Z:W + X)
        ## if have duplicates with X and W, need to drop from consideration
        drp <- is.na( coef(lm.tau) )
        cof <- coef(lm.tau)[ !drp ]
    } else {
        lm.tau <- lm( Y ~ Z + W + Z:W )
        cof <- coef( lm.tau )
    }
    te.model <- grep( "Z", names(cof) )

    stopifnot( all( rownames( vcov(lm.tau) ) == names(cof) ) )

    te.hat <- as.numeric( cof[te.model] )
    te.cov <- vcov(lm.tau)[te.model,te.model]

    te.grid <- te.hat
    if ( grid.size > 1 ) {
        ## sample points from our confidence region, focusing on points
        ## close to our point estimate
        te.grid <-  rbind( te.grid, rcrmtvnorm(grid.size - 1, mu = te.hat, sigma = te.cov, alpha = gamma ) )
        colnames( te.grid ) <- colnames( te.cov )
    }

    ## calculate individual treatment effects for different treatment impact models
    ## each row is individual and each column is specific treatment effect model
    te.mat <- t(te.grid %*% t(cbind(1, W)))

    Y0.mat <- Y*(1 - Z) + (Y - te.mat) * Z
    Y1.mat <- Y*Z + (Y + te.mat) * (1 - Z)

    list( te.grid=te.grid, te.mat=te.mat, Y0.mat=Y0.mat, Y1.mat=Y1.mat )
}

#' Generate permutation distributions for FRT tests
#'
#' Core engine used by \code{FRTCI} and \code{FRTCI_interact} to compute
#' the randomization distribution of a test statistic under permuted
#' treatment assignments.
#'
#' @param Y Outcome vector.
#' @param Z Treatment assignment vector (0/1).
#' @param test.stat Test statistic function taking (Y, Z, ...).
#' @param Y0.mat N x grid.size matrix of imputed control potential outcomes.
#' @param Y1.mat N x grid.size matrix of imputed treated potential outcomes.
#' @param B Number of permutations.
#' @param n.cores Number of cores for parallel computation.
#' @param get.z.star Optional function to generate permuted Z vectors.
#'   If NULL, uses \code{sample(Z)}.
#' @param verbose Whether to display a progress bar. Default is TRUE.
#' @param ... Additional arguments passed to \code{test.stat}.
#'
#' @return A list with components:
#'   \describe{
#'     \item{ks.obs}{Observed test statistic.}
#'     \item{ks.mat}{B x grid.size matrix of permuted test statistics.}
#'     \item{ci.p}{Vector of p-values for each grid point.}
#'   }
#'
#' @keywords internal
generate_permutations <- function( Y, Z, test.stat, Y0.mat, Y1.mat, B, n.cores, get.z.star=NULL, verbose = TRUE, ... ) {

    ## SET UP STORAGE MATRICES
    n.te.vec <- ncol(Y0.mat)

    ## CALCULATE OBSERVED TEST STATISTICS
    ks.obs <- test.stat( Y, Z, ... )

    if ( verbose ) {
        pb <- txtProgressBar(min = 0, max = B, style = 3)
    }

    ## DECLARE PARALLEL ENVIRONMENT
    if(n.cores == 1){
        '%oper%' <- foreach::'%do%'
    }else {
        '%oper%' <- foreach::'%dopar%'
        cl <- makePSOCKcluster(n.cores, outfile="")
        registerDoParallel(cl)
        on.exit(stopCluster(cl))
    }

    ## RANDOMIZATION DISTRIBUTION
    ks.mat <- foreach(b = 1:B, .combine = "rbind") %oper% {
        if ( verbose ) {
            setTxtProgressBar(pb, b)
        }
        if ( !is.null( get.z.star ) ) {
            Z.star <- get.z.star( Z, ... )
        } else {
            Z.star <- sample(Z)
        }

        ## CI method
        ci.out <- vapply(seq_len(n.te.vec),
                         function(i){
                             ## CALCULATE RANDOMIZED VALUE
                             Yobs.star <- Z.star*Y1.mat[,i] + (1 - Z.star)*Y0.mat[,i]

                             ## COMPUTE TEST STATISTICS
                             ks.star <- test.stat(Yobs.star, Z.star, ... )

                             ## OUTPUT
                             ks.star
                         }, numeric(1))

        return(as.numeric(ci.out))
    }

    if ( verbose ) {
        close(pb)
    }

    ## CALCULATE P-VALUES
    ci.p <- apply(ks.mat, 2, function(ks.star){
        sum(ks.star >= ks.obs)/B
    })

    list( ks.obs=ks.obs, ks.mat=ks.mat, ci.p=ci.p )
}
