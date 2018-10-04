## Set global variables for foreach
globalVariables('b') 

#' FRTCI
#'
#' Conduct the FRT CI Method on passed data
#' @usage FRTCI(Y, Z, test.stat, B, gamma, grid.gamma, grid.size,
#' te.vec, return.matrix, n.cores, verbose, ...)
#' 
#' @param Y  Observed outcome vector
#' @param Z  Treatment assigment vector (1=Tx, 0=Co)
#' @param test.stat  Test statistic function to use on the data. Default is shifted Kolmogorov-Smirnov statistic.
#' @param B  Number of permutations to take. Default is 500.
#' @param gamma How wide of a CI to make around tau-hat for search. Default is 0.0001.
#' @param grid.gamma Parameter to govern where the grid points are sampled.  Bigger values means more samples towards the estimated tau-hat. Default is 100*gamma.
#' @param grid.size Number of points in the grid. Default is 151.
#' @param te.vec Vector of taus to examine if you want to override generating ones automatically. Default is NULL.
#' @param return.matrix  Whether to return the matrix of all the imputed statistics.  Default is FALSE.
#' @param n.cores Number of cores to use to parallelize permutation step. Default is 1.
#' @param verbose  Whether to print out progress bar when fitting and other diagnostics. Default is TRUE.
#' @param ... Extra arguments passed to the generate.permutations function (and the test.stat function).
#'
#' @return A list with the value of the test statistic on the observed data, the value of the CI-adjusted p-value,
#' the plug-in p-value, and other information on the test.
#' @examples
#' Z <- rep(c(0, 1), 100)
#' tau <- 4
#' Y <- ifelse(Z, rnorm(100, tau), rnorm(100, 0))
#' tst <- FRTCI(Y, Z, B = 50, grid.size = 50)
#' @export
#' @importFrom graphics abline lines plot rug
#' @importFrom stats binom.test binomial coef confint ecdf glm lm lowess predict qchisq qnorm var vcov
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom quantreg rq
#' @importFrom plyr ddply summarize .
#' @importFrom mvtnorm rmvnorm
#' @importFrom foreach "%do%" "%dopar%" foreach
#' @importFrom parallel makePSOCKcluster stopCluster
#' @importFrom doParallel registerDoParallel
FRTCI <- function(Y, Z, test.stat = SKS.stat, B=500, 
                  gamma=0.0001, grid.gamma=100*gamma, 
                  grid.size=151,
                  te.vec=NULL, return.matrix=FALSE,
                  n.cores=1,
                  verbose=TRUE, ... ) {
    
    
    if ( is.null(te.vec) ) {
        if ( grid.size %% 2 == 0 ) {
            grid.size <- grid.size+1
        }
        
        te.vec <- get.tau.vector( Y, Z, gamma=gamma, grid.size=grid.size, grid.gamma=grid.gamma )
    } else {
        grid.size = length( te.vec )
    }
    te.hat <- attr( te.vec, "te.hat" )
    te.se <- attr( te.vec, "te.se" )
    te.MOE <- attr( te.vec, "te.MOE" )
    
    ## IMPUTE MISSING POTENTIAL OUTCOMES
    Y1.mat <- sapply(te.vec, function(te) ifelse(Z, Y, Y + te) )
    Y0.mat <- sapply(te.vec, function(te) ifelse(Z, Y - te, Y) )
    
    res <- generate.permutations( Y, Z=Z, test.stat=test.stat, Y0.mat=Y0.mat, Y1.mat=Y1.mat, B=B, n.cores, verbose=verbose, ... )
    
    ci.p = res$ci.p + gamma
    
    t = res$ks.obs
    p.value = max( ci.p )
    n=length(Y)
    method = paste( "FRT CI Test for Tx Effect Heterogeniety with ", substitute(test.stat), sep="" )
    DAT = paste( n, " observations", sep="")
    
    if ( !return.matrix ) {
        ks.mat=NULL
    } else {
        ks.mat = res$ks.mat
    }
    
    structure(list(statistic = t, p.value = p.value,
                   p.value.plug = ci.p[(grid.size+1)/2],
                   method=method,
                   data.name = DAT,
                   Y=Y, Z=Z, n=n, ci.p=ci.p, te.vec=te.vec,
                   te.hat=te.hat,
                   te.SE=te.se,
                   te.MOE = te.MOE,
                   B=B, gamma=gamma, ks.mat=ks.mat ),
              
              class = "FRTCI.test")
} 

#' FRTplug
#'
#' Return the plug-in p-value for plugging in tau-hat under the permutation method
#' Just calls FRTCI.test with the specific, single tau of tau-hat.
#'
#' @usage FRTplug(Y, Z, test.stat, tau.hat, ...)
#' 
#' @param Y Observed outcome vector
#' @param Z Treatment assigment vector (1=Tx, 0=Co)
#' @param test.stat Test statistic function to use on the data.
#' @param tau.hat The value of the plug-in treatment effect. Defaults to the sample average treatment effect.
#' @param ... Extra arguments passed to the generate.permutations function (and the test.stat function).
#'
#' @return A list with the value of the test statistic on the observed data, the value of the CI-adjusted p-value,
#' the plug-in p-value, and other information on the test.
#' @examples
#' Z <- rep(c(0, 1), 100)
#' tau <- 4
#' Y <- ifelse(Z, rnorm(100, tau), rnorm(100, 0))
#' tst <- FRTplug(Y, Z)
#' @export
FRTplug <- function( Y, Z, test.stat=SKS.stat, tau.hat=mean(Y[Z == 1]) - mean(Y[Z == 0]), ... ){
    mth = FRTCI( Y, Z, test.stat, te.vec=c(tau.hat), n.cores = 1, ...)
    mth$method = paste( "FRT Plug-in Test for Tx Effect Heterogeniety with ", substitute(test.stat), sep="" )
    mth
}

#' FRTCI.interact
#' 
#' Conduct the FRT CI Method, adjusting for covariates using a linear model, on passed data
#' for a model of effects W'beta with unknown beta.
#'
#' @usage FRTCI.interact(Y, Z, W, X, test.stat, B, gamma, grid.gamma, grid.size, return.matrix,
#' n.cores, verbose, ...)
#' 
#' @param Y  Observed outcome vector
#' @param Z  Treatment assigment vector (1=Tx, 0=Co)
#' @param W Additional pre-treatment covariates to interact with T to increase precision of treatment effect estimate.
#' @param X Additional pre-treatment covariates to adjust for in estimation, but not to interact with treatment.
#' @param test.stat  Test statistic function to use on the data. Default is shifted Kolmogorov-Smirnov statistic.
#' @param B  Number of permutations to take. Default is 500.
#' @param gamma How wide of a CI to make around tau-hat for search. Default is 0.0001.
#' @param grid.gamma Parameter to govern where the grid points are sampled.  Bigger values means more samples towards the estimated tau-hat. Default is 100*gamma.
#' @param grid.size Number of points in the grid. Default is 151.
#' @param return.matrix  Whether to return the matrix of all the imputed statistics.  Default is FALSE.
#' @param n.cores Number of cores to use to parallelize permutation step. Default is 1.
#' @param verbose  Whether to print out progress bar when fitting and other diagnostics. Default is TRUE.
#' @param ... Extra arguments passed to the generate.permutations function (and the test.stat function).
#'
#' @return A list with the value of the test statistic on the observed data, the value of the CI-adjusted p-value,
#' the plug-in p-value, and other information on the test.
#' @export
FRTCI.interact <- function( Y, Z, W, X=NULL, test.stat = SKS.stat.int.cov, B=500, 
                            gamma=0.0001, grid.gamma=100*gamma, 
                            grid.size=151, return.matrix=FALSE, 
                            n.cores=1, verbose=TRUE, ... ) {
    
    grid.info = get.testing.grid( Y, Z, W=W, X=X, gamma, grid.size )
    
    te.MOE = NA
    
    te.grid = grid.info$te.grid
    Y1.mat <- grid.info$Y1.mat
    Y0.mat <- grid.info$Y0.mat
    
    res <- generate.permutations( Y, Z, test.stat, Y0.mat, Y1.mat, B=B, n.cores=n.cores, verbose=verbose, X=X, W=W, ... )
    
    t = res$ks.obs
    ci.p = res$ci.p + gamma
    p.value = max( ci.p )
    n=length(Y)
    method = paste( "FRT CI Test for Tx Effect Heterogeniety Beyond a Systematic Model with ", substitute(test.stat), sep="" )
    DAT = paste( n, " observations", sep="")
    
    if ( !return.matrix ) {
        ks.mat=NULL
    } else {
        ks.mat = res$ks.mat
    }
    
    structure(list(statistic = t, p.value = p.value,
                   p.value.plug = ci.p[1],
                   method=method,
                   data.name = DAT,
                   Y=Y, Z=Z, n=n, ci.p=ci.p, te.grid=te.grid,
                   B=B, gamma=gamma, ks.mat=ks.mat,
                   W=W, X=X ),
              
              class = "FRTCI.interact.test")
}

