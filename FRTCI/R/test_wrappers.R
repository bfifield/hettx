## Set global variables for foreach
globalVariables('b')

#' fishpidetect
#'
#' Test for systematic treatment effect heterogeneity using
#' Fisherian permutation inference methods.
#' @usage fishpidetect(Y, Z, W, X, plugin, tau.hat, test.stat, te.vec,
#' B, gamma, grid.gamma, grid.size, return.matrix, n.cores, verbose, ...)
#'
#' @param Y  Observed outcome vector
#' @param Z  Treatment assigment vector 
#' @param W Additional pre-treatment covariates to interact with T to define linear model of treatment effects. Default is NULL.
#' @param X Additional pre-treatment covariates to adjust for in estimation, but not to interact with treatment. Default is NULL.
#' @param plugin Whether to calculate the plug-in p-value without sweeping over range of possible treatment effect magnitudes. Default is FALSE.
#' @param tau.hat The value of the plug-in treatment effect. Default is sample average treatment effect.
#' @param test.stat  Test statistic function to use on the data. Default is shifted Kolmogorov-Smirnov statistic.
#' @param te.vec Vector of taus to examine if you want to override generating ones automatically. Default is NULL.
#' @param B  Number of permutations to take. Default is 500.
#' @param gamma How wide of a CI to make around tau-hat for search. Default is 0.0001.
#' @param grid.gamma Parameter to govern where the grid points are sampled.  Bigger values means more samples towards the estimated tau-hat. Default is 100*gamma.
#' @param grid.size Number of points in the grid. Default is 151.
#' @param return.matrix  Whether to return the matrix of all the imputed statistics.  Default is FALSE.
#' @param n.cores Number of cores to use to parallelize permutation step. Default is 1.
#' @param verbose  Whether to print out progress bar when fitting and other diagnostics. Default is TRUE.
#' @param ... Extra arguments passed to the generate.permutations function and test.stat functions.
#'
#' @return If plug-in, the value of the test and the associated p-value. If not, a list with the value of
#' the test statistic on the observed data, the value of the CI-adjusted p-value, the plug-in p-value, and other information on the test.
#'
#' @examples
#' Z <- rep(c(0, 1), 100)
#' tau <- 4
#' Y <- ifelse(Z, rnorm(100, tau), rnorm(100, 0))
#' tst <- fishpidetect(Y, Z, B = 50, grid.size = 50)
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
fishpidetect <- function(Y, Z, W = NULL, X = NULL, plugin = FALSE, tau.hat = NULL,
                         test.stat = ifelse( is.null(W) & is.null(X), SKS.stat, ifelse( is.null(W), SKS.stat.cov, SKS.stat.int.cov ) ),
                         te.vec = NULL,
                         B = 500, gamma = 0.0001, grid.gamma = 100*gamma,
                         grid.size = 151, return.matrix = FALSE,
                         n.cores = 1, verbose = TRUE, ...){
    
    ## ------------------------
    ## Run checks on input data
    ## ------------------------
    ## NA and length checks on Y and Z
    if(sum(is.na(Y)) > 0){
        stop("You have NAs in your outcome vector.")
    }
    if(sum(is.na(Z)) > 0){
        stop("You have NAs in your treatment assignment vector.")
    }
    if(length(unique(Z)) > 2){
        stop("You have more than two unique values in your treatment assignment vector.")
    }
    if(length(Y) != length(Z)){
        stop("Y and Z must be the same length.")
    }

    ## Class checks on Y and Z
    if(!inherits(Y, "numeric")){
        stop("Y must be a numeric vector.")
    }
    if(inherits(Z, "logical")){
        Z <- ifelse(Z, 1, 0)
    }else if(inherits(Z, "factor") | inherits(Z, "character")){
        Z <- as.character(Z)
        tcat <- unique(Z)[1]
        Z <- ifelse(Z == tcat, 1, 0)
        cat("Z is factor or character class, setting", tcat, "as treatment category.\n")
    }

    ## Class and input checks on W and X
    if(plugin & !is.null(W)){
        warning("Running plug-in test but covariates provided for W, will not adjust for existing heterogeneity.")
    }
    if(plugin & is.null(tau.hat)){
        warning("Using plug-in test with no argument passed to tau.hat, will estimate tau.hat from the data.\n")
    }
    if(!is.null(te.vec) & !is.null(W)){
        warning("Adjusting for existing heterogeneity but treatment effect vector provided, ignoring provided treatment effect vector.")
    }
    if(!is.null(W)){
        if(!(identical(WSKS.t, test.stat) | identical(SKS.pool.t, test.stat))){
            if(!(inherits(W, "data.frame") | inherits(W, "matrix"))){
                stop("W must be either a data frame or a matrix.")
            }
        }else{
            if(!inherits(W, "factor")){
                stop("If using SKS.pool.t or WSKS.t as test statistics, W must be a factor vector.")
            }
        }
        if(!inherits(W, "factor")){
            na.check <- apply(W, 2, function(x){sum(is.na(x))})
        }else{
            na.check <- sum(is.na(W))
        }
        if(any(na.check > 0)){
            stop("You have NAs in your W matrix.")
        }
        if(!inherits(W, "factor")){
            unq.check <- apply(W, 2, function(x){length(unique(x))})
        }else{
            unq.check <- length(unique(W))
        }
        if(any(unq.check == 1)){
            stop("You have some columns in W with no variation.")
        }
        if(!inherits(W, "factor")){
            if(nrow(W) != length(Y)){
                stop("W must have as many observations as Y.")
            }
        }else{
            if(length(W) != length(Y)){
                stop("W must have as many observations as Y.")
            }
        }
    }
    if(!is.null(X)){
        if(!(inherits(X, "data.frame") | inherits(X, "matrix"))){
            stop("X must be either a data frame or a matrix.")
        }
        na.check <- apply(X, 2, function(x){sum(is.na(x))})
        if(any(na.check > 0)){
            stop("You have NAs in your X matrix.")
        }
        unq.check <- apply(X, 2, function(x){length(unique(x))})
        if(any(unq.check == 1)){
            stop("You have some columns in X with no variation.")
        }
        if(nrow(X) != length(Y)){
            stop("X must have as many observations as Y.")
        }
    }

    ## Checks on functions
    no.adj.funs <- c(KS.stat, SKS.stat, rq.stat)
    adj.funs <- c(SKS.stat.cov.pool, SKS.stat.cov, SKS.stat.cov.rq, rq.stat.cond.cov, rq.stat.uncond.cov)
    adj.int.funs <- c(SKS.stat.int.cov.pool, SKS.stat.int.cov)
    int.funs <- c(WSKS.t, SKS.pool.t)
    if(is.null(W) & is.null(X)){
        store.test <- rep(NA, length(no.adj.funs))
        for(i in 1:length(no.adj.funs)){
            store.test[i] <- identical(no.adj.funs[[i]], test.stat)
        }
        if(!any(store.test)){
            stop("You have provided an invalid test statistic when not adjusting for covariates or specifying interactions. Must provide one of KS.stat, SKS.stat, or rq.stat.")
        }
    }else if(is.null(W)){
        store.test <- rep(NA, length(adj.funs))
        for(i in 1:length(adj.funs)){
            store.test[i] <- identical(adj.funs[[i]], test.stat)
        }
        if(!any(store.test)){
            stop("You have provided an invalid test statistic when adjusting for covariates X but not specifying interactions. Must provide one of SKS.stat.cov.pool, SKS.stat.cov, SKS.stat.cov.rq, rq.stat.cond.cov, or rq.stat.uncond.cov.")
        }
    }else if(is.null(X)){
        just.int.funs <- c(adj.int.funs, int.funs)
        store.test <- rep(NA, length(just.int.funs))
        for(i in 1:length(just.int.funs)){
            store.test[i] <- identical(just.int.funs[[i]], test.stat)
        }
        if(!any(store.test)){
            stop("You have provided an invalid test statistic when specifying interactions W but not adjusting for covariates. Must provide one of SKS.stat.int.cov.pool, SKS.stat.int.cov, WSKS.t, or SKS.pool.t.")
        }
    }else{
        store.test <- rep(NA, length(adj.int.funs))
        for(i in 1:length(adj.int.funs)){
            store.test[i] <- identical(adj.int.funs[[i]], test.stat)
        }
        if(!any(store.test)){
            stop("You have provided an invalid test statistic when specifying interactions W but and adjusting for covariates X. Must provide one of SKS.stat.int.cov.pool or SKS.stat.int.cov.")
        }
    }
    
    ## ------------------------------------------------
    ## Detect whether to run plugin, FRTCI, or interact
    ## ------------------------------------------------
    if(plugin){ ## Plug-in test
        if(is.null(tau.hat)){
            tau.hat <- mean(Y[Z == 1]) - mean(Y[Z == 0])
        }
        fpi_out <- FRTplug(Y = Y, Z = Z, test.stat = test.stat, tau.hat = tau.hat, verbose= verbose, ...)
    }else if(is.null(W)){ ## FRTCI - with or without adjusting
        fpi_out <- FRTCI(Y = Y, Z = Z, X = X, test.stat = test.stat, B = B, gamma = gamma,
                         grid.gamma = grid.gamma, grid.size = grid.size,
                         te.vec = te.vec, return.matrix = return.matrix,
                         n.cores = n.cores, verbose = verbose, ...)
    }else{ ## FRTCI.interact
        fpi_out <- FRTCI.interact(Y = Y, Z = Z, W = W, X = X,
                                  test.stat = test.stat, B = B, gamma = gamma,
                                  grid.gamma = grid.gamma, grid.size = grid.size,
                                  return.matrix = return.matrix, n.cores = n.cores,
                                  verbose = verbose, ...)
    }

    return(fpi_out)

}

## Test sweeping over range of plausibe taus for conservative p-value
FRTCI <- function(Y, Z, X = NULL, test.stat = SKS.stat, B=500, 
                  gamma=0.0001, grid.gamma=100*gamma, 
                  grid.size=151,
                  te.vec=NULL, return.matrix=FALSE,
                  n.cores=1,
                  verbose=TRUE, ... ) {
    
    
    if ( is.null(te.vec) ) {
        if ( grid.size %% 2 == 0 ) {
            grid.size <- grid.size+1
        }

        if( is.null(X) ){
            te.vec <- get.tau.vector( Y, Z, gamma=gamma, grid.size=grid.size, grid.gamma=grid.gamma )
        }else{
            te.vec <- get.tau.vector( Y, Z, X, gamma=gamma, grid.size=grid.size, grid.gamma=grid.gamma )
        }
    } else {
        grid.size = length( te.vec )
    }
    te.hat <- attr( te.vec, "te.hat" )
    te.se <- attr( te.vec, "te.se" )
    te.MOE <- attr( te.vec, "te.MOE" )
    
    ## IMPUTE MISSING POTENTIAL OUTCOMES
    Y1.mat <- sapply(te.vec, function(te) ifelse(Z, Y, Y + te) )
    Y0.mat <- sapply(te.vec, function(te) ifelse(Z, Y - te, Y) )

    if( is.null(X) ){
        res <- generate.permutations( Y, Z=Z, test.stat=test.stat, Y0.mat=Y0.mat, Y1.mat=Y1.mat, B=B, n.cores, verbose=verbose, ... )
    }else{
        res <- generate.permutations( Y, Z=Z, test.stat=test.stat, Y0.mat=Y0.mat, Y1.mat=Y1.mat, B=B, n.cores, verbose=verbose, X=X, ... )
    }
    
    ci.p = res$ci.p + gamma
    
    t = res$ks.obs
    p.value = max( ci.p )
    n=length(Y)
    method = paste( "FRT CI Test for Tx Effect Heterogeneity with ", substitute(test.stat), sep="" )
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

## Test using plug-in sample average treatment effect
FRTplug <- function( Y, Z, test.stat=SKS.stat, tau.hat=mean(Y[Z == 1]) - mean(Y[Z == 0]), ... ){
    mth = FRTCI( Y, Z, test.stat=test.stat, te.vec=c(tau.hat), n.cores = 1, ...)
    mth$method = paste( "FRT Plug-in Test for Tx Effect Heterogeneity with ", substitute(test.stat), sep="" )
    mth
}

## Test with treatment-covariate interactions to adjust for known heterogeneity
FRTCI.interact <- function( Y, Z, W, X=NULL, test.stat = SKS.stat.int.cov, B=500, 
                            gamma=0.0001, grid.gamma=100*gamma, 
                            grid.size=151, return.matrix=FALSE, 
                            n.cores=1, verbose=TRUE, ... ) {
    
    grid.info = get.testing.grid( Y, Z, W=W, X=X, gamma, grid.size )
    
    te.MOE = NA
    
    te.grid = grid.info$te.grid
    Y1.mat <- grid.info$Y1.mat
    Y0.mat <- grid.info$Y0.mat

    if( is.null(X) ){
        res <- generate.permutations( Y, Z, test.stat, Y0.mat, Y1.mat, B=B, n.cores=n.cores, verbose=verbose, W=W, ... )
    }else{
        res <- generate.permutations( Y, Z, test.stat, Y0.mat, Y1.mat, B=B, n.cores=n.cores, verbose=verbose, X=X, W=W, ... )
    }
    
    t = res$ks.obs
    ci.p = res$ci.p + gamma
    p.value = max( ci.p )
    n=length(Y)
    method = paste( "FRT CI Test for Tx Effect Heterogeneity Beyond a Systematic Model with ", substitute(test.stat), sep="" )
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
              
              class = "FRTCI.test")
}

