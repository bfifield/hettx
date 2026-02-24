## Test sweeping over range of plausibe taus for conservative p-value
FRTCI <- function(Y, Z, X = NULL, test.stat = SKS_stat, B=500,
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
            te.vec <- get_tau_vector( Y, Z, gamma=gamma, grid.size=grid.size, grid.gamma=grid.gamma )
        }else{
            te.vec <- get_tau_vector( Y, Z, X, gamma=gamma, grid.size=grid.size, grid.gamma=grid.gamma )
        }
    } else {
        grid.size <- length( te.vec )
    }
    te.hat <- attr( te.vec, "te.hat" )
    te.se <- attr( te.vec, "te.se" )
    te.MOE <- attr( te.vec, "te.MOE" )

    ## IMPUTE MISSING POTENTIAL OUTCOMES
    Y1.mat <- vapply(te.vec, function(te) ifelse(Z, Y, Y + te), numeric(length(Y)) )
    Y0.mat <- vapply(te.vec, function(te) ifelse(Z, Y - te, Y), numeric(length(Y)) )

    if( is.null(X) ){
        res <- generate_permutations( Y, Z=Z, test.stat=test.stat, Y0.mat=Y0.mat, Y1.mat=Y1.mat, B=B, n.cores, verbose=verbose, ... )
    }else{
        res <- generate_permutations( Y, Z=Z, test.stat=test.stat, Y0.mat=Y0.mat, Y1.mat=Y1.mat, B=B, n.cores, verbose=verbose, X=X, ... )
    }

    ci.p <- res$ci.p + gamma

    t <- res$ks.obs
    p.value <- max( ci.p )
    n <- length(Y)
    method <- "FRT CI Test for Treatment Effect Heterogeneity"
    DAT <- paste( n, " observations", sep="")

    if ( !return.matrix ) {
        ks.mat <- NULL
    } else {
        ks.mat <- res$ks.mat
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
FRTplug <- function( Y, Z, test.stat=SKS_stat, tau.hat=mean(Y[Z == 1]) - mean(Y[Z == 0]), ... ){
    mth <- FRTCI( Y, Z, test.stat=test.stat, te.vec=c(tau.hat), n.cores = 1, ...)
    mth$method <- "FRT Plug-in Test for Treatment Effect Heterogeneity"
    mth
}

## Test with treatment-covariate interactions to adjust for known heterogeneity
FRTCI_interact <- function( Y, Z, W, X=NULL, test.stat = SKS_stat_int_cov, B=500,
                            gamma=0.0001, grid.gamma=100*gamma,
                            grid.size=151, return.matrix=FALSE,
                            n.cores=1, verbose=TRUE, ... ) {

    grid.info <- get_testing_grid( Y, Z, W=W, X=X, gamma, grid.size )

    te.MOE <- NA

    te.grid <- grid.info$te.grid
    Y1.mat <- grid.info$Y1.mat
    Y0.mat <- grid.info$Y0.mat

    if( is.null(X) ){
        res <- generate_permutations( Y, Z, test.stat, Y0.mat, Y1.mat, B=B, n.cores=n.cores, verbose=verbose, W=W, ... )
    }else{
        res <- generate_permutations( Y, Z, test.stat, Y0.mat, Y1.mat, B=B, n.cores=n.cores, verbose=verbose, X=X, W=W, ... )
    }

    t <- res$ks.obs
    ci.p <- res$ci.p + gamma
    p.value <- max( ci.p )
    n <- length(Y)
    method <- "FRT CI Test for Tx Effect Heterogeneity Beyond a Systematic Model"
    DAT <- paste( n, " observations", sep="")

    if ( !return.matrix ) {
        ks.mat <- NULL
    } else {
        ks.mat <- res$ks.mat
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
