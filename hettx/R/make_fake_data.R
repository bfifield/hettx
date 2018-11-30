##
## This script provides different DGP for making fake data
## for simulation studies.
##
## Miratrix, 2014

## library( plyr )
## library( MASS )




##########################################################################################
##
## Methods to make random science tables
##
##########################################################################################


#' Generate dataset according to a linear model.
#'
#' Given the parameters, generate a dataset and return a potential outcomes schedule (science table)
#' of synthetic potential outcomes.
#'
#' @param n Sample size
#' @param gamma.vec   Control outcome surface
#' @param gamma2.vec  Quadratic terms
#' @param beta.vec   Treatment effect surface
#' @param ideo.sd   Ideosyncratic residual variation
#' @param quad.tx  Quadratic treatment effects?
#' @param stand.norm.X   Xs shifted around 0 or no
#' @param corr.X  TRUE or FALSE.  Have Xs correlated or no.
#'
#' @family data_generators
#'
#' @return List of elements of data (not data frame)
#' @export
#' @importFrom MASS mvrnorm
#' @importFrom stats cov median model.matrix pchisq pnorm quantile rexp rnorm sd
make.linear.data = function( n, gamma.vec = c( 1, 2, 2, 1 ),
                             gamma2.vec = NULL,
                             beta.vec = c(-1,-1,1),
                             ideo.sd = 0,
                             quad.tx=FALSE,
                             stand.norm.X=FALSE,
                             corr.X=!stand.norm.X ) {

    J = length( gamma.vec )
    if ( !is.null( gamma2.vec ) ) {
        J2 = length( gamma2.vec )
    } else {
        J2 = 0
    }
    K = length( beta.vec )
    NC = max( J, K, J2 ) - 1
    quad = !is.null(gamma2.vec) || quad.tx

    if ( corr.X ) {
        sig = diag( NC ) * 0.5 + matrix( rep( 0.5, NC^2 ), nrow=NC )
    } else {
        sig = diag( rep(1,NC) )
    }

    if ( stand.norm.X ) {
        mu=rep(0, NC)
    } else {
        mu=rep(3, NC)
    }

    X = mvrnorm( n, mu = mu, Sigma = sig )
    if ( quad ) {
            X2 = ( sweep( X, 2, STATS=mu, FUN="-" ) )^2
    }

    colnames(X) = LETTERS[1:ncol(X)]
    X = cbind( 1, X )
    Y.0 = as.numeric( X[,1:J] %*% gamma.vec  + rnorm( n ) )
    if ( !is.null( gamma2.vec ) ) {
        Y.0 = Y.0 + as.numeric( X2[,1:J2] %*% gamma2.vec )
    }
    if ( quad.tx ) {
         tau =  as.numeric( X[,1:K] %*% beta.vec + X2 %*% beta.vec )
    } else {
        tau = as.numeric( X[,1:K] %*% beta.vec )
    }
    if ( ideo.sd > 0 ) {
        tau = tau + rnorm( n, mean=0, sd=ideo.sd )
    }
    Y.1 = Y.0 + tau

    #data.frame( Yobs= Yobs, Z=Z, X=X )
    list( X=X[,-1], Y.0=Y.0, Y.1=Y.1, tau=tau )

}



#' @describeIn make.linear.data Generate dataset according to a quadratic model
#' @export
make.quadradic.data = function( n, beta.vec = c(-1,-1,1) ) {
    X = mvrnorm( n, mu=c(1,2), Sigma = matrix( c(1,0.5,0.5,1), nrow=2 ) )
    colnames(X) = LETTERS[1:ncol(X)]
    Y.0 = 1 + 2 * X[,1]^2 + 3*sqrt(abs(X[,2])) + rnorm( n )
    tau = beta.vec[1] + beta.vec[2] * X[,1] + beta.vec[3] * X[,2]

    Y.1 = Y.0 + tau

    list( X=X, Y.0=Y.0, Y.1=Y.1, tau=tau )

}



#' @describeIn make.linear.data Generate dataset with a skew
#' @export
make.skew.data = function( n, beta.vec = c(-1,-1,1) ) {
    X = mvrnorm( n, mu=c(1,2), Sigma = matrix( c(1,0.5,0.5,1), nrow=2 ) )
    X[,1] = X[,1] + rexp( n )
    colnames(X) = LETTERS[1:ncol(X)]
    Y.0 = 1 + 2 * X[,1]^2 + 3*sqrt(abs(X[,2])) + rnorm( n )
    tau = beta.vec[1] + beta.vec[2] * X[,1] + beta.vec[3] * X[,2]

    Y.1 = Y.0 + tau

    Y.1 = Y.1 + rexp( n ) - 1

    list( X=X, Y.0=Y.0, Y.1=Y.1, tau=tau )
}



##########################################################################################
##
## Generate random experiment data
##
##########################################################################################


#' Make fake data for simulations
#'
#' Randomize a science table to get observed outcomes and treatment assignment
#'
#' @param n  Sample size
#' @param p Proportion treated
#' @param science.table.generator   Data generator
#' @param include.POs  TRUE/FALSE. Keep POs
#' @param as.data.frame  TRUE/FALSE.  Return as dataframe or as list of elements.
#' @param ...  Additional to be passed to science.table.generator
#'
#' @return Either a list of elements or a dataframe.
#' @export
make.randomized.dat = function( n, p = 0.6, science.table.generator = make.linear.data, include.POs = TRUE, as.data.frame=TRUE, ... ) {

    dat = science.table.generator( n=n, ... )

    dat$Z = 0 + (sample(n) <= n*p)
    dat$Yobs = with( dat, ifelse( Z, Y.1, Y.0 ) )

    if ( as.data.frame ) {
        dat = as.data.frame( do.call( cbind, dat ) )
    }

    if ( !include.POs ) {
        dat$Y.0 = dat$Y.1 = dat$tau = NULL
    }

    dat
}


#' Generate fake data with noncompliance.
#'
#' This will generate and randomize a science table to get observed outcomes and treatment assignment
#'
#' @param n  Sample size
#' @param p Proportion treated
#' @param science.table.generator   Method to generate potential outcomes
#' @param include.POs  Preserve potential outcomes in returned value
#' @param ...  To be passed to science.table.generator
#'
#' @return Data frame with data randomized to tx and control, and compliers, etc.
#' @export
#' @seealso make.randomized.dat
make.randomized.compliance.dat = function( n, p = 0.6, science.table.generator = make.linear.data,
                                           include.POs = TRUE, ... ) {

    dat = make.randomized.dat( n, p, science.table.generator, include.POs = TRUE, as.data.frame = TRUE, ... )

    # will take treatment
    S = ifelse( dat$A > quantile( dat$A, 0.15 ), "C", "NT" )

    # will always take treatment
    S[ dat$B > quantile( dat$B, 0.8 ) & dat$A > median( dat$A ) ] = "AT"
    dat$S = S
    dat$D = ifelse( dat$Z == 1, ifelse( dat$S == "NT", 0, 1 ),
                                ifelse( dat$S == "AT", 1, 0 ) )

    dat$Y.1[ dat$S != "C" ] = dat$Y.0[ dat$S != "C" ]
    dat$Yobs = ifelse( dat$Z, dat$Y.1, dat$Y.0 )
    dat$tau = dat$Y.1 - dat$Y.0

    if ( ! include.POs ) {
        dat$S = dat$Y.1 = dat$Y.0 = dat$tau = NULL
    }

    dat
}

##########################################################################################
##
## Testing code
##
##########################################################################################

## if ( FALSE ) {

##     dat = make.linear.data( 10 )
##     dat
##     str(dat)


##     dat = make.skew.data( 500 )
##     hist( dat$Y.0 )
##     hist( dat$Y.1 )
##     plot( Y.1 - Y.0 ~ Y.0, data=dat )
##     par(mfrow=c(1,1) )
##     plot( data.frame( dat$X ) )


##     make.randomized.dat( 10 )
##     dt = make.randomized.dat( 10, as.data.frame=TRUE )
##     str( dt )

##     dt = make.randomized.dat( 10, as.data.frame=TRUE, include.POs = FALSE)
##     dt


##     make.randomized.dat( 10, gamma.vec=c(1,1,1,1,1,1), beta.vec=c(1,-1,1) )

## }


