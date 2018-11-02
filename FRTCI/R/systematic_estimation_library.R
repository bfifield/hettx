##
## Library of functions for the estimands, etc. from the JASA paper
##
## (c) 2016 - Miratrix, Deng, Feller
##


scat = function( str, ... ) {
    cat( sprintf( str, ... ) )
}


#' Calculate systematic effects model with full potential outcomes
#'
#' Calculate the actual beta for our sample using the full
#' potential outcome schedule (if available)
#'
#' This function is useful for sanity checks and simulation studies
#' to compare real estimators to optimal baselines.
#'
#' Not useful for observed data as it relys on complete knowledge of all potential outcomes.
#'
#' @param Y1 Full list of Y(1)
#' @param Y0 Full list of Y(0)
#' @param Z Realized treatment assignment
#' @param formula  Formula for the covariates X
#' @param data  Dataframe with covariates
#' @param method   RI or OLS.
#' @export
#' @seealso make.randomized.dat
calc.beta.oracle = function( Y1, Y0, Z, formula, data, method=c("RI","OLS") ) {

    Y1 = eval( substitute( Y1 ), data )
    Y0 = eval( substitute( Y0 ), data )

    X = model.matrix( formula, data )
    stopifnot( nrow(X) == length( Y1 ) )
    stopifnot( nrow(X) == length( Y0 ) )

    n = nrow(X)
    S.x0 = t(X) %*% Y0 / n
    S.x1 = t(X) %*% Y1 / n
    Sxx = t(X) %*% X / n
    Sxx.inv = solve( Sxx )
    gamma0 = Sxx.inv %*% S.x0
    gamma1 = Sxx.inv %*% S.x1
    beta.vec = gamma1 - gamma0
    nms = rownames(beta.vec)
    beta.vec = as.numeric( beta.vec )
    names(beta.vec) = nms

    #' covariance between X and tau
    tauX = (Y1 - Y0) * X
    S.xtau = apply( tauX, 2, mean )

    Z = eval( substitute( Z ), data )
    if ( !is.null( Z ) ) {
        stopifnot( length(Z) == nrow( X ) )

        #' covariance between X and Y(0)
        Y0X = Y0 * X
        S.x0 = apply(Y0X, 2, mean)

        #' covariance between X and Y(1)
        Y1X = Y1 * X
        S.x1 = apply(Y1X, 2, mean)

        Sxx = t(X) %*% X / n
        Sxx.inv = solve( Sxx )
        SE = Sxx.inv %*% ( cov( Y1X ) / sum(Z!=0) + cov( Y0X )/sum(Z==0) - cov( tauX ) / n ) %*% Sxx.inv

        SE.cons = Sxx.inv %*% ( cov( Y1X ) / sum(Z!=0) + cov( Y0X )/sum(Z==0)  ) %*% Sxx.inv
    } else {
        SE = NULL
    }

    #' Calculate individual systematic effects
    tau = X %*% beta.vec
    e1 = Y1 - X %*% gamma1
    e0 = Y0 - X %*% gamma0

    epsilon = (Y1 - Y0)  - tau
    stopifnot ( all( round( epsilon - (e1 - e0), digits=10 ) == 0 ) )

    #' return results
    res = list(beta.hat   = beta.vec,
               gamma1.hat = gamma1,
               gamma0.hat = gamma0,
               cov.beta   = SE,
               cov.beta.cons = SE.cons,
               tau.hat    = tau,
               epsilon    = epsilon,
               e1         = e1,
               e0         = e0,
               cov.tauX   = cov( tauX ),
               Sxx        = Sxx,
               chisq.stat = NA,
               p.value    = NA )

    res$method = "Oracle RI"
    res$call = match.call()
    res$X = X
    res$Y = ifelse( Z, Y1, Y0 )
    res$Y1 = Y1
    res$Y0 = Y0
    res$Z = Z
    res$tau = Y1-Y0

    class( res ) = "RI.regression.result"

    res
}





#' Estimate Systematic Treatment Effect Variation
#'
#' Estimate our treatment effect model given by the formula "~ X"
#' Will do the different estimators listed in our paper given the 'method' flag
#'
#' @param Yobs  Outcome
#' @param Z  Treatment assignment vector
#' @param formula   Formula describing what variables to model treatment effect for.
#' @param control.formula  Formula for what variables to adjust for.  If we have a control formula, we adjust
#' @param data  Dataframe with variables listed in formula and control.formula
#' @param method RI or OLS.  method=OLS is shorthand for setting the empirical.Sxx variable to TRUE, nothing more.
#' @param empirical.Sxx    Estimate seperate Sxx for treatment and control if TRUE, use known Sxx if not.
#'
#' @export
#' @seealso print.RI.regression.result
est.beta = function( Yobs, Z, formula, control.formula=NULL, data,
					method = c( "RI", "OLS" ),
					empirical.Sxx = FALSE ) {

    method = match.arg(method)
    if ( method == "OLS" ) {
        empirical.Sxx = TRUE
    }

    # if we have a control formula, we adjust
    adjust.Stx = !is.null( control.formula )

    # Get our data loaded into our environment
    Yobs = eval( substitute( Yobs ), data )
    Z = eval( substitute( Z ), data )
    X = model.matrix( formula, data )

    # Make sure our data is all of right dimension and form
    stopifnot( is.numeric( Yobs ) )
    stopifnot( is.numeric( Z ) )
    stopifnot( nrow(X) == length( Z ) )
    stopifnot( length( Yobs ) == length( Z ) )
    stopifnot( all( sort( unique( Z ) ) == c( 0, 1 ) ) )


    #' sample size
    N   = length(Z)
    N1  = sum(Z)
    N0  = N - N1

    #' dimension of X
    K   = ncol(X)

    #' treatment group
    X1  = X[Z==1, ]
    Y1  = Yobs[Z==1]

    #' control group
    X0  = X[Z==0, ]
    Y0  = Yobs[Z==0]


    Y1X      = Y1*X1
    Y0X      = Y0*X0

    S1x.hat  = apply(Y1X, 2, mean)
    S0x.hat  = apply(Y0X, 2, mean)

    # Calculate sample covariance between components of Y(1)X and also for components of Y(0)X
    if ( adjust.Stx ) {
        # We will use additional control variables for increased precision.
        W = model.matrix( control.formula, data )
        stopifnot( nrow( W ) == nrow( X ) )
        J = ncol( W )

        # utility function to estimate S.tx with model adjustment using covariate block W
        predict.YX.t = function( Y.t, X.t, W.t, W.bar ) {
            n.t = nrow(X.t)

            S.ww.t = t(W.t) %*% W.t / n.t
            YX.t = Y.t * X.t  # I.e., repeat Y.t so equiv to lhs of matrix( rep( Y.t, J ), ncol=J )

            # Regress YX.t onto W.t
            B.hat.t = solve( S.ww.t ) %*% t( W.t ) %*% YX.t / n.t

            W.bar.t = apply( W.t, 2, mean )

            # Note: this is negative of B_t'(bar{W} - W_i) from paper
            W.t.cent = W.t - rep( W.bar, each=nrow(W.t))
            YXt.hat = W.t.cent %*% B.hat.t
            YXt.hat
        }
        W1 = W[Z!=0,]
        W0 = W[Z==0,]

        W.bar = apply( W, 2, mean )

        Y1X.hat = predict.YX.t( Y1, X1, W1, W.bar )
        S1x.hat = S1x.hat - apply( Y1X.hat, 2, mean )

        Y0X.hat = predict.YX.t( Y0, X0, W0, W.bar )
        S0x.hat = S0x.hat - apply( Y0X.hat, 2, mean )
    }

    #' calculate covariance matrix of X and estimate our gammas
    if ( empirical.Sxx ) {
        Sxx.1 = t( X1 ) %*% X1 / nrow( X1 )
        Sxx.0 = t( X0 ) %*% X0 / nrow( X0 )
        Sxx.1.inv = solve( Sxx.1 )
        Sxx.0.inv = solve( Sxx.0 )
        gamma1.hat = Sxx.1.inv %*% S1x.hat
        gamma0.hat = Sxx.0.inv %*% S0x.hat
    } else {
        Sxx = t( X ) %*% X / N
        Sxx.inv = solve( Sxx )
        gamma1.hat = Sxx.inv %*% S1x.hat
        gamma0.hat = Sxx.inv %*% S0x.hat
    }
	gamma1.hat = as.numeric( gamma1.hat )
	gamma0.hat = as.numeric( gamma0.hat )


	#' residuals
    e1       = Y1 - as.vector(X1%*%gamma1.hat)
    e0       = Y0 - as.vector(X0%*%gamma0.hat)


    #' estimate beta
    beta.hat = gamma1.hat - gamma0.hat
    beta.hat = as.vector( beta.hat )
    names( gamma1.hat ) = names( gamma0.hat ) = names(beta.hat) = colnames(X)


    #' final estimator for tau, the individual systematic treatments
    tau.hat  = X%*%beta.hat
    tau.hat  = as.vector(tau.hat)


    #' covariance
    if ( empirical.Sxx ) {
        E1 = e1 * X1
        E0 = e0 * X0
    } else {
        E1 = Y1X
        E0 = Y0X
        Sxx.1.inv = Sxx.0.inv = Sxx.inv
    }

    if ( adjust.Stx ) {
        E1 = E1 - Y1X.hat
        E0 = E0 - Y0X.hat
    }

    cov.beta = Sxx.1.inv %*% ( cov(E1)/N1 ) %*% Sxx.1.inv + Sxx.0.inv %*% ( cov(E0)/N0 ) %*% Sxx.0.inv


    #' testing whether any of the non-intercept terms are nonzero
    beta1.hat   = beta.hat[2:K]
    cov.beta1   = cov.beta[2:K, 2:K]

    chisq.stat  = t(beta1.hat) %*% solve(cov.beta1, beta1.hat)

    chisq.pv    = pchisq(chisq.stat,   df = K-1, lower.tail = FALSE)

    # Some usefull stuff
    ATE = mean(Y1) - mean(Y0)
    SE.ATE = sqrt( var( Y1 ) / N1 + var( Y0 ) / N0 )


    #' return results
    res = list(Yobs       = Yobs,
               Z          = Z,
               X          = X,
               gamma1.hat = gamma1.hat,
               gamma0.hat = gamma0.hat,
               e1         = e1,
               e0         = e0,
               beta.hat   = beta.hat,
               cov.beta   = cov.beta,
               chisq.stat = chisq.stat,
               p.value    = as.numeric(chisq.pv),
               ATE        = ATE,
               SE.ATE     = SE.ATE,
               SD.Y0      = sd( Y0 ),
               SD.Y1      = sd( Y1 ) )

    res$method = paste( ifelse( empirical.Sxx, "OLS", "RI" ), ifelse( adjust.Stx, "Adjusted","Unadjusted"), sep="-" )

    res$call = match.call()
    res$X = X
    res$Y = Yobs
    res$Z = Z

    class( res ) = "RI.regression.result"

    return(res)
}













#' using the OLS output to bound the R2 measure
#'
#' @param RI.result  RI regression result (from est.beta, e.g.,) to process.
#' @param rho.step  Grid size for sensitivity analysis on values of rho.
#'
#' @export
#'
#' @return RI.R2.result object.
#' @seealso print.RI.R2.result
R2.ITT = function( RI.result, rho.step=0.05 ) {
    stopifnot( is.RI.regression.result(RI.result) )

    with( RI.result, {

        #' systematic component
        delta     = X%*%beta.hat
        Sdd       = as.numeric( var(delta) )

        #' idiosyncratic component
        V1        = var(e1)
        V0        = var(e0)

        #' approximate the intergral by summation--step quantile functions
        q.step    = seq(from = 0, to = 1, length.out = min(length(Z), 5000))
        delta.q   = q.step[2] - q.step[1]

        q.e1      = quantile(e1, prob = q.step)
        q.e0      = quantile(e0, prob = q.step)
        See.lower = sum( (q.e1 - q.e0)^2 )*delta.q
        See.upper = sum( (q.e1 - rev(q.e0))^2 )*delta.q

        #' R2: bounds and sensitivity analysis
        R2.lower  = Sdd/(Sdd + See.upper)
        R2.middle = Sdd/(Sdd + V1 + V0)
        R2.upper  = Sdd/(Sdd + See.lower)
        rho       = seq(0, 1, rho.step)
        R2.sensitivity = Sdd/(Sdd + rho*See.lower + (1 - rho)*(V1 + V0))

        #' results
        res = list(Sdd            = Sdd,
                   See.lower      = See.lower,
                   See.upper0     = V0 + V1,
                   See.upper      = See.upper,
                   R2.lower       = R2.lower,
                   R2.lower0      = R2.middle,
                   R2.upper       = R2.upper,
                   rho            = rho,
                   R2.sensitivity = R2.sensitivity)

        class( res ) = "RI.R2.result"
        res$type = "ITT"

        res
    } )
}






#####################################################################################
##
## Our LATE estimators
##
#####################################################################################


#' Estimate systematic effects in presence of noncompliance
#'
#' Estimate a model of effects for compliers under exclusion restriction for
#' never- and always-takers and under monotonicity.
#'
#' @inheritParams est.beta
#' @param D treatment received
#' @export
#'
#' @seealso est.beta
est.beta.LATE = function( Z, D, Yobs, formula, data, method=c( "RI", "2SLS" ) ) {
    method = match.arg( method )

    Yobs = eval( substitute( Yobs ), data )
    Z = as.numeric( eval( substitute( Z ), data ) )
    X = model.matrix( formula, data )
    D = as.numeric( eval( substitute( D ), data ) )

    stopifnot( is.numeric( Yobs ) )
    stopifnot( is.numeric( Z ) )
    stopifnot( all( sort( unique( Z ) ) == c( 0, 1 ) ) )
    stopifnot( all( sort( unique( D ) ) == c( 0, 1 ) ) )
    stopifnot( nrow(X) == length( Yobs ) )
    stopifnot( nrow(X) == length( D ) )
    stopifnot( nrow(X) == length( Z ) )


    #' sample size, and sub-sample sizes classified by (Z, D)
    N        = length(Z)
    N1       = sum(Z)
    N0       = N - N1
    K        = ncol(X)

    index11  = (1:N)[Z==1&D==1]
    index10  = (1:N)[Z==1&D==0]
    index01  = (1:N)[Z==0&D==1]
    index00  = (1:N)[Z==0&D==0]
    n11      = length(index11)
    n10      = length(index10)
    n01      = length(index01)
    n00      = length(index00)
    stopifnot( min( n11, n10, n01, n00 ) > 1 )

    Y11      = Yobs[index11]
    Y10      = Yobs[index10]
    Y01      = Yobs[index01]
    Y00      = Yobs[index00]
    X11      = X[index11, ]
    X10      = X[index10, ]
    X01      = X[index01, ]
    X00      = X[index00, ]

    #' estimation of the proportions of the latent strata
    pi.n     = n10/N1
    pi.a     = n01/N0
    pi.c     = n11/N1 - n01/N0

    Sxxc1     = t(X11)%*%X11/N1 - t(X01)%*%X01/N0
    Sxxc1.inv = solve(Sxxc1)

    Sxxc0     = t(X00)%*%X00/N0 - t(X10)%*%X10/N1
    Sxxc0.inv = solve(Sxxc0)

    if ( method=="2SLS" ) {
        #' TSLS
        Sxx      = t(X)%*%X/N
        Sxx.D    = t(X[D==1, ])%*%(X[D==1, ])/N
        Sxx.T    = t(X[Z==1, ])%*%(X[Z==1, ])/N
        Sxx.TD   = t(X11)%*%X11/N
        S.inv    = solve(  rbind( cbind(Sxx, Sxx.D), cbind(Sxx.T, Sxx.TD) )   )
        rhs      = c( as.vector( apply(Yobs*X, 2, mean) ), as.vector( apply(Z*Yobs*X, 2, mean) )  )
        TSLS     = as.vector( S.inv%*%rhs )
        gamma.hat = TSLS[1:K]
        beta.hat  = TSLS[(K+1):(2*K)]

        #' residual
        e        = Yobs - as.vector(X%*%gamma.hat) - D*as.vector(X%*%beta.hat)
        e        = as.vector(e)

    } else {

        Sx1.11 = apply(Y11*X11, 2, sum)/N1
        Sx0.01 = apply(Y01*X01, 2, sum)/N0

        #' gamma for compliers: under treatment
        gamma1.hat  = Sxxc1.inv %*% ( Sx1.11 - Sx0.01 )
        gamma1.hat  = as.vector(gamma1.hat)

        Sx0.00 = apply(Y00*X00, 2, sum)/N0
        Sx1.10 = apply(Y10*X10, 2, sum)/N1

        #' gamma for compliers: under control
        gamma0.hat  = Sxxc0.inv %*% (Sx0.00 - Sx1.10)
        gamma0.hat  = as.vector(gamma0.hat)


        #' unbiased estimator for beta
        beta.hat = gamma1.hat - gamma0.hat

        #' residual
        e        = Yobs - ifelse( D, as.vector(X%*%gamma1.hat), as.vector(X%*%gamma0.hat) )
        e        = as.vector(e)
    }

    #' covariance of beta
    X1       = X[Z==1, ]
    e1       = e[Z==1]
    X0       = X[Z==0, ]
    e0       = e[Z==0]

    cov.beta = Sxxc1.inv%*%( cov(e1*X1)/N1 )%*%Sxxc1.inv + Sxxc0.inv%*%(  cov(e0*X0)/N0  )%*%Sxxc0.inv
    rownames( cov.beta ) = colnames( cov.beta ) = colnames( X )

        #' testing
        beta1.hat   = beta.hat[2:K]
        cov.beta1   = cov.beta[2:K, 2:K]

        chisq.stat  = t(beta1.hat)%*%solve(cov.beta1, beta1.hat)

        chisq.pv    = pchisq(chisq.stat,   df = K-1, lower.tail = FALSE)

        beta.hat = as.numeric( beta.hat )
        names( beta.hat ) = colnames( X )

    if ( method=="2SLS" ) {
        gamma.hat = as.numeric( gamma.hat )
        names( gamma.hat ) = colnames( X )
        gamma1.hat = gamma0.hat = NA
    } else {
        gamma1.hat = as.numeric( gamma1.hat )
        gamma0.hat = as.numeric( gamma0.hat )
        names( gamma0.hat ) = names( gamma1.hat ) = colnames( X )
        gamma.hat = NA
    }

    #' return
    res = list(Yobs       = Yobs,
               Z          = Z,
               D          = D,
               X          = X,
               N          = N,
               N1         = N1,
               N0         = N0,
               index11    = index11,
               index10    = index10,
               index01    = index01,
               index00    = index00,
               n11        = n11,
               n10        = n10,
               n01        = n01,
               n00        = n00,
               pi.a       = pi.a,
               pi.n       = pi.n,
               pi.c       = pi.c,
               gamma.hat  = gamma.hat,
               gamma1.hat = gamma1.hat,
               gamma0.hat = gamma0.hat,
               e          = e,
               e11        = e[index11],
               e10        = e[index10],
               e01        = e[index01],
               e00        = e[index00],
               beta.hat   = beta.hat,
               cov.beta   = cov.beta,
               chisq.stat = chisq.stat,
               p.value         = as.numeric(chisq.pv))

    res$Y = Yobs
    res$Z = Z
    res$D = D
    res$X = X

    res$method = paste( "LATE RI-", method, sep="" )
    res$call = match.call()

    class( res ) = "RI.regression.result"

    return(res)
}




#' using the LATE estimation output to bound the R2 measures for compliers, etc.
#'
#' @inheritParams R2.ITT
#'
#' @export
#' @seealso R2.ITT
R2.LATE = function( RI.result, rho.step=0.05 ) {

    with( RI.result, {
        #' tau.c: Wald estimator (LATE estimate)
        ITT = mean(Y[Z==1]) - mean(Y[Z==0])
        tau.c = ITT/pi.c

        #' systematic component by U
        Stautau.U = pi.c*(1 - pi.c)*tau.c^2

        #' systematic component by X for compliers
        delta     = as.vector(  X%*%beta.hat  )
        mean.c    = (   mean(delta) -  sum(delta[index10])/N1 - sum(delta[index01])/N0    ) /pi.c
        deltatau2 = (   delta - mean.c )^2


        Sdd.c     = (   mean(deltatau2) -  sum(deltatau2[index10])/N1 - sum(deltatau2[index01])/N0    ) /pi.c
        Sdd.c     = max(0, Sdd.c)


        #' evaluate the ECDFs at the unique values of residuals = Yunique
        Yunique   = sort( unique(  e ) )
        F1        = (ecdf(e11)(Yunique)*n11/N1 - ecdf(e01)(Yunique)*n01/N0)/pi.c
        F0        = (ecdf(e00)(Yunique)*n00/N0 - ecdf(e10)(Yunique)*n10/N1)/pi.c


        #' approximate the intergral by summation--step quantile functions
        L.approx  = min(length(Z), 5000)
        q.step    = seq( from = 0, to = 1, length.out = L.approx )
        delta.q   = q.step[2] - q.step[1]

        #' initial values of the integrands
        Quantile1 = q.step
        Quantile2 = q.step
        Quantile3 = q.step

        for(ll in 1:L.approx)
        {
            Quantile1[ll]     = Yunique[ which.max(F1 >= q.step[ll]) ]
            Quantile2[ll]     = Yunique[ which.max(F0 >= q.step[ll]) ]
            Quantile3[ll]     = Yunique[ which.max(F0 >= 1 - q.step[ll]) ]
        }


        mean1       = sum(Quantile1)*delta.q
        mean0       = sum(Quantile2)*delta.q
        See.c.lower = sum( (Quantile1 - Quantile2)^2 )*delta.q - mean1^2 - mean0^2
        See.c.middle= sum( (Quantile1)^2 )*delta.q + sum( (Quantile2)^2 )*delta.q - mean1^2 - mean0^2
        See.c.upper = sum( (Quantile1 - Quantile3)^2 )*delta.q - mean1^2 - mean0^2


        #' three R2: bounds and sensitivity analysis
        #' R2.U
        R2.U.lower  = Stautau.U/(  Stautau.U + pi.c*Sdd.c + pi.c*See.c.upper  )
        R2.U.middle = Stautau.U/(  Stautau.U + pi.c*Sdd.c + pi.c*See.c.middle  )
        R2.U.upper  = Stautau.U/(  Stautau.U + pi.c*Sdd.c + pi.c*See.c.lower  )

        rho              = seq(0, 1, rho.step )
        R2.U.sensitivity = Stautau.U/(  Stautau.U + pi.c*Sdd.c + pi.c*( rho*See.c.lower + (1-rho)*See.c.middle )  )

        #' R2.tau.c
        R2.lower  = Sdd.c/(  Sdd.c + See.c.upper  )
        R2.middle = Sdd.c/(  Sdd.c + See.c.middle  )
        R2.upper  = Sdd.c/(  Sdd.c + See.c.lower  )

        R2.sensitivity = Sdd.c/(  Sdd.c + rho*See.c.lower + (1-rho)*See.c.middle  )

        #' R2.UX
        R2.UX.lower     = (Stautau.U + pi.c*Sdd.c)/(  Stautau.U + pi.c*Sdd.c + pi.c*See.c.upper  )
        R2.UX.middle    = (Stautau.U + pi.c*Sdd.c)/(  Stautau.U + pi.c*Sdd.c + pi.c*See.c.middle  )
        R2.UX.upper     = (Stautau.U + pi.c*Sdd.c)/(  Stautau.U + pi.c*Sdd.c + pi.c*See.c.lower  )

        R2.UX.sensitivity    = (Stautau.U + pi.c*Sdd.c)/(  Stautau.U + pi.c*Sdd.c
                                                           + pi.c*( rho*See.c.lower + (1-rho)*See.c.middle )  )


        #' results
        res = list(pi.c          = pi.c,
                   tau.c         = tau.c,
                   ITT           = ITT,
                   Stautau.U     = Stautau.U,
                   Sdd           = Sdd.c,
                   See.lower   = See.c.lower,
                   See.upper0  = See.c.middle,
                   See.upper   = See.c.upper,

                   Yunique       = Yunique,
                   F1            = F1,
                   F0            = F0,

                   rho           = rho,

                   R2.U.lower    = R2.U.lower,
                   R2.U.lower0   = R2.U.middle,
                   R2.U.upper    = R2.U.upper,
                   R2.U.sensitivity = R2.U.sensitivity,

                   R2.lower  = R2.lower,
                   R2.lower0 = R2.middle,
                   R2.upper  = R2.upper,
                   R2.sensitivity = R2.sensitivity,

                   R2.UX.lower     = R2.UX.lower,
                   R2.UX.lower0    = R2.UX.middle,
                   R2.UX.upper     = R2.UX.upper,
                   R2.UX.sensitivity = R2.UX.sensitivity )

        res$type = "LATE"
        class( res ) = "RI.R2.result"

        return(res)
    } )
}



#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'
#'
#' Clean printing of our result objects
#'
#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'#'


#' Check if object is a RI.regression.result class
#' @param x Object to test.
#'
#' @export
is.RI.regression.result = function( x ) {
    inherits(x, "RI.regression.result")
}

#' Extract coefficients of a fit RI regression model.
#' @param object A RI.regression.result object.
#' @param ... Unused
#'
#' @export
coef.RI.regression.result = function( object, ... ) {
    object$beta.hat
}


#' Print out result of systematic treatment effect estimation.
#'
#' @param x  Result from est.beta
#'
#' @param digits  Number of digits for rounding
#' @param ...  Unused
#'
#' @export
print.RI.regression.result = function( x, digits = max(3L, getOption("digits") - 3L), ... ) {

    scat(  "RI het tx regression (%s):\n\t%s\n\n", x$method, paste( deparse( x$call ), sep="\n\t", collapse="\n\t" ) )
    #browser()
    #scat(  "RI het tx regression (%s):\n\n", x$method )

    if (length(coef(x))) {
        cat("Coefficients:\n")
        print.default(format( coef(x), digits = digits), print.gap = 2L,
                      quote = FALSE)
    }
    else cat("No coefficients\n")
    cat("\nVar-Cor Matrix:\n")
    print( x$cov.beta )
    scat(  "\nChi-squared test for systematic variation: X^2=%.2f ; p=%.3f\n", x$chisq.stat, x$p.value )

    scat( "\nDetails:  ATE = %.3f +/- %.3f     SD( Y(0) ) = %.3f   SD( Y(1) ) = %.3f", x$ATE, 1.96 * x$SE.ATE, x$SD.Y0, x$SD.Y1 )

    invisible( x )
}

#' Get vcov() from object.
#'
#' @param object est.beta object
#' @param ... unused
#'
#' @export
vcov.RI.regression.result = function( object, ... ) {
    object$cov.beta
}

#' Extract the standard errors from a var-cov matrix.
#'
#' @param object est.beta object
#' @param ... unused
#'
#' @export
SE = function( object, ... ) {
    sqrt( diag( vcov( object ) ) )
}

#' Print out results on estimating treatment effect R2
#'
#' Print various R2 measures including the bounds due to the sensitivity analysis.
#'
#' @param x  Object to print
#' @param ...  Ignored
#'
#' @export
print.RI.R2.result = function( x, ... ) {
    scat(  "Heterogenous Tx Effect R^2 (%s)\n", x$type )

    if ( x$type == "LATE" ) {
        scat( "\tR^2 (within Compliers):\t\t %.3f | %.3f -- %.3f\n", x$R2.lower, x$R2.lower0, x$R2.upper )

        scat( "\tR^2 for noncompliance alone:\t %.3f | %.3f -- %.3f\n", x$R2.U.lower, x$R2.U.lower0, x$R2.U.upper )

        scat( "\tR^2 for covariates & compliance: %.3f | %.3f -- %.3f\n", x$R2.UX.lower, x$R2.UX.lower0, x$R2.UX.upper )

        scat( "\n  Variances:\n" )
        scat( "\tSystematic Tx Var (Compliers):\t %.3f\n", x$Sdd  )
        scat( "\tSystematic Tx Var (Strata):\t %.3f\n", x$Stautau.U )
        totsys = x$pi.c * x$Sdd + x$Stautau.U
        scat( "\tTotal Systematic Var:\t\t %.3f\n", totsys )
        scat( "\tIdeosyncratic Tx Var (Comp):\t %.3f -- %.3f | %0.3f\n",  x$See.lower, x$See.upper0, x$See.upper )
        scat( "\tTotal variation:\t\t %.3f -- %.3f | %.3f\n", totsys + x$pi.c*x$See.lower, totsys+x$pi.c*x$See.upper0, totsys+x$pi.c*x$See.upper )
        scat( "\n  Details: \tLATE = %.3f   ITT = %.3f    Proportion compliers = %.3f\n", x$tau.c, x$ITT, x$pi.c )
    } else {
        scat( "\tR^2:\t\t\t %.3f | %.3f -- %.3f\n", x$R2.lower, x$R2.lower0, x$R2.upper  )
        scat( "\n  Treatment Effect Variances:\n" )
        scat( "\tSystematic:\t %.3f\n", x$Sdd  )
        scat( "\tIdeosyncratic:\t %.3f -- %.3f | %0.3f\n", x$See.lower, x$See.upper0, x$See.upper )
        scat( "\tTotal:\t %.3f -- %.3f | %0.3f\n", x$Sdd + x$See.lower, x$Sdd  +x$See.upper0, x$Sdd + x$See.upper )


    }
}


#' Make a plot of the treatment effect R2 estimates
#'
#'
#' @param x Results from est.beta, etc.
#' @param main Title for plot
#' @param ADD TRUE if add to existing plot. FALSE make a new plot.
#' @param ...  Arguments to pass to plotting of points.
#'
#' @export
#' @seealso calc.beta
plot.RI.R2.result = function( x, main=paste( "R2 for Het Tx (", x$type, ")", sep=""),
                              ADD=FALSE, ... ) {
    with( x, {
        if ( !ADD ) {
            plot( R2.sensitivity ~ rho, type="l", ylab="R2 (Tx variation explained)", ylim=c(0,1), xlim=c(0,1), lwd=2, main=main, ... )
        }
        points( R2.sensitivity ~ rho, type="l", lwd=2, ... )
        points( range( R2.sensitivity ) ~ range( rho ), pch=19 )
    } )

    #abline( v=0, lwd=1, col="red" )
    abline( h=x$R2.lower, col="red" )
    invisible( x )
}





#' Variance ratio test
#'
#' Given vector of observed outcomes and treatment vector, test to see if there is evidence
#' the variances are different (taking kurtosis into account).
#'
#' @param Yobs  Outcome
#' @param Z  Treatment assignment vector
#' @param data  Dataframe with variables listed in formula and control.formula
#' @export
#' @importFrom moments kurtosis
variance.ratio.test = function(Yobs, Z, data= NULL)
{
  if (!is.null( data ) ) {
    Yobs = eval( substitute( Yobs ), data )
    Z = eval( substitute( Z ), data )
  }
  Y1 = Yobs[Z==1]
  Y0 = Yobs[Z==0]

  N1 = length(Y1)
  N0 = length(Y0)

  log.varR = log(var(Y1)/var(Y0))

  asy.se   = sqrt(  (kurtosis(Y1) - 1)/N1 + (kurtosis(Y0) - 1)/N0   )

  pvalue   = as.numeric((1 - pnorm(abs(log.varR), 0, asy.se)))

  res = data.frame( pvalue = pvalue, var1 = var( Y1 ), var0 = var( Y0 ))
  res$ratio = res$var1 / res$var0
  res$log.ratio = log( res$ratio )
  res$asy.se = asy.se
  res$z = res$log.ratio / res$asy.se
  return( res )
}




