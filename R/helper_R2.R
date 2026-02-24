R2_ITT <- function( RI.result, rho.step=0.05 ) {
    stopifnot( is_RI_regression_result(RI.result) )

    with( RI.result, {

        ## systematic component
        delta     <- X%*%beta.hat
        Sdd       <- as.numeric( var(delta) )

        ## idiosyncratic component
        V1        <- var(e1)
        V0        <- var(e0)

        ## approximate the intergral by summation--step quantile functions
        q.step    <- seq(from = 0, to = 1, length.out = min(length(Z), 5000))
        delta.q   <- q.step[2] - q.step[1]

        q.e1      <- quantile(e1, prob = q.step)
        q.e0      <- quantile(e0, prob = q.step)
        See.lower <- sum( (q.e1 - q.e0)^2 )*delta.q
        See.upper <- sum( (q.e1 - rev(q.e0))^2 )*delta.q

        ## R2: bounds and sensitivity analysis
        R2.lower  <- Sdd/(Sdd + See.upper)
        R2.middle <- Sdd/(Sdd + V1 + V0)
        R2.upper  <- Sdd/(Sdd + See.lower)
        rho       <- seq(0, 1, rho.step)
        R2.sensitivity <- Sdd/(Sdd + rho*See.lower + (1 - rho)*(V1 + V0))

        ## results
        res <- list(Sdd            = Sdd,
                   See.lower      = See.lower,
                   See.upper.sharp     = V0 + V1,
                   See.upper      = See.upper,
                   R2.lower       = R2.lower,
                   R2.lower.sharp      = R2.middle,
                   R2.upper       = R2.upper,
                   rho            = rho,
                   R2.sensitivity = R2.sensitivity)

        class( res ) <- "RI.R2.result"
        res$type <- "ITT"

        res
    } )
}

R2_LATE <- function( RI.result, rho.step=0.05 ) {

    with( RI.result, {
        ## tau.c: Wald estimator (LATE estimate)
        ITT <- mean(Y[Z==1]) - mean(Y[Z==0])
        tau.c <- ITT/pi.c

        ## systematic component by U
        Stautau.U <- pi.c*(1 - pi.c)*tau.c^2

        ## systematic component by X for compliers
        delta     <- as.vector(  X%*%beta.hat  )
        mean.c    <- (   mean(delta) -  sum(delta[index10])/N1 - sum(delta[index01])/N0    ) /pi.c
        deltatau2 <- (   delta - mean.c )^2


        Sdd.c     <- (   mean(deltatau2) -  sum(deltatau2[index10])/N1 - sum(deltatau2[index01])/N0    ) /pi.c
        Sdd.c     <- max(0, Sdd.c)


        ## evaluate the ECDFs at the unique values of residuals = Yunique
        Yunique   <- sort( unique(  e ) )
        F1        <- (ecdf(e11)(Yunique)*n11/N1 - ecdf(e01)(Yunique)*n01/N0)/pi.c
        F0        <- (ecdf(e00)(Yunique)*n00/N0 - ecdf(e10)(Yunique)*n10/N1)/pi.c


        ## approximate the intergral by summation--step quantile functions
        L.approx  <- min(length(Z), 5000)
        q.step    <- seq( from = 0, to = 1, length.out = L.approx )
        delta.q   <- q.step[2] - q.step[1]

        ## initial values of the integrands
        Quantile1 <- q.step
        Quantile2 <- q.step
        Quantile3 <- q.step

        for(ll in 1:L.approx)
        {
            Quantile1[ll]     <- Yunique[ which.max(F1 >= q.step[ll]) ]
            Quantile2[ll]     <- Yunique[ which.max(F0 >= q.step[ll]) ]
            Quantile3[ll]     <- Yunique[ which.max(F0 >= 1 - q.step[ll]) ]
        }


        mean1       <- sum(Quantile1)*delta.q
        mean0       <- sum(Quantile2)*delta.q
        See.c.lower <- sum( (Quantile1 - Quantile2)^2 )*delta.q - mean1^2 - mean0^2
        See.c.middle<- sum( (Quantile1)^2 )*delta.q + sum( (Quantile2)^2 )*delta.q - mean1^2 - mean0^2
        See.c.upper <- sum( (Quantile1 - Quantile3)^2 )*delta.q - mean1^2 - mean0^2


        ## three R2: bounds and sensitivity analysis
        ## R2.U
        R2.U.lower  <- Stautau.U/(  Stautau.U + pi.c*Sdd.c + pi.c*See.c.upper  )
        R2.U.middle <- Stautau.U/(  Stautau.U + pi.c*Sdd.c + pi.c*See.c.middle  )
        R2.U.upper  <- Stautau.U/(  Stautau.U + pi.c*Sdd.c + pi.c*See.c.lower  )

        rho              <- seq(0, 1, rho.step )
        R2.U.sensitivity <- Stautau.U/(  Stautau.U + pi.c*Sdd.c + pi.c*( rho*See.c.lower + (1-rho)*See.c.middle )  )

        ## R2.tau.c
        R2.lower  <- Sdd.c/(  Sdd.c + See.c.upper  )
        R2.middle <- Sdd.c/(  Sdd.c + See.c.middle  )
        R2.upper  <- Sdd.c/(  Sdd.c + See.c.lower  )

        R2.sensitivity <- Sdd.c/(  Sdd.c + rho*See.c.lower + (1-rho)*See.c.middle  )

        ## R2.UX
        R2.UX.lower     <- (Stautau.U + pi.c*Sdd.c)/(  Stautau.U + pi.c*Sdd.c + pi.c*See.c.upper  )
        R2.UX.middle    <- (Stautau.U + pi.c*Sdd.c)/(  Stautau.U + pi.c*Sdd.c + pi.c*See.c.middle  )
        R2.UX.upper     <- (Stautau.U + pi.c*Sdd.c)/(  Stautau.U + pi.c*Sdd.c + pi.c*See.c.lower  )

        R2.UX.sensitivity    <- (Stautau.U + pi.c*Sdd.c)/(  Stautau.U + pi.c*Sdd.c
                                                           + pi.c*( rho*See.c.lower + (1-rho)*See.c.middle )  )


        ## results
        res <- list(pi.c          = pi.c,
                   tau.c         = tau.c,
                   ITT           = ITT,
                   Stautau.U     = Stautau.U,
                   Sdd           = Sdd.c,
                   See.lower   = See.c.lower,
                   See.upper.sharp  = See.c.middle,
                   See.upper   = See.c.upper,

                   Yunique       = Yunique,
                   F1            = F1,
                   F0            = F0,

                   rho           = rho,

                   R2.U.lower    = R2.U.lower,
                   R2.U.lower.sharp   = R2.U.middle,
                   R2.U.upper    = R2.U.upper,
                   R2.U.sensitivity = R2.U.sensitivity,

                   R2.lower  = R2.lower,
                   R2.lower.sharp = R2.middle,
                   R2.upper  = R2.upper,
                   R2.sensitivity = R2.sensitivity,

                   R2.UX.lower     = R2.UX.lower,
                   R2.UX.lower.sharp    = R2.UX.middle,
                   R2.UX.upper     = R2.UX.upper,
                   R2.UX.sensitivity = R2.UX.sensitivity )

        res$type <- "LATE"
        class( res ) <- "RI.R2.result"

        return(res)
    } )
}


