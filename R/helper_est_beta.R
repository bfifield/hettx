calc_beta_oracle <- function( formula, data, interaction.formula, method=c("RI","OLS"), na.rm = FALSE ) {

    ## ---------------------------
    ## Get variables from formulas
    ## and check formulas
    ## ---------------------------
    if(inherits(data, "data.frame")) {
        data <- as.data.frame(data)
    }else{
        stop("The input data must be of class tbl_df, data.table, or data.frame.")
    }

    if(length(lhs_vars(formula)) != 2 | length(rhs_vars(formula)) != 1){
        stop("The formula argument must be of the form treated_outcome + control_outcome ~ treatment.")
    }
    main.vars <- get_vars(formula,data=data)
    if(any(!(main.vars %in% colnames(data)))){
        stop("Some variables in formula are not present in your data.")
    }

    ## Interaction formula
    if(length(lhs_vars(interaction.formula)) != 0){
        stop("Do not provide an outcome variable in interaction.formula.")
    }
    interaction.vars <- get_vars(interaction.formula,data=data)
    if(any(!(interaction.vars %in% colnames(data)))){
        stop("Some variables in interaction.formula are not present in your data.")
    }
    if(any(main.vars %in% interaction.vars)){
        stop("Either your outcome variable or your treatment variable is present in your interaction formula.")
    }

    ## NA handling
    vars <- c(main.vars, interaction.vars)
    if(na.rm == TRUE){
        data <- na.omit(data[,vars])
    }else{
        if(sum(is.na(data[,vars])) > 0){
            stop("You have missing values in your input data. Either set na.rm = TRUE or check your input data for missingness.")
        }
    }

    ## Extract data
    Y1 <- data[,main.vars[1]]
    Y0 <- data[,main.vars[2]]
    Z <- data[,main.vars[3]]
    X <- model.matrix(interaction.formula, data)

    stopifnot( nrow(X) == length( Y1 ) )
    stopifnot( nrow(X) == length( Y0 ) )

    n <- nrow(X)
    S.x0 <- t(X) %*% Y0 / n
    S.x1 <- t(X) %*% Y1 / n
    Sxx <- t(X) %*% X / n
    Sxx.inv <- solve( Sxx )
    gamma0 <- Sxx.inv %*% S.x0
    gamma1 <- Sxx.inv %*% S.x1
    beta.vec <- gamma1 - gamma0
    nms <- rownames(beta.vec)
    beta.vec <- as.numeric( beta.vec )
    names(beta.vec) <- nms

    ## covariance between X and tau
    tauX <- (Y1 - Y0) * X
    S.xtau <- apply( tauX, 2, mean )

    Z <- eval( substitute( Z ), data )
    if ( !is.null( Z ) ) {
        stopifnot( length(Z) == nrow( X ) )

        ## covariance between X and Y(0)
        Y0X <- Y0 * X
        S.x0 <- apply(Y0X, 2, mean)

        ## covariance between X and Y(1)
        Y1X <- Y1 * X
        S.x1 <- apply(Y1X, 2, mean)

        Sxx <- t(X) %*% X / n
        Sxx.inv <- solve( Sxx )
        SE <- Sxx.inv %*% ( cov( Y1X ) / sum(Z!=0) + cov( Y0X )/sum(Z==0) - cov( tauX ) / n ) %*% Sxx.inv

        SE.cons <- Sxx.inv %*% ( cov( Y1X ) / sum(Z!=0) + cov( Y0X )/sum(Z==0)  ) %*% Sxx.inv
    } else {
        SE <- NULL
    }

    ## Calculate individual systematic effects
    tau <- X %*% beta.vec
    e1 <- Y1 - X %*% gamma1
    e0 <- Y0 - X %*% gamma0

    epsilon <- (Y1 - Y0)  - tau
    stopifnot ( all( round( epsilon - (e1 - e0), digits=10 ) == 0 ) )

    ## return results
    res <- list(beta.hat   = beta.vec,
               gamma1.hat = gamma1,
               gamma0.hat = gamma0,
               cov.beta   = SE,
               cov.beta.cons = SE.cons,
               tau.hat    = tau,
               epsilon    = epsilon,
               e1         = as.vector( e1 ),
               e0         = as.vector( e0 ),
               cov.tauX   = cov( tauX ),
               Sxx        = Sxx,
               chisq.stat = NA,
               p.value    = NA )

    res$method <- "Oracle RI"
    res$X <- X
    res$Y <- ifelse( Z, Y1, Y0 )
    res$Y1 <- Y1
    res$Y0 <- Y0
    res$Z <- Z
    res$tau <- Y1-Y0

    class( res ) <- c("RI.regression.result", "RI.regression.result.oracle")

    res
}


est_beta_ITT <- function( formula, data, interaction.formula, control.formula=NULL,
                        method = c( "RI", "OLS" ),
                        na.rm = FALSE) {

    empirical.Sxx <- FALSE
    method <- match.arg(method)
    if ( method == "OLS" ) {
        empirical.Sxx <- TRUE
    }

    # if we have a control formula, we adjust
    adjust.Stx <- !is.null( control.formula )

    ## -------------------------
    ## Set up data with formulas
    ## -------------------------
    ## Data to data.frame
    if(inherits(data, "data.frame")) {
        data <- as.data.frame(data)
    }else{
        stop("The input data must be of class tbl_df, data.table, or data.frame.")
    }

    ## Main formula
    if(length(lhs_vars(formula)) != 1 | length(rhs_vars(formula)) != 1){
        stop("The formula argument must be of the form outcome ~ treatment.")
    }
    main.vars <- get_vars(formula,data=data)
    if(any(!(main.vars %in% colnames(data)))){
        stop("Some variables in formula are not present in your data.")
    }

    ## Interaction formula
    if(length(lhs_vars(interaction.formula)) != 0){
        stop("Do not provide an outcome variable in interaction.formula.")
    }
    interaction.vars <- get_vars(interaction.formula,data=data)
    if(any(!(interaction.vars %in% colnames(data)))){
        stop("Some variables in interaction.formula are not present in your data.")
    }
    if(any(main.vars %in% interaction.vars)){
        stop("Either your outcome variable or your treatment variable is present in your interaction formula.")
    }

    ## Control formula
    if(!is.null(control.formula)){
        if(length(lhs_vars(control.formula)) != 0){
            stop("Do not provide an outcome variable in control.formula.")
        }
        control.vars <- get_vars(control.formula,data=data)
        if(any(!(control.vars %in% colnames(data)))){
            stop("Some variables in control.formula are not present in your data.")
        }
        if(any(main.vars %in% interaction.vars)){
            stop("Either your outcome variable or your treatment variable is present in your control formula.")
        }
    }else{
        control.vars <- NULL
    }

    ## NA handling
    vars <- c(main.vars, interaction.vars, control.vars)
    if(na.rm == TRUE){
        data <- na.omit(data[,vars])
    }else{
        if(sum(is.na(data[,vars])) > 0){
            stop("You have missing values in your input data. Either set na.rm = TRUE or check your input data for missingness.")
        }
    }

    ## Extract data
    Yobs <- data[,main.vars[1]]
    Z <- data[,main.vars[2]]
    X <- model.matrix(interaction.formula, data)

    # Make sure our data is all of right dimension and form
    stopifnot( is.numeric( Yobs ) )
    stopifnot( is.numeric( Z ) )
    stopifnot( nrow(X) == length( Z ) )
    stopifnot( length( Yobs ) == length( Z ) )
    if( !(all( sort(unique(Z)) == c(0,1)) ) ) {
        stop("Your treatment variable must only take values of 0 and 1.")
    }


    ## sample size
    N   <- length(Z)
    N1  <- sum(Z)
    N0  <- N - N1

    ## dimension of X
    K   <- ncol(X)

    ## treatment group
    X1  <- X[Z==1, ]
    Y1  <- Yobs[Z==1]

    ## control group
    X0  <- X[Z==0, ]
    Y0  <- Yobs[Z==0]


    Y1X      <- Y1*X1
    Y0X      <- Y0*X0

    S1x.hat  <- apply(Y1X, 2, mean)
    S0x.hat  <- apply(Y0X, 2, mean)

    # Calculate sample covariance between components of Y(1)X and also for components of Y(0)X
    if ( adjust.Stx ) {
        # We will use additional control variables for increased precision.
        W <- model.matrix(control.formula, data)
        stopifnot( nrow( W ) == nrow( X ) )
        J <- ncol( W )

        # utility function to estimate S.tx with model adjustment using covariate block W
        predict.YX.t <- function( Y.t, X.t, W.t, W.bar ) {
            n.t <- nrow(X.t)

            S.ww.t <- t(W.t) %*% W.t / n.t
            YX.t <- Y.t * X.t  # I.e., repeat Y.t so equiv to lhs of matrix( rep( Y.t, J ), ncol=J )

            # Regress YX.t onto W.t
            B.hat.t <- solve( S.ww.t ) %*% t( W.t ) %*% YX.t / n.t

            W.bar.t <- apply( W.t, 2, mean )

            # Note: this is negative of B_t'(bar{W} - W_i) from paper
            W.t.cent <- W.t - rep( W.bar, each=nrow(W.t))
            YXt.hat <- W.t.cent %*% B.hat.t
            YXt.hat
        }
        W1 <- W[Z!=0,]
        W0 <- W[Z==0,]

        W.bar <- apply( W, 2, mean )

        Y1X.hat <- predict.YX.t( Y1, X1, W1, W.bar )
        S1x.hat <- S1x.hat - apply( Y1X.hat, 2, mean )

        Y0X.hat <- predict.YX.t( Y0, X0, W0, W.bar )
        S0x.hat <- S0x.hat - apply( Y0X.hat, 2, mean )
    }

    ## calculate covariance matrix of X and estimate our gammas
    if ( empirical.Sxx ) {
        Sxx.1 <- t( X1 ) %*% X1 / nrow( X1 )
        Sxx.0 <- t( X0 ) %*% X0 / nrow( X0 )
        Sxx.1.inv <- solve( Sxx.1 )
        Sxx.0.inv <- solve( Sxx.0 )
        gamma1.hat <- Sxx.1.inv %*% S1x.hat
        gamma0.hat <- Sxx.0.inv %*% S0x.hat
    } else {
        Sxx <- t( X ) %*% X / N
        Sxx.inv <- solve( Sxx )
        gamma1.hat <- Sxx.inv %*% S1x.hat
        gamma0.hat <- Sxx.inv %*% S0x.hat
    }
    gamma1.hat <- as.numeric( gamma1.hat )
    gamma0.hat <- as.numeric( gamma0.hat )

    ## residuals
    e1       <- Y1 - as.vector(X1%*%gamma1.hat)
    e0       <- Y0 - as.vector(X0%*%gamma0.hat)


    ## estimate beta
    beta.hat <- gamma1.hat - gamma0.hat
    beta.hat <- as.vector( beta.hat )
    names( gamma1.hat ) <- names( gamma0.hat ) <- names(beta.hat) <- colnames(X)


    ## final estimator for tau, the individual systematic treatments
    tau.hat  <- X%*%beta.hat
    tau.hat  <- as.vector(tau.hat)


    ## covariance
    if ( empirical.Sxx ) {
        E1 <- e1 * X1
        E0 <- e0 * X0
    } else {
        E1 <- Y1X
        E0 <- Y0X
        Sxx.1.inv <- Sxx.0.inv <- Sxx.inv
    }

    if ( adjust.Stx ) {
        E1 <- E1 - Y1X.hat
        E0 <- E0 - Y0X.hat
    }

    cov.beta <- Sxx.1.inv %*% ( cov(E1)/N1 ) %*% Sxx.1.inv + Sxx.0.inv %*% ( cov(E0)/N0 ) %*% Sxx.0.inv


    ## testing whether any of the non-intercept terms are nonzero
    beta1.hat   <- beta.hat[2:K]
    cov.beta1   <- cov.beta[2:K, 2:K]

    chisq.stat  <- t(beta1.hat) %*% solve(cov.beta1, beta1.hat)

    chisq.pv    <- pchisq(chisq.stat,   df = K-1, lower.tail = FALSE)

    # Some usefull stuff
    ATE <- mean(Y1) - mean(Y0)
    SE.ATE <- sqrt( var( Y1 ) / N1 + var( Y0 ) / N0 )


    ## return results
    res <- list(Yobs       = Yobs,
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

    res$method <- paste( ifelse( empirical.Sxx, "OLS", "RI" ),
                        ifelse( adjust.Stx, "Adjusted","Unadjusted"), sep="-" )

    res$X <- X
    res$Y <- Yobs
    res$Z <- Z
    res$control.vars <- control.vars
    res$main.vars <- main.vars
    res$interaction.vars <- interaction.vars

    class( res ) <- c("RI.regression.result", "RI.regression.result.ITT")

    return(res)
}


est_beta_LATE <- function(formula, data, interaction.formula,
                          method=c("RI", "2SLS"), na.rm = TRUE){
    method <- match.arg( method )

    ## -------------------------
    ## Set up data with formulas
    ## -------------------------
    ## Data to data.frame
    if(inherits(data, "data.frame")) {
        data <- as.data.frame(data)
    }else{
        stop("The input data must be of class tbl_df, data.table, or data.frame.")
    }

    ## Main formula
    formula.char <- paste( deparse(formula), collapse=" " )
    formula.char <- gsub( "\\s+", " ", formula.char, perl=FALSE )
    if(grepl("\\|", formula.char)){
        formula <- as.formula(gsub("[>|]", "+", formula.char))
    }else{
        stop("The formula must be of the form outcome ~ treatment | instrument.")
    }
    if(length(lhs_vars(formula)) != 1 | length(rhs_vars(formula)) != 2){
        stop("The formula argument must be of the form outcome ~ treatment | instrument.")
    }
    main.vars <- get_vars(formula,data=data)
    if(any(!(main.vars %in% colnames(data)))){
        stop("Some variables in formula are not present in your data.")
    }

    ## Interaction formula
    if(length(lhs_vars(interaction.formula)) != 0){
        stop("Do not provide an outcome variable in interaction.formula.")
    }
    interaction.vars <- get_vars(interaction.formula,data=data)
    if(any(!(interaction.vars %in% colnames(data)))){
        stop("Some variables in interaction.formula are not present in your data.")
    }
    if(any(main.vars %in% interaction.vars)){
        stop("Either your outcome variable, your treatment variable, or your instrument variable is present in your interaction formula.")
    }

    ## NA handling
    vars <- c(main.vars, interaction.vars)
    if(na.rm == TRUE){
        data <- na.omit(data[,vars])
    }else{
        if(sum(is.na(data[,vars])) > 0){
            stop("You have missing values in your input data. Either set na.rm = TRUE or check your input data for missingness.")
        }
    }

    Yobs <- data[,main.vars[1]]
    D <- data[,main.vars[2]]
    Z <- data[,main.vars[3]]
    X <- model.matrix(interaction.formula, data)

    if( !(all( sort(unique(Z)) == c(0,1) ) ) ){
        stop("Your instrument variable must only take values of 0 and 1.")
    }
    if(!(all( sort(unique(D)) == c(0,1)) ) ){
        stop("Your treatment variable must only take values of 0 and 1.")
    }

    stopifnot( is.numeric( Yobs ) )
    stopifnot( is.numeric( Z ) )
    stopifnot( nrow(X) == length( Yobs ) )
    stopifnot( nrow(X) == length( D ) )
    stopifnot( nrow(X) == length( Z ) )


    ## sample size, and sub-sample sizes classified by (Z, D)
    N        <- length(Z)
    N1       <- sum(Z)
    N0       <- N - N1
    K        <- ncol(X)

    index11  <- (1:N)[Z==1&D==1]
    index10  <- (1:N)[Z==1&D==0]
    index01  <- (1:N)[Z==0&D==1]
    index00  <- (1:N)[Z==0&D==0]
    n11      <- length(index11)
    n10      <- length(index10)
    n01      <- length(index01)
    n00      <- length(index00)
    stopifnot( min( n11, n10, n01, n00 ) > 1 )

    Y11      <- Yobs[index11]
    Y10      <- Yobs[index10]
    Y01      <- Yobs[index01]
    Y00      <- Yobs[index00]
    X11      <- X[index11, ]
    X10      <- X[index10, ]
    X01      <- X[index01, ]
    X00      <- X[index00, ]

    ## estimation of the proportions of the latent strata
    pi.n     <- n10/N1
    pi.a     <- n01/N0
    pi.c     <- n11/N1 - n01/N0

    Sxxc1     <- t(X11)%*%X11/N1 - t(X01)%*%X01/N0
    Sxxc1.inv <- solve(Sxxc1)

    Sxxc0     <- t(X00)%*%X00/N0 - t(X10)%*%X10/N1
    Sxxc0.inv <- solve(Sxxc0)

    if ( method=="2SLS" ) {
        ## TSLS
        Sxx      <- t(X)%*%X/N
        Sxx.D    <- t(X[D==1, ])%*%(X[D==1, ])/N
        Sxx.T    <- t(X[Z==1, ])%*%(X[Z==1, ])/N
        Sxx.TD   <- t(X11)%*%X11/N
        S.inv    <- solve(  rbind( cbind(Sxx, Sxx.D), cbind(Sxx.T, Sxx.TD) )   )
        rhs      <- c( as.vector( apply(Yobs*X, 2, mean) ), as.vector( apply(Z*Yobs*X, 2, mean) )  )
        TSLS     <- as.vector( S.inv%*%rhs )
        gamma.hat <- TSLS[1:K]
        beta.hat  <- TSLS[(K+1):(2*K)]

        ## residual
        e        <- Yobs - as.vector(X%*%gamma.hat) - D*as.vector(X%*%beta.hat)
        e        <- as.vector(e)

    } else {

        Sx1.11 <- apply(Y11*X11, 2, sum)/N1
        Sx0.01 <- apply(Y01*X01, 2, sum)/N0

        ## gamma for compliers: under treatment
        gamma1.hat  <- Sxxc1.inv %*% ( Sx1.11 - Sx0.01 )
        gamma1.hat  <- as.vector(gamma1.hat)

        Sx0.00 <- apply(Y00*X00, 2, sum)/N0
        Sx1.10 <- apply(Y10*X10, 2, sum)/N1

        ## gamma for compliers: under control
        gamma0.hat  <- Sxxc0.inv %*% (Sx0.00 - Sx1.10)
        gamma0.hat  <- as.vector(gamma0.hat)


        ## unbiased estimator for beta
        beta.hat <- gamma1.hat - gamma0.hat

        ## residual
        e        <- Yobs - ifelse( D, as.vector(X%*%gamma1.hat), as.vector(X%*%gamma0.hat) )
        e        <- as.vector(e)
    }

    ## covariance of beta
    X1       <- X[Z==1, ]
    e1       <- e[Z==1]
    X0       <- X[Z==0, ]
    e0       <- e[Z==0]

    cov.beta <- Sxxc1.inv%*%( cov(e1*X1)/N1 )%*%Sxxc1.inv + Sxxc0.inv%*%(  cov(e0*X0)/N0  )%*%Sxxc0.inv
    rownames( cov.beta ) <- colnames( cov.beta ) <- colnames( X )

        ## testing
        beta1.hat   <- beta.hat[2:K]
        cov.beta1   <- cov.beta[2:K, 2:K]

        chisq.stat  <- t(beta1.hat)%*%solve(cov.beta1, beta1.hat)

        chisq.pv    <- pchisq(chisq.stat,   df = K-1, lower.tail = FALSE)

        beta.hat <- as.numeric( beta.hat )
        names( beta.hat ) <- colnames( X )

    if ( method=="2SLS" ) {
        gamma.hat <- as.numeric( gamma.hat )
        names( gamma.hat ) <- colnames( X )
        gamma1.hat <- gamma0.hat <- NA
    } else {
        gamma1.hat <- as.numeric( gamma1.hat )
        gamma0.hat <- as.numeric( gamma0.hat )
        names( gamma0.hat ) <- names( gamma1.hat ) <- colnames( X )
        gamma.hat <- NA
    }

    ## return
    res <- list(Yobs       = Yobs,
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

    res$Y <- Yobs
    res$Z <- Z
    res$control.vars <- NA
    res$main.vars <- main.vars
    res$interaction.vars <- interaction.vars

    res$D <- D
    res$X <- X

    res$method <- paste( "LATE RI-", method, sep="" )

    class( res ) <- c("RI.regression.result", "RI.regression.result.LATE")

    return(res)
}
