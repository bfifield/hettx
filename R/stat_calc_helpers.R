#' test.stat.info
#'
#' A list of test statistics for detect.idiosyncratic(), and information on use cases when each is appropriate.
#'
#' @usage test.stat.info()
#'
#' @export
test.stat.info <- function(){
    cat("
## -----------------------------------------------------
## Available test statistics for detect.idiosyncratic():
## -----------------------------------------------------
Test statistics when not adjusting for covariates or specifying interactions:
\t1) KS.stat - Calculate classic (not shifted) KS statistic. If tau passed, Y1 will be shifted by tau.
\t2) SKS.stat - Shifted KS statistic. Calculate KS distance between Y0 and Y1 shfited by tau.
\t3) rq.stat - Classic KS statistic via quantile regression with covariates\n
Test statistics when adjusting for covariates but not specifying interactions:
\t1) SKS.stat.cov.pool - Shifted KS statistic with covariates to increase precision.
\t2) SKS.stat.cov - Shifted KS statistic with covariates with model for outcomes calculated on control group only.
\t   This avoids 'splitting' the treatment variation between treatment and control groups. We recommend this method
\t   over the 'pool' method in SKS.stat.cov.pool.
\t3) SKS.stat.cov.rq - Shifted KS statistic via quantile regression with covariates.
\t4) rq.stat.cond.cov - KS statistic via quantile regression with covariates. Conditional approach; see Koenker
\t    and Xiao (2002) for more information.
\t5) rq.stat.uncond.cov - KS statistic via quantile regression with covariates. Unconditional approach; see
\t    Firpo (2007) for more information.\n
Test statistics when specifying interactions, with or without covariate adjustment:
\t1) SKS.stat.int.cov.pool - Shifted KS statistic with a linear treatment effect model and optional covariate
\t   adjustment. This will attempt to remove any systematic variation corresponding to the specified interaction
\t   model and then return an SKS statistic on the residuals to measure any variation 'left over'. This is the
\t   test statistic used in Ding, Feller, and Miratrix (2016), JRSS-B.
\t2) SKS.stat.int.cov - Similar to SKS.stat.int.cov.pool, but this method first adjusts for baseline and then
\t   models treatment effects on the residuals to not split treatment effects. We recommend this method over
\t   the 'pool' method in SKS.stat.int.cov.pool.\n
Test statistics when specifying interactions, without covariate adjustment:
\t1) WSKS.t - Calculates the shifted KS statistic within each group specified in the interaction model, and
\t   then aggregates together as a weighted average. Should be used when the interaction model is a single
\t   categorical or factor covariate.
\t2) SKS.pool.t - Subtract off group level treatment effect estimates and then look at KS statistic on residuals.
\t   Should be used when the interaction model is a single categorical or factor covariate.\n")}

#' KS.stat
#'
#' Calculate classic (not shifted) KS statistic, code modified version of R's ks.test.
#' If tau passed, Y1 will be shifted by tau.
#'
#' @usage KS.stat(Y, Z, tau, alternative)
#'
#' @param Y Observed outcome vector
#' @param Z Treatment assigment vector
#' @param tau Value of treatment effect for shifting Y1. Default is NULL (Y1 not shifted).
#' @param alternative Direction of test ("two.sided", "less", "greater")
#'
#' @return The value of the test.
#' @export
KS.stat <- function( Y, Z, tau = NULL, alternative = c("two.sided", "less", "greater") ) {
    x = Y[Z==1]
    y = Y[Z==0]
    if ( !is.null( tau ) ) {
        x = x - tau
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

#' SKS.stat
#' 
#' Shifted kolmogorov-smirnov statistic. Calculate KS distance between Y0 and Y1
#' shifted by sample tau.
#'
#' @usage SKS.stat(Y, Z)
#'
#' @param Y Observed outcome vector
#' @param Z Treatment assigment vector
#'
#' @return The value of the test.
#'
#' @export
SKS.stat <- function(Y, Z)
{
    Y1 = Y[Z==1]
    Y0 = Y[Z==0]
    
    Y1.star   = Y1 - mean(Y1)
    Y0.star   = Y0 - mean(Y0)
    
    unique.points = c(Y1.star, Y0.star)
    
    Fn1 = ecdf(Y1.star)
    Fn0 = ecdf(Y0.star)
    
    difference = Fn1(unique.points) - Fn0(unique.points)
    
    return(max(abs(difference)))
    
}

#' SKS.stat.cov.pool
#' 
#' Shifted kolmogorov-smirnov statistic with covariates
#' to increase precision.
#'
#' This is the test statistic used Ding, Feller, and Miratrix (2016), JRSS-B.
#'
#' @usage SKS.stat.cov.pool(Y, Z, X)
#'
#' @param Y  Observed outcome vector
#' @param Z  Treatment assigment vector 
#' @param X Additional pre-treatment covariates to adjust for in estimation, but not to interact with treatment. Default is NULL.
#'
#' @return The value of the test.
#'
#' @export
SKS.stat.cov.pool <- function(Y, Z, X)
{
    this.lm <- lm(Y ~ Z + X)
    Y1.star <- this.lm$res[Z == 1]
    Y0.star <- this.lm$res[Z == 0]
    
    unique.points = c(Y1.star, Y0.star)
    
    Fn1 = ecdf(Y1.star)
    Fn0 = ecdf(Y0.star)
    
    difference = Fn1(unique.points) - Fn0(unique.points)
    
    return(max(abs(difference)))
}

#' SKS.stat.cov
#' 
#' Shifted kolmogorov-smirnov statistic with covariates
#' with model for outcomes calculated on control group only.
#' This avoids "splitting" the treatment variation between tx 
#' and co groups.  
#' We recommend this method over the "pool" method.
#'
#' @usage SKS.stat.cov(Y, Z, X)
#'
#' @inheritParams SKS.stat.cov.pool
#'
#' @return The value of the test.
#'
#' @export
SKS.stat.cov <- function(Y, Z, X)
{
    this.lm <- lm(Y ~ X, subset = Z == 0)
    
    Y0.hat <- predict(this.lm, newdata = as.data.frame(X))
    Y0.res <- Y - Y0.hat
    
    Y1.star <- Y0.res[Z == 1] - mean(Y0.res[Z == 1])
    Y0.star <- Y0.res[Z == 0] - mean(Y0.res[Z == 0])
    
    unique.points = c(Y1.star, Y0.star)
    
    Fn1 = ecdf(Y1.star)
    Fn0 = ecdf(Y0.star)
    
    difference = Fn1(unique.points) - Fn0(unique.points)
    
    return(max(abs(difference)))
}

#' SKS.stat.int.cov.pool
#'
#' Shifted kolmogorov-smirnov statistic with a linear treatment
#' effect model defined by W
#'
#' This will attempt to remove any systematic variation corresponding
#' to W and then return a SKS statistic on the residuals to measure
#' any variation "left over".
#'
#' X are _additional_ covariates to adjust for beyond those involved
#' in treatment effect model.  It will automatically adjust for W as
#' well.  Do not put a covariate in for both X and W.
#'
#' This is the test statistic used in Ding, Feller, and Miratrix (2016), JRSS-B.
#'
#' @usage SKS.stat.int.cov.pool(Y, Z, W, X)
#'
#' @param Y  Observed outcome vector
#' @param Z  Treatment assigment vector
#' @param W Additional pre-treatment covariates to interact with T to define linear model of treatment effects. 
#' @param X Additional pre-treatment covariates to adjust for in estimation, but not to interact with treatment. Default is NULL.
#'
#' @return The value of the test.
#'
#' @export
SKS.stat.int.cov.pool <- function( Y, Z, W, X=NULL )
{ 
    if ( !is.null( X ) ) {
        this.lm <- lm( Y ~ Z + X + W + Z:W )
    } else {
        this.lm <- lm( Y ~ Z + W + Z:W )        
    }
    Y1.star <- this.lm$res[Z == 1]
    Y0.star <- this.lm$res[Z == 0]
    
    unique.points = c(Y1.star, Y0.star)
    
    Fn1 = ecdf(Y1.star)
    Fn0 = ecdf(Y0.star)
    
    difference = Fn1(unique.points) - Fn0(unique.points)
    
    return(max(abs(difference)))
}

#' SKS.stat.int.cov
#' 
#' Shifted kolmogorov-smirnov statistic with a linear treatment
#' effect model defined by W
#'
#' This will attempt to remove any systematic variation corresponding
#' to W and then return a SKS statistic on the residuals to measure
#' any variation "left over".
#'
#' X are _additional_ covariates to adjust for beyond those involved
#' in treatment effect model.  It will automatically ajust for W as
#' well.  Do not put a covariate in for both X and W.
#'
#' This is the test statistic used in Ding, Feller, and Miratrix (2016), JRSS-B.
#'
#' This method first adjusts for baseline and then models treatment effect
#' on the residuals to not split treatment effects.
#'
#' We recommend this method over the "pool" method.
#'
#' @usage SKS.stat.int.cov(Y, Z, W, X)
#'
#' @inheritParams SKS.stat.int.cov.pool
#'
#' @return The value of the test.
#'
#' @export
SKS.stat.int.cov <- function( Y, Z, W, X=NULL )
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
    
    unique.points = c(Y1.star, Y0.star)
    
    Fn1 = ecdf(Y1.star)
    Fn0 = ecdf(Y0.star)
    
    difference = Fn1(unique.points) - Fn0(unique.points)
    
    return(max(abs(difference)))
}



##
## Other possible test statistics
##

#' SKS.stat.cov.rq
#' 
#' Shifted kolmogorov-smirnov statistic with covariates and rq.
#'
#' @usage SKS.stat.cov.rq(Y, Z, X)
#'
#' @inheritParams SKS.stat.cov.pool
#'
#' @return The value of the test.
#' @export
SKS.stat.cov.rq <- function(Y, Z, X)
{

    this.lm <- lm(Y ~ Z + X)
    this.rq <- suppressWarnings(rq(Y ~ Z + X, tau = seq(0.05, 0.95, by = 0.05)))
    
    z.rq <- coef(this.rq)["Z",]
    return( max(z.rq - coef(this.lm)["Z"]) )
    
}

#' rq.stat
#' 
#' Kolmogorov-smirnov statistic via quantile regression with covariates
#'
#' @usage rq.stat(Y, Z, rq.pts)
#'
#' @param Y  Observed outcome vector
#' @param Z  Treatment assigment vector
#' @param rq.pts Sequence of quantile points at which to evaluate the test. Default is seq(.1, .9, by = .1). Should not go beyond 0 and 1.
#'
#' @return The value of the test.
#' @export
rq.stat <- function(Y, Z, rq.pts = seq(0.1, 0.9, by = 0.1))
{

    if(min(rq.pts) <= 0 | max(rq.pts) >= 1){
        stop("All values of rq.pts must be strictly greater than 0 and strictly less than 1.")
    }
    
    this.lm <- lm(Y ~ Z)
    this.rq <- suppressWarnings(rq(Y ~ Z, tau = rq.pts))
    
    z.rq <- coef(this.rq)["Z",]
    return( max(z.rq - coef(this.lm)["Z"]) )
    
}

#' rq.stat.cond.cov
#' 
#' Kolmogorov-smirnov statistic via quantile regression with covariates
#' Conditional approach; see Koenker and Xiao (2002)
#'
#' @usage rq.stat.cond.cov(Y, Z, X, rq.pts)
#'
#' @param Y  Observed outcome vector
#' @param Z  Treatment assigment vector
#' @param X Additional pre-treatment covariates to adjust for in estimation, but not to interact with treatment. 
#' @param rq.pts Sequence of quantile points at which to evaluate the test. Default is seq(.1, .9, by = .1). Should not go beyond 0 and 1.
#'
#' @return The value of the test.
#' @export
rq.stat.cond.cov <- function(Y, Z, X, rq.pts = seq(0.1, 0.9, by = 0.1))
{

    if(min(rq.pts) <= 0 | max(rq.pts) >= 1){
        stop("All values of rq.pts must be strictly greater than 0 and strictly less than 1.")
    }
    
    this.lm <- lm(Y ~ Z + X)
    this.rq <- suppressWarnings(rq(Y ~ Z + X, tau = rq.pts))
    
    z.rq <- coef(this.rq)["Z",]
    return( max(z.rq - coef(this.lm)["Z"]) )
    
}


#' rq.stat.uncond.cov
#' 
#' Kolmogorov-smirnov statistic via quantile regression with covariates
#' Unconditional approach; see Firpo (2007)
#'
#' @usage rq.stat.uncond.cov(Y, Z, X, rq.pts)
#'
#' @inheritParams rq.stat.cond.cov
#'
#' @return The value of the test.
#' @export
rq.stat.uncond.cov <- function(Y, Z, X, rq.pts = seq(0.1, 0.9, by = 0.1))
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

#' WSKS.t
#'
#' Weighted average of the group-level SKS statistics
#'
#' @usage WSKS.t(Y, Z, W)
#'
#' @param Y Observed outcome vector
#' @param Z Treatment assigment vector
#' @param W A a factor or categorical covariate.
#'
#' @return The value of the test.
#' @export
WSKS.t <- function( Y, Z, W ) {
    
    dd = ddply( data.frame(Y=Y,Z=Z,W=W), "W", summarize, 
               t.sks = SKS.stat( Y, Z ),
               n.k = length(Y) )
    n = length(Y)
    return( sum( dd$t.sks * dd$n.k / n ) )
}

#' SKS.pool.t
#' 
#' Subtract off group level treatment effect estimates and then look
#' at KS statistic on residuals.  
#'
#' Distinct from the interacted lm in that the control units are not
#' shifted and centered with respect to eachother.
#'
#' @usage SKS.pool.t(Y, Z, W)
#'
#' @inheritParams WSKS.t
#'
#' @return The value of the test.
#' @export
SKS.pool.t <- function( Y, Z, W ) {
    
    dat <- data.frame( Y=Y, Z=Z, W=W )
    mns = ddply( dat, .(W, Z), summarize, mean=mean(Y) )
    taus = with( mns, mean[Z==1] - mean[Z==0] )
    return( KS.stat( dat$Y - dat$Z*taus[dat$W], dat$Z ) )
}
