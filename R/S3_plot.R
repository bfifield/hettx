#' plot.FRTCI.test
#' 
#' Plot curve from FRTCI.test object.
#'
#' @usage \method{plot}{FRTCI.test}(x, true.tau, xlab, ylab, true.tau.col, plot.envelope, ci.line.col, ...)
#'
#' @param x An object of class \code{FRTCI.test}
#' @param true.tau The true value of tau, if known. Default is NULL.
#' @param xlab X-axis label. Default is tau.
#' @param ylab Y-axis label. Default is "p-value".
#' @param true.tau.col Color to plot true tau value, if provided. Default is red.
#' @param plot.envelope Plot envelope around tested values of tau. Default is TRUE.
#' @param ci.line.col Color to plot confidence interval around estimated treatment effect. Default is blue.
#' @param ... Further arguments to be passed to \code{print.FRTCI.test()}
#'
#' @examples
#' Z <- rep(c(0, 1), 100)
#' tau <- 4
#' Y <- ifelse(Z, rnorm(100, tau), rnorm(100, 0))
#' df <- data.frame(Y=Y, Z=Z)
#' tst <- detect_idiosyncratic(Y ~ Z, df, B = 50, grid.size = 50)
#' plot(tst)
#' 
#' @export
#' @method plot FRTCI.test
plot.FRTCI.test <- function( x, true.tau=NULL, 
                             xlab=expression(tau), ylab="p-value", true.tau.col="red",
                             plot.envelope=TRUE, ci.line.col="blue", ... ) {
    
  if ( !is.null( x$W ) ) {
    stop( "No default plot for test results beyond a systematic model" )
  }
  
    cnts <- (x$ci.p - x$gamma) * x$B
    bts <- vapply( cnts, function( cnt ) {
        bt <- binom.test( cnt, x$B )
        as.numeric( bt$conf.int )
    }, numeric(2) )
    
    plot( x$te.vec, x$ci.p, ylim=range(bts), type="l", xlab=xlab, ylab=ylab, ...  )
    abline( v=x$te.hat, col= ci.line.col )
    abline( v=c(x$te.hat-x$te.MOE, x$te.hat+x$te.MOE), col= ci.line.col, lty=2 )
    
    if ( plot.envelope ) {
        lines( x$te.vec, bts[1,], lty=3, col="grey" )
        lines( x$te.vec, bts[2,], lty=3, col="grey" )
        lines( lowess(x$te.vec, x$ci.p, f=1/6),lty=3,lwd=2,col="darkgray" )
    }
    
    rug( x$te.vec, ticksize=0.01 )
    
    if (!is.null(true.tau) ) {
        abline( v=true.tau, lwd=2, col= true.tau.col )
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
#' @examples
#' df <- make.randomized.dat( 1000, gamma.vec=c(1,1,1,2), beta.vec=c(-1,-1,1,0) )
#' es <- estimate_systematic( Yobs ~ Z,  interaction.formula = ~ A + B, data = df )
#' r2_out <- R2(es)
#' plot(r2_out)
#' 
#' @export
#' @seealso calc_beta_oracle
plot.RI.R2.result <- function( x, main=paste( "R2 for Het Tx (", x$type, ")", sep=""),
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

