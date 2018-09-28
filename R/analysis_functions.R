#' print.FRTCI.test
#' 
#' Print out FRTCI.test object.
#'
#' @usage \method{print}{FRTCI.test}(x, ...)
#'
#' @param x An object of class \code{FRTCI.test}
#' @param ... Further arguments to be passed to \code{print.FRTCI.test()}
#'
#' @export
#' @method print FRTCI.test
print.FRTCI.test = function( x, ... ) {
    cat( "\t", x$method, "\n" )
    cat( "data: ", x$data.name, "\n" )
    cat( "t =", x$statistic, ", p-value =", x$p.value, " (plug = ", x$p.value.plug, ")\n" )
    rng = range(x$ci.p)
    cat( "\tp-value range =", rng[1], "-", rng[2] , "\n" )
    
    if ( !is.null( x$te.vec ) ) {    
        rng = range( x$te.vec )
        cat( "CI range = ", rng[1], "-",rng[2], "\tMOE =", x$te.MOE, "\n" )
        cat( "\tNeyman Difference point estimate of average =", x$te.hat, " SE =", x$te.SE,  "\n" )
    }
    cat( "\t# tested points =", length( x$ci.p ), "\n" )
    cat( "\tB = ", x$B, "\tgamma = ", x$gamma, "\n" )
    
    if ( !is.null( x$W ) ) {
        cat( "\tCorrected for", paste( colnames(x$W), sep=", " ) )        
    }

    if ( !is.null( x$X ) ) {
        cat( "\tAdjusted for", paste( colnames(x$X), sep=", " ) )        
    }
    
    
}

## Give confidence bounds (from monte carlo simulation error) for the p-values returned by a test
get.p.value <- function( tst ) {
    cnts = (tst$ci.p - tst$gamma) * tst$B
    bts = sapply( cnts, function( cnt ) {
        bt = binom.test( cnt, tst$B )
        confint( bt )[2:3]
    } )
    stopifnot( tst$p.value == max( tst$ci.p ) )
    c( p.value=max( tst$ci.p ), min.p= min( bts[1,] ), max.p=max( bts[2,] ), plug=tst$p.value.plug )
}

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
#' @export
#' @method plot FRTCI.test
plot.FRTCI.test <- function( x, true.tau=NULL, 
                             xlab=expression(tau), ylab="p-value", true.tau.col="red",
                             plot.envelope=TRUE, ci.line.col="blue", ... ) {
    
    cnts = (x$ci.p - x$gamma) * x$B
    bts = sapply( cnts, function( cnt ) {
        bt = binom.test( cnt, x$B )
        bt$conf.int
    } )
    
    plot( x$te.vec, x$ci.p, ylim=range(bts), type="l", xlab=xlab, ylab=ylab, ...  )
    abline( v=x$te.hat, col= ci.line.col )
    abline( v=c(x$te.hat-x$te.MOE,x$te.hat+x$te.MOE), col= ci.line.col, lty=2 )
    
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
