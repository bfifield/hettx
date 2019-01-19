#' @export
summary.FRTCI.test <- function(object, ...){

    ## Create data frame
    rng <- range(object$ci.p)
    df <- data.frame(object$statistic, object$p.value, object$p.value.plug, rng[1], rng[2])
    names(df) <- c("Estimate", "P-Value (Sweep)", "P-Value (Plug-In)",
                   "P-Value Lower CI", "P-Value Upper CI")
    rownames(df) <- NULL
    
    ## Create output
    out <- vector(mode = "list")
    out$call <- object$call
    out$estimates <- df
    out$test.stat <- object$test.stat
    out$B <- object$B
    out$gamma <- object$gamma
    class(out) <- "summary.FRTCI.test"
    return(out)   
}

#' @export
print.summary.FRTCI.test <- function(x, ...){
    
    cat("\n")
    cat("Call:\n")
    print(x$call)
    cat("\n")
    cat("Test Statistic:", x$test.stat, "\n")
    cat("\n")
    cat("Estimates:\n")
    print(x$estimates)
    
}

#' @export
print.FRTCI.test = function( x, ... ) {
    print(summary(x)$estimates)
}
 
#' get p-value along with uncertainty on p-value
#' 
#' Give confidence bounds (from monte carlo simulation error) for the p-values returned by a test
#' 
#' @param tst A FRTCI.test object from fishpidetect
#' @return p-value and range of p-values due to monte carlo error.
#' @export
get.p.value <- function( tst ) {
    cnts = (tst$ci.p - tst$gamma) * tst$B
    bts = sapply( cnts, function( cnt ) {
        bt = binom.test( cnt, tst$B )
        as.numeric( bt$conf.int )
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
