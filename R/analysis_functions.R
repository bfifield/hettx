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

