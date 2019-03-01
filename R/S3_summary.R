#' @export
summary.FRTCI.test <- function(object, ...){

    ## Create data frame
    df <- data.frame(object$statistic, object$p.value, object$p.value.plug)
    names(df) <- c("Statistic", "P-Value (Sweep)", "P-Value (Plug-In)" )
    rownames(df) <- NULL
    
    ## Create output
    out <- vector(mode = "list")
    out$call <- object$call
    out$estimates <- df
    out$test.stat <- object$test.stat
    out$B <- object$B
    out$gamma <- object$gamma
    out$method = object$method
    out$FRTCI = object
    
    if ( !is.null( object$W ) ) {
      out$grid.range = apply( object$te.grid, 2, range )
    } else {
    
    cnts = (object$ci.p - object$gamma) * object$B
    bts = sapply( cnts, function( cnt ) {
      bt = binom.test( cnt, object$B )
      bt$conf.int
    } )
    minp = min( bts[1,] )
    maxp = max( bts[2,] )
    out$p.value.CI = c( minp, maxp ) + object$gamma 
    }
    
    
    class(out) <- "summary.FRTCI.test"
    return(out)   
}

#' @export
summary.RI.regression.result <- function(object, ...){
    out <- vector(mode = "list")
    out$method <- object$method
    out$call <- object$call
    if(length(coef(object))){
        out$coefficients <- coef(object)
    }
    out$vcov <- vcov(object)
    out$chisq.stat <- object$chisq.stat
    out$p.value <- object$p.value
    out$ATE <- object$ATE
    out$SE.ATE <- object$SE.ATE
    out$SD.Y0 <- object$SD.Y0
    out$SD.Y1 <- object$SD.Y1
    class(out) <- "summary.RI.regression.result"

    return(out)
}

#' @export
summary.RI.R2.result <- function(object, ...){

    ## Type
    out <- vector(mode = "list")
    out$method <- object$type

    if(object$type == "ITT"){

        ## -------
        ## For ITT
        ## -------
        
        ## Data frame for R2
        df_hte_r2 <- data.frame(
            object$R2.lower, object$R2.lower.sharp, object$R2.upper
        )
        names(df_hte_r2) <- c(
            "R2 Lower Bound", "R2 Lower Bound (Sharp)", "R2 Upper Bound"
        )
        rownames(df_hte_r2) <- NULL
        out$hte_r2 <- df_hte_r2
        
        ## Store systematic treatment effect variance
        out$hte_variance_systematic <- object$Sdd
        
        ## Data frame for idiosyncratic variance
        df_hte_idio <- data.frame(
            object$See.lower, object$See.upper.sharp, object$See.upper
        )
        names(df_hte_idio) <- c(
            "Lower Bound", "Upper Bound (Sharp)", "Upper Bound"
        )
        rownames(df_hte_idio) <- NULL
        out$hte_variance_idiosyncratic <- df_hte_idio

        ## Data frame for total variance
        df_hte_total <- data.frame(
            object$Sdd + object$See.lower, object$Sdd + object$See.upper.sharp,
            object$Sdd + object$See.upper
        )
        names(df_hte_total) <- c(
            "Lower Bound", "Upper Bound (Sharp)", "Upper Bound"
        )
        out$hte_variance_total <- df_hte_total
    }else{

        ## --------
        ## For LATE
        ## --------

        ## Data frame for R2
        df_hte_r2 <- data.frame(
            c(object$R2.lower, object$R2.U.lower, object$R2.UX.lower),
            c(object$R2.lower.sharp, object$R2.U.lower.sharp,
              object$R2.UX.lower.sharp),
            c(object$R2.upper, object$R2.U.upper, object$R2.UX.upper)
        )
        names(df_hte_r2) <- c(
            "R2 Lower Bound", "R2 Lower Bound (Sharp)", "R2 Upper Bound"
        )
        rownames(df_hte_r2) <- c(
            "Compliers", "Noncompliers", "Covariates and compliers"
        )
        out$hte_r2 <- df_hte_r2

        ## Variances
        out$hte_variance_systematic_compliers <- object$Sdd
        out$hte_variance_systematic_strata <- object$Stautau.U
        totsys = object$pi.c * object$Sdd + object$Stautau.U
        out$hte_variance_systematic_total <- totsys

        ## Variance Data frames
        df_hte_idio <- data.frame(
            object$See.lower, object$See.upper.sharp, object$See.upper
        )
        names(df_hte_idio) <- c(
            "Lower Bound", "Upper Bound (Sharp)", "Upper Bound"
        )
        rownames(df_hte_idio) <- NULL
        out$hte_variance_idiosyncratic <- df_hte_idio

        df_hte_total <- data.frame(
            totsys + object$pi.c * object$See.lower,
            totsys + object$pi.c * object$See.upper.sharp,
            totsys + object$pi.c * object$See.upper
        )
        names(df_hte_total) <- c(
            "Lower Bound", "Upper Bound (Sharp)", "Upper Bound"
        )
        out$hte_variance_total <- df_hte_total

        out$LATE <- object$tau.c
        out$ITT <- object$ITT
        out$prop_compliers <- object$pi.c
        
    }
    
    class(out) <- "summary.RI.R2.result"
    return(out)
    
}

