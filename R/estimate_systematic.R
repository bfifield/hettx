##
## Library of functions for the estimands, etc. from the JASA paper
##
## (c) 2016 - Miratrix, Deng, Feller
##

scat <- function( str, ... ) {
    cat( sprintf( str, ... ) )
}

#' Calculate systematic effects model using LATE, ITT, or full potential outcomes.
#'
#' Implements the systematic effects model proposed in Ding, Feller, and Miratrix
#' (2018). Can estimate an ITT or LATE model, or the actual beta in cases where
#' full potential outcomes schedule is available.
#'
#' The OLS method differs from the RI method only by how the Sxx matrix is handled.
#' In the OLS case, seperate Sxx for treatment and control are calculated for each 
#' treatment arm. For RI the known Sxx based on all units is used.
#' 
#' @usage estimate.systematic(formula, data, interaction.formula, control.formula,
#' method, na.rm)
#'
#' @param formula An object of class formula, as in lm(). For ITT estimation, specify as Y ~ Z with only the treatment variable on the right-hand side. For LATE estimation, specify as Y ~ D | Z with only the endogenous variable (D) and the instrument (Z) on the right-hand side separated by a vertical bar (|). For oracle estimation (where full potential outcome schedule is known), specify as Y(1) + Y(0) ~ Z with only the treatment variable on the right-hand side and the variables indicating the outcome under treatment and the outcome under control on the left-hand-side. The first variable on the left-hand-side will be treated as the outcome under treatment, and the second variable on the right-hand-side will be treated as the outcome under control.
#' @param data A data.frame, tbl_df, or data.table with the input data.
#' @param interaction.formula A right-sided formula with pre-treatment covariates to model treatment effects for on the right hand side, such as ~ x1 + x2 + x3. 
#' @param control.formula A right-sided formula with pre-treatment covariates to adjust for on the right hand side, such as ~ x1 + x2 + x3. Default is NULL (no variables adjusted for). Will be ignored for LATE estimation and oracle estimation. Default is NULL
#' @param method RI or OLS (for ITT and oracle), RI or 2SLS (for LATE). method=OLS is shorthand for setting the empirical.Sxx variable to TRUE, nothing more.
#'
#' @param na.rm A logical flag indicating whether to list-wise delete missing data. The function will report an error if missing data exist. Default is FALSE.
#'
#' @export
#' @importFrom stats as.formula
estimate.systematic <- function( formula, data, interaction.formula, control.formula=NULL,
                    method=c("RI","OLS","2SLS"), na.rm = FALSE) {

    ## Check formula
    formula.char <- paste( deparse(formula), collapse=" " )
    method <- match.arg(method)
    if(length(lhs.vars(formula)) == 2){
        if(!is.null(control.formula)){
            cat("control.formula specified, ignoring for oracle estimation.\n")
        }
        eb.out <- calc.beta.oracle(formula=formula, data=data,
                                   interaction.formula=interaction.formula,
                                   method=method, na.rm=na.rm)
    }else if(grepl("\\|", formula.char)){
        if(!is.null(control.formula)){
            cat("control.formula specified, ignoring for LATE estimation.\n")
        }
        eb.out <- est.beta.LATE(formula=formula, data=data,
                                interaction.formula=interaction.formula,
                                method=method, na.rm=na.rm)
    }else if(length(lhs.vars(formula)) == 1 & length(rhs.vars(formula)) == 1){
        eb.out <- est.beta.ITT(formula=formula, data=data,
                               interaction.formula=interaction.formula,
                               control.formula=control.formula,
                               method=method,
                               na.rm=na.rm)
    }else{
        stop("Please input a proper argument for formula. For ITT estimation, use Y ~ Z, for LATE estimation, use Y ~ T | Z, and for oracle estimation, use Y1 + Y0 ~ Z.")
    }
    eb.out$call <- match.call()
    return(eb.out)
    
}

#' Bounds the R2 measure using OLS output for ITT from est.beta, and bounds R2 measure for compliers using LATE estimation from est.beta. 
#'
#' @usage R2(est.beta, rho.step)
#'
#' @param est.beta The output from `est.beta()`. Either an estimate of overall systematic effect variation, or systematic effect variation for compliers.
#' @param rho.step Grid size for sensitivity analysis on values of rho. Default is 0.05
#'
#' @export
#'
#' @return RI.R2.result object.
#' @seealso print.RI.R2.result
R2 <- function( est.beta, rho.step=0.05 ) {
    stopifnot( is.RI.regression.result(est.beta) )
    if( inherits(est.beta, "RI.regression.result.LATE") ){
        r2 <- R2.LATE(RI.result=est.beta, rho.step=rho.step)
    }else{
        r2 <- R2.ITT(RI.result=est.beta, rho.step=rho.step)
    }
    return(r2)
}

## ------------------------------------
## Clean printing of our result objects
## ------------------------------------
is.RI.regression.result <- function( x ) {
    inherits(x, "RI.regression.result")
}

#' Extract coefficients of a fit RI regression model.
#' @param object A RI.regression.result object.
#' @param ... Unused
#'
#' @export
coef.RI.regression.result <- function( object, ... ) {
    object$beta.hat
}

#' Get vcov() from object.
#'
#' @param object est.beta object
#' @param ... unused
#'
#' @export
vcov.RI.regression.result <- function( object, ... ) {
    object$cov.beta
}

#' Extract the standard errors from a var-cov matrix.
#'
#' @param object est.beta object
#' @param ... unused
#'
#' @export
SE <- function( object, ... ) {
    sqrt( diag( vcov( object ) ) )
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
variance.ratio.test <- function(Yobs, Z, data= NULL)
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




