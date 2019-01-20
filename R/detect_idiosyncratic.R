## Set global variables for foreach
globalVariables('b')

#' detect.idiosyncratic
#'
#' Test for systematic treatment effect heterogeneity using
#' Fisherian permutation inference methods.
#' @usage detect.idiosyncratic(formula, data, interaction.formula, control.formula,
#' plugin, tau.hat, test.stat, te.vec, B, gamma, grid.gamma, grid.size,
#' return.matrix, na.rm, n.cores, verbose, ...)
#'
#' @param formula An object of class formula, as in lm(), such as Y ~ Z with only the treatment variable on the right-hand side.
#' @param data A data.frame, tbl_df, or data.table with the input data.
#' @param interaction.formula A right-sided formula with pre-treatment covariates to model treatment effects for on the right hand side, such as ~ x1 + x2 + x3. Defaultis NULL (no interactions modeled)
#' @param control.formula A right-sided formula with pre-treatment covariates to adjust for on the right hand side, such as ~ x1 + x2 + x3. Default is NULL (no variables adjusted for)
#' @param plugin Whether to calculate the plug-in p-value without sweeping over range of possible treatment effect magnitudes. Default is FALSE.
#' @param tau.hat The value of the plug-in treatment effect. Default is sample average treatment effect.
#' @param test.stat  Test statistic function to use on the data. Default is shifted Kolmogorov-Smirnov statistic.
#' @param te.vec Vector of taus to examine if you want to override generating ones automatically. Default is NULL.
#' @param B  Number of permutations to take. Default is 500.
#' @param gamma How wide of a CI to make around tau-hat for search. Default is 0.0001.
#' @param grid.gamma Parameter to govern where the grid points are sampled.  Bigger values means more samples towards the estimated tau-hat. Default is 100*gamma.
#' @param grid.size Number of points in the grid. Default is 151.
#' @param return.matrix  Whether to return the matrix of all the imputed statistics.  Default is FALSE.
#' @param na.rm A logical flag indicating whether to list-wise delete missing data. The function will report an error if missing data exist. Default is FALSE.
#' @param n.cores Number of cores to use to parallelize permutation step. Default is 1.
#' @param verbose  Whether to print out progress bar when fitting and other diagnostics. Default is TRUE.
#' @param ... Extra arguments passed to the generate.permutations function and test.stat functions.
#'
#' @return If plug-in, the value of the test and the associated p-value. If not, a list with the value of
#' the test statistic on the observed data, the value of the CI-adjusted p-value, the plug-in p-value, and other information on the test.
#'
#' @examples
#' Z <- rep(c(0, 1), 100)
#' tau <- 4
#' Y <- ifelse(Z, rnorm(100, tau), rnorm(100, 0))
#' df <- data.frame(Y=Y, Z=Z)
#' tst <- detect.idiosyncratic(Y ~ Z, df, B = 50, grid.size = 50)
#' @export
#' @importFrom graphics abline lines plot rug
#' @importFrom stats binom.test binomial coef confint ecdf glm lm lowess predict qchisq qnorm var vcov na.omit
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom quantreg rq
#' @importFrom plyr ddply summarize .
#' @importFrom mvtnorm rmvnorm
#' @importFrom foreach "%do%" "%dopar%" foreach
#' @importFrom parallel makePSOCKcluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom formula.tools get.vars lhs.vars rhs.vars
detect.idiosyncratic <- function(formula, data,
                                 interaction.formula = NULL, control.formula = NULL,
                                 plugin = FALSE, tau.hat = NULL,
                                 test.stat = ifelse( is.null(W) & is.null(X), "SKS.stat", ifelse( is.null(W), "SKS.stat.cov", "SKS.stat.int.cov" ) ),
                                 te.vec = NULL,
                                 B = 500, gamma = 0.0001, grid.gamma = 100*gamma,
                                 grid.size = 151, return.matrix = FALSE, na.rm = FALSE,
                                 n.cores = 1, verbose = TRUE, ...){
    
    ## ---------------------------
    ## Get variables from formulas
    ## and check formulas
    ## ---------------------------
    if(any(class(data) %in% c("tbl_df", "data.table", "data.frame"))) {
        data <- as.data.frame(data)
    }else{
        stop("The input data must be of class tbl_df, data.table, or data.frame.")
    }
    
    if(length(lhs.vars(formula)) != 1 | length(rhs.vars(formula)) != 1){
        stop("The formula argument must be of the form outcome ~ treatment.")
    }
    main.vars <- get.vars(formula)
    if(any(!(main.vars %in% colnames(data)))){
        stop("Some variables in formula are not present in your data.")
    }
    
    if(!is.null(interaction.formula)){
        if(length(lhs.vars(interaction.formula)) != 0){
            stop("Do not provide an outcome variable in interaction.formula.")
        }
        interaction.vars <- get.vars(interaction.formula)
        if(any(!(interaction.vars %in% colnames(data)))){
            stop("Some variables in interaction.formula are not present in your data.")
        }
        if(any(main.vars %in% interaction.vars)){
            stop("Either your outcome variable or your treatment variable is present in your interaction formula.")
        }
    }else{
        interaction.vars <- NULL
    }

    if(!is.null(control.formula)){
        if(length(lhs.vars(control.formula)) != 0){
            stop("Do not provide an outcome variable in control.formula.")
        }
        control.vars <- get.vars(control.formula)
        if(any(!(control.vars %in% colnames(data)))){
            stop("Some variables in control.formula are not present in your data.")
        }
        if(any(main.vars %in% interaction.vars)){
            stop("Either your outcome variable or your treatment variable is present in your control formula.")
        }
        if(!is.null(interaction.formula)){
            if(any(control.vars %in% interaction.vars)){
                stop("You have variables in your interaction formula that are also present in your control formula.")
            }
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
    Y <- data[,main.vars[1]]
    Z <- data[,main.vars[2]]
    if(!is.null(interaction.formula)){
        W <- as.matrix(data[,interaction.vars])
    }else{
        W <- NULL
    }
    if(!is.null(control.formula)){
        X <- as.matrix(data[,control.vars])
    }else{
        X <- NULL
    }
    
    ## ------------------------
    ## Run checks on input data
    ## ------------------------
    ## Class checks on Y and Z
    if(!inherits(Y, "numeric")){
        stop("Y must be a numeric vector.")
    }
    if(inherits(Z, "logical")){
        Z <- ifelse(Z, 1, 0)
    }else if(inherits(Z, "factor") | inherits(Z, "character")){
        Z <- as.character(Z)
        tcat <- unique(Z)[1]
        Z <- ifelse(Z == tcat, 1, 0)
        cat("Z is factor or character class, setting", tcat, "as treatment category.\n")
    }

    ## Class and input checks on W and X
    if(plugin & !is.null(interaction.formula)){
        warning("Running plug-in test but covariates provided for interaction.formula, will not adjust for specified heterogeneity.")
    }
    if(plugin & is.null(tau.hat)){
        warning("Using plug-in test with no argument passed to tau.hat, will estimate tau.hat from the data.\n")
    }
    if(!is.null(te.vec) & !is.null(W)){
        warning("Adjusting for existing heterogeneity but treatment effect vector provided, ignoring provided treatment effect vector.")
    }
    if(!is.null(W)){
        if(test.stat %in% c("WSKS.t", "SKS.pool.t")){
            if(ncol(W) > 1){
                stop("Can only adjust for a single (categorical) covariate when using WSKS.t or SKS.pool.t as test statistics.")
            }
            W <- as.factor(W)
        }
        if(!inherits(W, "factor")){
            unq.check <- apply(W, 2, function(x){length(unique(x))})
        }else{
            unq.check <- length(unique(W))
        }
        if(any(unq.check == 1)){
            stop("Some variables specified for the interaction model in interaction.formula have no variation.")
        }
    }
    if(!is.null(X)){
        if(!(inherits(X, "data.frame") | inherits(X, "matrix"))){
            stop("X must be either a data frame or a matrix.")
        }
        unq.check <- apply(X, 2, function(x){length(unique(x))})
        if(any(unq.check == 1)){
            stop("Some variables specified for adjustment in control.formula have no variation.")
        }
    }

    ## Checks on functions
    no.adj.funs <- c("KS.stat", "SKS.stat", "rq.stat")
    adj.funs <- c("SKS.stat.cov.pool", "SKS.stat.cov", "SKS.stat.cov.rq", "rq.stat.cond.cov", "rq.stat.uncond.cov")
    adj.int.funs <- c("SKS.stat.int.cov.pool", "SKS.stat.int.cov")
    int.funs <- c("WSKS.t", "SKS.pool.t")
    if(is.null(W) & is.null(X)){
        if(!(test.stat %in% no.adj.funs)){
            stop("You have provided an invalid test statistic when not adjusting for covariates or specifying interactions. Must provide one of KS.stat, SKS.stat, or rq.stat. Please see test.stat.info() for more information.")
        }
    }else if(is.null(W)){
        if(!(test.stat %in% adj.funs)){
            stop("You have provided an invalid test statistic when adjusting for covariates in control.formula but not specifying interactions. Must provide one of SKS.stat.cov.pool, SKS.stat.cov, SKS.stat.cov.rq, rq.stat.cond.cov, or rq.stat.uncond.cov. Please see test.stat.info() for more information.")
        }
    }else if(is.null(X)){
        just.int.funs <- c(adj.int.funs, int.funs)
        if(!(test.stat %in% just.int.funs)){
            stop("You have provided an invalid test statistic when specifying interactions in interaction.formula but not adjusting for covariates. Must provide one of SKS.stat.int.cov.pool, SKS.stat.int.cov, WSKS.t, or SKS.pool.t. Please see test.stat.info() for more information.")
        }
    }else{
        if(!(test.stat %in% adj.int.funs)){
            stop("You have provided an invalid test statistic when specifying interactions in interaction.formula and adjusting for covariates in control.formula. Must provide one of SKS.stat.int.cov.pool or SKS.stat.int.cov. Please see test.stat.info() for more information.")
        }
    }
    
    ## ------------------------------------------------
    ## Detect whether to run plugin, FRTCI, or interact
    ## ------------------------------------------------
    if(plugin){ ## Plug-in test
        if(is.null(tau.hat)){
            tau.hat <- mean(Y[Z == 1]) - mean(Y[Z == 0])
        }
        fpi_out <- FRTplug(Y = Y, Z = Z, test.stat = match.fun(test.stat), tau.hat = tau.hat, verbose= verbose, ...)
    }else if(is.null(W)){ ## FRTCI - with or without adjusting
        fpi_out <- FRTCI(Y = Y, Z = Z, X = X, test.stat = match.fun(test.stat), B = B, gamma = gamma,
                         grid.gamma = grid.gamma, grid.size = grid.size,
                         te.vec = te.vec, return.matrix = return.matrix,
                         n.cores = n.cores, verbose = verbose, ...)
    }else{ ## FRTCI.interact
        fpi_out <- FRTCI.interact(Y = Y, Z = Z, W = W, X = X,
                                  test.stat = match.fun(test.stat), B = B, gamma = gamma,
                                  grid.gamma = grid.gamma, grid.size = grid.size,
                                  return.matrix = return.matrix, n.cores = n.cores,
                                  verbose = verbose, ...)
    }

    ## Add to output
    fpi_out$call <- match.call()
    fpi_out$test.stat <- test.stat

    return(fpi_out)

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

