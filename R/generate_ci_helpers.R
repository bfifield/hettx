## Generate truncated multivariate normal draws
## where we don't return anything in the tails defined by alpha.
rcrmtvnorm <- function(n, mu, sigma, alpha)  {
    ## dimension
    K = length(mu)
    
    ## return a matrix
    Mat = matrix(0, n, K)
    
    ## critical value of chisq distribution of dof = K
    ctv.chisq = qchisq(p = 1 - alpha, df = K)
    
    count = 0
    repeat{
        
        draw.try = rmvnorm(n = 1, mean = mu, sigma = sigma)
        draw.try = as.vector(draw.try)
        
        ##within confidence region or not?
        chisq.try = t(draw.try - mu)%*%solve(sigma, draw.try - mu)
        if(chisq.try <= ctv.chisq)
        {
            count = count + 1
            Mat[count, ] = draw.try
            
            if(count == n) break       	
        }## if      	
        
    }## repeat	
    
    return(Mat)
}

## Generate a sequence of tau values in the confidence interval of 
## effects to search over with the permutation test.
get.tau.vector = function( Y, Z, X = NULL, gamma=0.001, grid.size=21, grid.gamma = 100*gamma ) {

    if ( is.null( X ) ) {
        
        te.hat <- mean(Y[Z == 1]) - mean(Y[Z == 0])
        te.se <- sqrt( var(Y[Z == 1])/sum(Z) + var(Y[Z == 0])/sum(1 - Z))
        
    } else {
        lm.tau <- lm( Y ~ Z + X )
        
        te.hat <- as.numeric(coef(lm.tau)["Z"])
        te.se <- as.numeric(sqrt( diag(vcov(lm.tau))["Z"]))
    }

    ## oversample points near the estimated tau-hat
    te.MOE <- qnorm(1 - gamma/2)*te.se
    te.vec <- te.hat + (te.MOE/qnorm(1-grid.gamma)) * qnorm( seq( grid.gamma, 1-grid.gamma, length.out=grid.size ) )
    te.vec
    
    attr( te.vec, "te.hat" ) <- te.hat
    attr( te.vec, "te.se" ) <- te.se
    attr( te.vec, "te.MOE" ) <- te.MOE
    
    te.vec
}

## Calculate a set of different models for effects based on a linear model of W, controlling
## for X (and W).  These models all correspond to estimated nusiance parameters beta.
## Each possible beta gives a collection of different models of effects.
##
## I.e., estimate the coefficients a, d for
## Y ~ a Z + b X + c W + d W:Z
##
## Then for each model, calculate individual imputed treatment effects for all observations and the
## associated science tables of imputed potential outcomes.
##
## @param  Y, Z, X, W: outcome, treatment assignment, covariates, and treatment varying covariates
##                  X and W can be the same for a fully interacted linear model.
## @return
##        te.grid    grid.size x p sized matrix with each row j corresponding to a different model of
##                   effects with specific beta.  p is number of columns in W.
##        Y0.mat, Y1.mat   Two N x grid.size matrices.  Each column j corresponds to a specific model of 
##                   effects with specific beta value.  Each column j corresponds to the same row in te.grid
##         te.mat    N x grid.size matrix with each row i corresponding to treatment effect for 
##                   unit i under model of effects j.
get.testing.grid = function( Y, Z, W, X=NULL, gamma=0.0001, grid.size=150 ) {
    
    ## get sample of treatment effect models to calculate p-values for
    if ( !is.null(X) ) {
        lm.tau <- lm(Y ~ Z + W + Z:W + X)
        ## if have duplicates with X and W, need to drop from consideration
        drp = is.na( coef(lm.tau) )
        cof <- coef(lm.tau)[ !drp ]                
    } else {
        lm.tau <- lm( Y ~ Z + W + Z:W )
        cof <- coef(lm.tau )
    }
    te.model = grep( "Z", names(cof) )
    
    stopifnot( all( rownames( vcov(lm.tau) ) == names(cof) ) )
    
    te.hat <- as.numeric( cof[te.model] )    
    te.cov <- vcov(lm.tau)[te.model,te.model]
    
    te.grid = te.hat
    if ( grid.size > 1 ) {
        ## sample points from our confidence region, focusing on points
        ## close to our point estimate
        te.grid <-  rbind( te.grid, rcrmtvnorm(grid.size - 1, mu = te.hat, sigma = te.cov, alpha = gamma ) )
    }
    
    ## calculate individual treatment effects for different treatment impact models
    ## each row is individual and each column is specific treatment effect model
    te.mat <- t(te.grid %*% t(cbind(1, W))) 
    
    Y0.mat <- Y*(1 - Z) + (Y - te.mat) * Z
    Y1.mat <- Y*Z + (Y + te.mat) * (1 - Z)
    
    list( te.grid=te.grid, te.mat=te.mat, Y0.mat=Y0.mat, Y1.mat=Y1.mat )
}

## Utility function used by FRTCI and FRTCI.interact to actually generate the 
## permutations.  Don't use directly.
generate.permutations = function( Y, Z, test.stat, Y0.mat, Y1.mat, B, get.z.star=NULL, verbose = TRUE, ... ) {
    
    ## SET UP STORAGE MATRICES
    n.te.vec <- ncol(Y0.mat)
    ks.mat <- matrix(NA, nrow = B, ncol = n.te.vec)
    
    ## CALCULATE OBSERVED TEST STATISTICS
    ks.obs <- test.stat( Y, Z, ... )
    
    if ( verbose ) {
        pb <- txtProgressBar(min = 0, max = B, style = 3)
    }
    
    ## RANDOMIZATION DISTRIBUTION
    for(b in 1:B){
        if ( verbose ) {
            setTxtProgressBar(pb, b)
        }
        if ( !is.null( get.z.star ) ) {
            Z.star = get.z.star( Z, ... )
        } else {
            Z.star <- sample(Z)
        }
        
        ## CI method
        ci.out <- sapply(1:n.te.vec, 
                         function(i){ 
                             ## CALCULATE RANDOMIZED VALUE
                             Yobs.star <- Z.star*Y1.mat[,i] + (1 - Z.star)*Y0.mat[,i]
                             
                             ## COMPUTE TEST STATISTICS
                             ks.star <- test.stat(Yobs.star, Z.star, ... )
                             
                             ## OUTPUT
                             ks.star
                         })  
        
        ks.mat[b,] <- as.numeric(ci.out)
    }
    if ( verbose ) {
        close(pb)
    }
    
    ## CALCULATE P-VALUES
    ci.p <- apply(ks.mat, 2, function(ks.star){
        sum(ks.star >= ks.obs)/B
    })
    
    list( ks.obs=ks.obs, ks.mat=ks.mat, ci.p=ci.p )
}
