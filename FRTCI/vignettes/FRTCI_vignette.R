## ---- echo=FALSE, message=FALSE, warning=FALSE---------------------------
library( plyr )
library( mvtnorm )
library( tidyverse )
library( FRTCI )
par( mgp=c(1.8,0.8,0), mar=c(2.5, 2.5,2,0.1) )
knitr::opts_chunk$set(fig.width = 8)

## ---- message=FALSE, warning=FALSE---------------------------------------
library( plyr )
library( mvtnorm )
library( tidyverse )
library( FRTCI )

## ---- echo=TRUE----------------------------------------------------------
head( ToyData )
td = gather( ToyData, x1, x2, x3, x4, key="X", value="value" )
td = gather( td, Y, tau, key="outcome", value="value2" )
ggplot( td, aes( x=value, y=value2, col=as.factor(Z) ) ) +
        facet_grid( outcome ~ X, scales="free" ) +
        geom_point( alpha=0.5) + 
        geom_smooth( method="loess", se=FALSE ) +
        labs( x="Covariates", y="" )

## ---- echo=TRUE----------------------------------------------------------
mean( ToyData$tau )

## ---- echo=TRUE----------------------------------------------------------
par( mfrow=c(1,2) )
ll0 = lm( Y ~ Z, data=ToyData )
plot( ecdf( resid(ll0)[ToyData$Z==1] ), pch=".", main="Marginal CDFs of \n treatment and control")
plot( ecdf( resid(ll0)[ToyData$Z==0] ), pch=".", col="red", add=TRUE )

ll1 = lm( Y ~ Z + x1 + x2 + x3 + x4, data=ToyData )
plot( ecdf( resid(ll1)[ToyData$Z==1] ), pch=".", main="Residual CDFs of \n treatment and control" )
plot( ecdf( resid(ll1)[ToyData$Z==0] ), pch=".", col="red", add=TRUE )

## ---- echo=TRUE----------------------------------------------------------
M0 <- lm( Y ~ Z * (x1+x2+x3), data=ToyData )
round( coef( M0 ), digits=1 )

## ---- echo=TRUE, results='asis'------------------------------------------
B <- 100
grid.size = 21

## ---- echo=TRUE, cache=TRUE----------------------------------------------
tst1 = FRTCI( Y=ToyData$Y, Z=ToyData$Z, B=B, grid.size = grid.size, verbose=FALSE )
print( tst1 )

## ---- echo=TRUE, cache=TRUE----------------------------------------------
tst1b = FRTCI( Y, Z, data=ToyData, B=B, grid.size = grid.size, verbose=FALSE )

## ---- echo=TRUE----------------------------------------------------------
Xfull <- model.matrix(Y ~ 0 + x1 + x2 + x3 + x4, data=ToyData )

# get points to check in the confidence interval
te.vec = get.tau.vector( Y=ToyData$Y, Z=ToyData$Z, X=Xfull, grid.size=grid.size )

# do permutation test for all these points and maximize
tst2 = FRTCI( Y=ToyData$Y, Z=ToyData$Z, B=B, te.vec = te.vec, test.stat=SKS.stat.cov, 
              X=Xfull, verbose=FALSE )
print( tst2 )

## ---- echo=TRUE----------------------------------------------------------
cat( "Adjusting just using the non-treatment related covariates\n" )
X34 <- model.matrix( Y ~ 0 + x3 + x4, data=ToyData )

te.vec = get.tau.vector( Y=ToyData$Y, Z=ToyData$Z, X=X34, grid.size = grid.size )

tst2b = FRTCI( Y=ToyData$Y, Z=ToyData$Z, B=B, te.vec = te.vec, test.stat=SKS.stat.cov, 
               X=X34, verbose=FALSE )
print( tst2b )

## ---- echo=TRUE----------------------------------------------------------
N = length(ToyData$Y)
Xfake = matrix( rnorm( N ), nrow=N )

tst1b = FRTCI( Y=ToyData$Y, Z=ToyData$Z, B=B, test.stat=SKS.stat.cov, X=Xfake,
               grid.size = grid.size, verbose=FALSE )
print( tst1b )

## ----ideo_beyond_systematic, echo=TRUE-----------------------------------
B = 20 

W1 <- model.matrix(Y ~ 0 + x1, data=ToyData)

tst3a1 <- FRTCI.interact( Y=ToyData$Y, Z=ToyData$Z, W=W1, B=B, verbose=FALSE )
print( tst3a1 )

## ---- echo=TRUE----------------------------------------------------------
tst3a2 <- FRTCI.interact( Y=ToyData$Y, Z=ToyData$Z, W=W1, X=X34, B=B, verbose=FALSE )
print( tst3a2 )

## ---- echo=TRUE----------------------------------------------------------
tst3a2b <- FRTCI.interact( Y=ToyData$Y, Z=ToyData$Z, W=W1, X=Xfull, B=B, verbose=FALSE )
print( tst3a2b )

## ----beyond_x1_x2, echo=TRUE---------------------------------------------
W12 <- model.matrix( Y ~ 0 + x1 + x2, data=ToyData )

tst3b <- FRTCI.interact( Y=ToyData$Y, Z=ToyData$Z, W=W12, B=B, verbose=FALSE )
print( tst3b)

## ---- echo=TRUE----------------------------------------------------------
tst3c <- FRTCI.interact( Y=ToyData$Y, Z=ToyData$Z, W=W12, X=X34, B=B, verbose=FALSE )
print( tst3c )

## ----full_correction, echo=TRUE------------------------------------------
Wfull <- model.matrix( Y ~ 0 + x1 + x2 + x3 + x4, data=ToyData)

tst3d <- FRTCI.interact( Y=ToyData$Y, Z=ToyData$Z, W=Wfull, B=B, verbose=FALSE )
print( tst3d )

## ----display, echo=TRUE--------------------------------------------------
tests = list( no_cov=tst1, useless_cov=tst1b, all_covariates=tst2, non_tx_covariates_only=tst2b, het_beyond_x1 = tst3a1, het_beyond_x1_with_x3_x4_cov=tst3a2, het_beyond_x1_with_all_cov=tst3a2, het_beyond_x1_x2=tst3b, het_beyond_x1_x2_with_cov=tst3c, het_beyond_all=tst3d )

agg.res = map_df( tests, get.p.value, .id = "test" )
agg.res

## ----cautionary_tale, echo=TRUE------------------------------------------
ll1 = lm( Y ~ Z + x1 + x2 + x3 + x4, data=ToyData )
print( summary( ll1 ) )

## ----cautionary_tale_2, echo=TRUE----------------------------------------
plot( ecdf( resid(ll1)[ToyData$Z==1] ), pch=".", main="Residual CDFs of treatment and control" )
plot( ecdf( resid(ll1)[ToyData$Z==0] ), pch=".", col="red", add=TRUE )

