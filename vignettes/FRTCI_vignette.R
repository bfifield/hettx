## ---- echo=FALSE, message=FALSE, warning=FALSE---------------------------
par( mgp=c(1.8,0.8,0), mar=c(2.5, 2.5,2,0.1) )
knitr::opts_chunk$set(fig.width = 8)

## ---- message=FALSE, warning=FALSE---------------------------------------
library( mvtnorm )
library( tidyverse )
library( hettx )
data( ToyData )

## ---- echo=TRUE----------------------------------------------------------
data( ToyData )
head( ToyData )
td = gather( ToyData, x1, x2, x3, x4, key="X", value="value" )
td = gather( td, Y, tau, key="outcome", value="value2" )
ggplot( td, aes( x=value, y=value2, col=as.factor(Z) ) ) +
        facet_grid( outcome ~ X, scales="free" ) +
        geom_point( alpha=0.5, size=0.5) + 
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
tst1 = detect.idiosyncratic( Y ~ Z, data=ToyData, B=B, grid.size = grid.size, verbose=FALSE )
print( tst1 )

## ---- echo=TRUE, cache=TRUE----------------------------------------------
tst2 = detect.idiosyncratic( Y ~ Z, data=ToyData, control.formula = ~ x1 + x2 + x3 + x4, B=B, test.stat=SKS.stat.cov, 
              verbose=FALSE )
print( tst2 )

## ---- echo=TRUE, cache=TRUE----------------------------------------------
tst2b = detect.idiosyncratic( Y ~ Z, data=ToyData, control.formula = ~ x3 + x4, B=B, test.stat=SKS.stat.cov, 
              verbose=FALSE )
print( tst2b )

## ---- echo=TRUE, cache=TRUE----------------------------------------------
N = nrow(ToyData)
ToyData$Xfake = rnorm( N )
tst1b = detect.idiosyncratic( Y ~ Z, data=ToyData, control.formula = ~ Xfake, B=B, test.stat=SKS.stat.cov, 
              verbose=FALSE )
print( tst1b )

## ----ideo_beyond_systematic, echo=TRUE, cache=TRUE-----------------------
B = 20 

tst3a1 = detect.idiosyncratic( Y ~ Z, data=ToyData, interaction.formula = ~ x1, B=B, test.stat=SKS.stat.int.cov, 
              verbose=FALSE )
print( tst3a1 )

## ---- echo=TRUE, cache=TRUE----------------------------------------------
tst3a2 <- detect.idiosyncratic( Y ~ Z, data=ToyData, 
                        interaction.formula = ~ x1, 
                        control.formula = ~ x3 + x4,
                        B=B, test.stat=SKS.stat.int.cov, 
              verbose=FALSE )
print( tst3a2 )

## ---- echo=TRUE, cache=TRUE----------------------------------------------
tst3a2b <- detect.idiosyncratic( Y ~ Z, data=ToyData, control.formula = ~ x2 + x3 + x4, 
                         interaction.formula = ~ x1, B=B, test.stat=SKS.stat.int.cov, 
              verbose=FALSE )
print( tst3a2b )

## ----beyond_x1_x2, echo=TRUE, cache=TRUE---------------------------------
tst3b <- detect.idiosyncratic( Y ~ Z, data=ToyData,
                         interaction.formula = ~ x1 + x2, B=B, test.stat=SKS.stat.int.cov, 
              verbose=FALSE )
print( tst3b )

## ---- echo=TRUE----------------------------------------------------------
tst3c <- detect.idiosyncratic( Y ~ Z, data=ToyData,
                         interaction.formula = ~ x1 + x2, 
                       control.formula = ~ x3 + x4,
                       B=B, test.stat=SKS.stat.int.cov, 
              verbose=FALSE )
print( tst3c )

## ----full_correction, echo=TRUE, cache=TRUE------------------------------
tst3d <- detect.idiosyncratic( Y ~ Z, data=ToyData,
                         interaction.formula = ~ x1 + x2 + x3 + x4, 
                       B=B, test.stat=SKS.stat.int.cov, 
              verbose=FALSE )
print( tst3d )

## ----display, echo=TRUE--------------------------------------------------
tests = list( no_cov=tst1, useless_cov=tst1b, all_covariates=tst2, 
              non_tx_covariates_only=tst2b, het_beyond_x1 = tst3a1,
              het_beyond_x1_with_x3_x4_cov=tst3a2, het_beyond_x1_with_all_cov=tst3a2,
              het_beyond_x1_x2=tst3b, 
              het_beyond_x1_x2_with_cov=tst3c, het_beyond_all=tst3d )

agg.res = purrr::map( tests, get.p.value  ) %>%
  purrr::map( as.list )
agg.res = bind_rows( agg.res, .id = "test" )
agg.res

## ----cautionary_tale, echo=TRUE------------------------------------------
ll1 = lm( Y ~ Z + x1 + x2 + x3 + x4, data=ToyData )
print( summary( ll1 ) )

## ----cautionary_tale_2, echo=TRUE----------------------------------------
plot( ecdf( resid(ll1)[ToyData$Z==1] ), pch=".", main="Residual CDFs of treatment and control" )
plot( ecdf( resid(ll1)[ToyData$Z==0] ), pch=".", col="red", add=TRUE )

