## ---- echo=FALSE, message=FALSE, warning=FALSE---------------------------
par( mgp=c(1.8,0.8,0), mar=c(2.5, 2.5,2,0.1) )
knitr::opts_chunk$set(fig.width = 8)

## ---- message=FALSE, warning=FALSE---------------------------------------
library( mvtnorm )
library( ggplot2 )
library( hettx )
library( dplyr )
library( tidyr )
library( purrr )
data( ToyData )

## ---- echo=TRUE----------------------------------------------------------
data( ToyData )
head( ToyData )
td <- pivot_longer( ToyData, cols=c(x1, x2, x3, x4), names_to="X", values_to="value" )
td <- pivot_longer( td, cols=c(Y, tau), names_to="outcome", values_to="value2" )
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
B <- 20
grid.size = 11

## ---- echo=TRUE----------------------------------------------------------
tst3c <- detect_idiosyncratic( Y ~ Z, data=ToyData,
                         interaction.formula = ~ x1 + x2, 
                       control.formula = ~ x3 + x4,
                       B=B, test.stat="SKS.stat.int.cov", 
                       verbose=FALSE )
summary( tst3c )

## ------------------------------------------------------------------------
get.p.value( tst1b )

## ----display, echo=TRUE--------------------------------------------------
tests = list( no_cov=tst1, useless_cov=tst1b, all_covariates=tst2, 
              non_tx_covariates_only=tst2b, het_beyond_x1 = tst3a1,
              het_beyond_x1_with_x3_x4_cov=tst3a2, het_beyond_x1_with_all_cov=tst3a2,
              het_beyond_x1_x2=tst3b, 
              het_beyond_x1_x2_with_cov=tst3c, het_beyond_all=tst3d )

agg.res <- purrr::map( tests, get.p.value ) %>%
  purrr::map( as.list )
agg.res <- bind_rows( agg.res, .id = "test" )
agg.res

## ----cautionary_tale, echo=TRUE------------------------------------------
ll1 = lm( Y ~ Z + x1 + x2 + x3 + x4, data=ToyData )
print( summary( ll1 ) )

## ----cautionary_tale_2, echo=TRUE----------------------------------------
plot( ecdf( resid(ll1)[ToyData$Z==1] ), pch=".", main="Residual CDFs of treatment and control" )
plot( ecdf( resid(ll1)[ToyData$Z==0] ), pch=".", col="red", add=TRUE )

## ------------------------------------------------------------------------
variance.ratio.test( ToyData$Y, ToyData$Z )

## ------------------------------------------------------------------------
test.stat.info()

