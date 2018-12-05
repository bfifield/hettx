## ------------------------------------------------------------------------
library( hettx )

## ------------------------------------------------------------------------
df = make.randomized.dat( 10000, gamma.vec=c(1,1,1,2), beta.vec=c(-1,-1,1,0) )
str( df )

## ------------------------------------------------------------------------
rs = estimate.systematic( Yobs ~ Z,  interaction.formula = ~ A + B, data=df )
rs

## ------------------------------------------------------------------------
vcov( rs )
SE( rs )

## ------------------------------------------------------------------------
confint( rs )

## ------------------------------------------------------------------------
M.ols.ours = estimate.systematic( Yobs ~ Z, ~ A + B, data=df, method="OLS" )
M.ols.ours
M.ols.ours$beta.hat

## ------------------------------------------------------------------------
M0 = lm( Yobs ~ (A+B) * Z, data=df )
M0

## ------------------------------------------------------------------------
M.ols.ours$beta - coef(M0)[4:6]

## ------------------------------------------------------------------------
estimate.systematic( Yobs ~ Z, interaction.formula = ~ A + B, 
          control.formula = ~ C, data=df )

## ------------------------------------------------------------------------
rsA2 = estimate.systematic( Yobs ~ Z,  ~ A + B, ~ A + B + C, data=df )
coef( rsA2 )

## ------------------------------------------------------------------------
rsB = estimate.systematic( Yobs ~ Z,  ~ A + B, ~ C, data=df, method = "OLS" )
coef( rsB )
rsB2 = estimate.systematic( Yobs ~ Z,  ~ A + B, ~ A + B + C, data=df, method = "OLS" )
coef( rsB2 )

## ------------------------------------------------------------------------
rsB.lm = lm( Yobs ~ Z * (A+B) + C, data=df )
coef( rsB.lm )

# Doesn't work?
#library( texreg )
#screenreg( rsB, rsB2)

cbind( C.only=coef( rsB ), ABC=coef( rsB2 ), lmC=coef( rsB.lm )[c(2,6,7)])

## ------------------------------------------------------------------------
Moracle = estimate.systematic( Y.1 + Y.0 ~ Z, ~ A + B, data=df )
Moracle
SE( Moracle )

## ------------------------------------------------------------------------
df = make.randomized.dat( 1000, beta.vec=c(-1,0,0) )

estimate.systematic( Y.1 + Y.0 ~ Z, ~ A + B, data=df )$beta.hat
estimate.systematic( Yobs ~ Z, ~ A + B, data=df )$beta.hat
estimate.systematic( Yobs ~ Z, ~ A + B, data=df, method="OLS" )$beta.hat

## ------------------------------------------------------------------------
df = make.randomized.dat( 1000, beta.vec=c(-1,1,1) )
rs = estimate.systematic( Yobs ~ Z, ~ A + B, data=df, method="OLS" )
r2 = R2( rs )
r2    

## ------------------------------------------------------------------------
df = make.randomized.dat( 1000, beta.vec=c(-1,1,1), ideo.sd=3 )
rs = estimate.systematic( Yobs ~ Z, ~ A + B, data=df, method="OLS" )
r2 = R2( rs )
r2    

## ------------------------------------------------------------------------
plot( r2 )
plot( r2, ADD=TRUE, col="green" )

## ------------------------------------------------------------------------
df = make.randomized.dat( 1000, beta.vec=c(-1,1,0) )
rs = estimate.systematic( Yobs ~ Z, ~ A + B, data=df, method="OLS" )
r2 = R2( rs )
r2    
plot( r2 )

## ------------------------------------------------------------------------
plot( df$tau ~ df$A )

## ------------------------------------------------------------------------
df = make.randomized.dat( 1000, beta.vec=c(-1,1,1), ideo.sd=1 )

rs = estimate.systematic( Yobs ~ Z, ~ A + B, data=df )
r2 = R2( rs )
r2    
plot( r2, col="green" )

# adjusted
rs = estimate.systematic( Yobs ~ Z, ~ A + B, ~ C, data=df )
r2 = R2( rs )
r2    
plot( r2, ADD=TRUE )

# adjusted + OLS
rs = estimate.systematic( Yobs ~ Z, ~ A + B, ~ C, data=df, method = "OLS" )
r2 = R2( rs )
r2    
plot( r2, ADD=TRUE, col="blue" )

## ------------------------------------------------------------------------
beta = c(-1,6,0)
n = 1000

data = make.randomized.compliance.dat( n, beta.vec=beta )
names(data)
head( data )
zd = with( data, interaction( Z, D, sep="-" ) )

## ----observed_subgroups--------------------------------------------------
boxplot( Yobs ~ zd, data=data )

## ------------------------------------------------------------------------
par( mfrow=c(1,2), mgp=c(1.8,0.8,0), mar=c(3,3,0.5,0.5) )
plot( Y.1 - Y.0 ~ A, data=data, col=as.factor(data$S), pch=19, cex=0.5 )
plot( Y.1 - Y.0 ~ B, data=data, col=as.factor(data$S), pch=19, cex=0.5 )
legend( "topleft", legend=levels( as.factor( data$S ) ), pch=19, col=1:3 )

## ------------------------------------------------------------------------
rs = estimate.systematic( Yobs ~ D | Z, ~ A + B, data=data )
rs
rs$beta.hat
SE( rs )

# error
rs$beta.hat - beta

r2 = R2( rs )
r2
plot( r2 )

## ------------------------------------------------------------------------
rs = estimate.systematic( Yobs ~ Z | D,  ~ A + B, data=data, method="2SLS" )
rs
SE( rs )
rs$beta.hat
rs$beta.hat - beta

