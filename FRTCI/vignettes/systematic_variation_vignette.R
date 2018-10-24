## ------------------------------------------------------------------------
library( FRTCI )

## ------------------------------------------------------------------------
df = make.randomized.dat( 100, gamma.vec=c(1,1,1,2), beta.vec=c(-1,-1,1) )
str( df )

## ------------------------------------------------------------------------
est.beta( Yobs, Z,  ~ A + B, data=df )

## ------------------------------------------------------------------------
df2 = df[c("A","B","Yobs","Z") ]
names(df2) = c( "AA", "BB", "myY", "myZ" )
rs = est.beta( myY, myZ, myY ~ AA + BB, data=df2 )
rs$beta.hat

## ------------------------------------------------------------------------
vcov( rs )
SE( rs )

## ------------------------------------------------------------------------
confint( rs )

## ------------------------------------------------------------------------
M.ols.ours = est.beta( Yobs, Z, ~ A + B, data=df, method="OLS" )
M.ols.ours
M.ols.ours$beta.hat

## ------------------------------------------------------------------------
M0 = lm( Yobs ~ (A+B) * Z, data=df )
M0

M.ols.ours$beta - coef(M0)[4:6]

## ------------------------------------------------------------------------
est.beta( Yobs, Z, ~ A + B, ~ A + B + C, data=df )

## ------------------------------------------------------------------------
est.beta( Yobs, Z,  ~ A + B, ~ A + B + C, data=df )$beta.hat
est.beta( Yobs, Z,  ~ A + B, ~ A + B + C, data=df, empirical.Sxx = TRUE )$beta.hat

## ------------------------------------------------------------------------
est.beta( Yobs, Z,  ~ A + B, ~ A + B + C, data=df, method="OLS" )$beta.hat

## ------------------------------------------------------------------------
Moracle = calc.beta.oracle( Y.1, Y.0,  Z, ~ A + B, data=df )
Moracle
SE( Moracle )

## ------------------------------------------------------------------------
df = make.randomized.dat( 1000, beta.vec=c(-1,0,0) )

calc.beta.oracle( Y.1, Y.0, Z, ~ A + B, data=df )$beta.hat
est.beta( Yobs, Z, ~ A + B, data=df )$beta.hat
est.beta( Yobs, Z, ~ A + B, data=df, method="OLS" )$beta.hat

## ------------------------------------------------------------------------
df = make.randomized.dat( 1000, beta.vec=c(-1,1,1) )
rs = est.beta( Yobs, Z, ~ A + B, data=df, method="OLS" )
r2 = R2.ITT( rs )
r2    
plot( r2 )
df = make.randomized.dat( 1000, beta.vec=c(-1,1,1), ideo.sd=3 )
rs = est.beta( Yobs, Z, ~ A + B, data=df, method="OLS" )
r2 = R2.ITT( rs )
r2    
plot( r2, ADD=TRUE, col="green" )

## ------------------------------------------------------------------------
df = make.randomized.dat( 1000, beta.vec=c(-1,1,0) )
rs = est.beta( Yobs, Z, ~ A + B, data=df, method="OLS" )
r2 = R2.ITT( rs )
r2    
plot( r2 )

## ------------------------------------------------------------------------
plot( df$tau ~ df$A )

## ------------------------------------------------------------------------
df = make.randomized.dat( 1000, beta.vec=c(-1,1,1), ideo.sd=1 )

rs = est.beta( Yobs, Z, ~ A + B, data=df )
r2 = R2.ITT( rs )
r2    
plot( r2, col="green" )

# adjusted
rs = est.beta( Yobs, Z, ~ A + B, ~ A + B, data=df )
r2 = R2.ITT( rs )
r2    
plot( r2, ADD=TRUE )

# adjusted + OLS
rs = est.beta( Yobs, Z, ~ A + B, ~ A + B + C, data=df, empirical.Sxx = TRUE )
r2 = R2.ITT( rs )
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
rs = est.beta.LATE( Z, D, Yobs, ~ A + B, data=data )
rs
rs$beta.hat
SE( rs )

# error
rs$beta.hat - beta

r2 = R2.LATE( rs )
r2
plot( r2 )

## ------------------------------------------------------------------------
rs = est.beta.LATE( Z, D, Yobs,  ~ A + B, data=data, method="2SLS" )
rs
SE( rs )
rs$beta.hat
rs$beta.hat - beta

