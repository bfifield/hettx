library(FRTCI)

context("Systematic estimation methods")

test_that("estimate oracle", {

  df = make.randomized.dat( 100, beta.vec=c(-1,-1,1) )

  orc = calc.beta.oracle( Y.1, Y.0, Z,  ~ A + B, data=df )
  expect_equivalent( coef( orc ), c(-1, -1, 1) )
  coef( orc )

  ess = est.beta( Yobs, Z,  ~ A + B, data=df )
  coef( ess )

  df2 = df[c("A","B","Yobs","Z") ]
  names(df2) = c( "AA", "BB", "myY", "myZ" )
  rs = est.beta( myY, myZ, myY ~ AA + BB, data=df2 )
  coef( rs )

  expect_equivalent( coef( ess ), coef( rs ) )

} )



test_that( "OLS estimation corresponds to lm", {

  df = make.randomized.dat( 100, beta.vec=c(-1,-1,1) )


  # The simple interaction approach using OLS
  M0 = lm( Yobs ~ (A+B) * Z, data=df )
  M0

  # same in our wrapper.  Sanity check
  M.ols = est.beta( Yobs, Z, Yobs ~ A + B, data=df, method="OLS" )
  M.ols

  expect_equivalent( coef( M.ols ), coef(M0)[4:6] )

} )


