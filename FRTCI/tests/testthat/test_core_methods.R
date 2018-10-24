library( testthat )
library(FRTCI)

context("Core methods")

test_that("FRTCI Runs", {
  data(ToyData)
  # Set parameters
  B <- 20
  grid.size = 11

  # Test for ideosyncratic treatment effect variation without covariates
  tst = FRTCI(ToyData$Y, ToyData$Z, B=B, test.stat=SKS.stat, grid.size = grid.size, verbose=FALSE)

  expect_false( is.null( tst ) )
  expect_is( tst, "FRTCI.test" )
})



test_that( "Example code from documentation, copied over", {
  data(ToyData)

  B <- 20
  grid.size = 51
  tst = FRTCI(Y, Z, ToyData, B=B, test.stat=SKS.stat, grid.size = grid.size)
  tst

  tst = FRTplug(Y, Z, ToyData, B=B, test.stat=SKS.stat, grid.size = grid.size)

  lmW <- lm( Y ~ x1, ToyData )
  W1 <- model.matrix(lmW)[,-1]
  tst <- FRTCI.interact( Y, Z, W=W1, data=ToyData, B=B )
  tst

  plot.FRTCI.curve( tst )
})


test_that( "Variance ratio test works", {
  data( ToyData )
  variance.ratio.test( Y, Z, data=ToyData )
})
