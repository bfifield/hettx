context("Core methods")

test_that("FRTCI Runs", {
  data(ToyData)
  # Set parameters
  B <- 20
  grid.size = 11

  # Test for ideosyncratic treatment effect variation without covariates
  tst = fishpidetect(ToyData$Y, ToyData$Z, B=B, grid.size = grid.size, verbose=FALSE)

  expect_false( is.null( tst ) )
  expect_is( tst, "FRTCI.test" )
})



test_that( "Example code from documentation, copied over", {
  data(ToyData)

  B <- 20
  grid.size = 51
  tst = fishpidetect(ToyData$Y, ToyData$Z, B=B, grid.size = grid.size)
  tst
  plot( tst )
  
  tst = fishpidetect(ToyData$Y, ToyData$Z, plugin = TRUE)

  lmW <- lm( Y ~ x1, ToyData )
  W1 <- model.matrix(lmW)[,-1]
  tst <- fishpidetect(ToyData$Y, ToyData$Z, W=as.matrix(W1), B=B )
  tst

})


test_that( "Variance ratio test works", {
  data( ToyData )
  variance.ratio.test( Y, Z, data=ToyData )
})
