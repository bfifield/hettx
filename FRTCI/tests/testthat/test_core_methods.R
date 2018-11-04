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

test_that( "Every test statistic works", {

    data(ToyData)

    B <- 20
    grid.size = 51

    ## -------------
    ## Test defaults
    ## -------------
    X <- model.matrix(~ x1 + x2, data = ToyData)[,-1]
    W <- model.matrix(~ x3 + x4, data = ToyData)[,-1]
    W.fact <- as.factor(sample(c("A", "B"), nrow(W), replace = TRUE))
    
    tst = fishpidetect(ToyData$Y, ToyData$Z, B=B, grid.size=grid.size)
    tst = fishpidetect(ToyData$Y, ToyData$Z, X=X, B=B, grid.size=grid.size)
    tst = fishpidetect(ToyData$Y, ToyData$Z, W=W, B=B, grid.size=grid.size)
    tst = fishpidetect(ToyData$Y, ToyData$Z, W=W, X=X, B=B, grid.size=grid.size)

    ## ------------------------------------------------
    ## Test statistics for no adjustment or interaction
    ## ------------------------------------------------
    tst = fishpidetect(ToyData$Y, ToyData$Z, test.stat=KS.stat, B=B, grid.size=grid.size)
    tst = fishpidetect(ToyData$Y, ToyData$Z, test.stat=SKS.stat, B=B, grid.size=grid.size)
    tst = fishpidetect(ToyData$Y, ToyData$Z, test.stat=rq.stat, B=B, grid.size=grid.size)

    ## ----------------------------------------------
    ## Test statistics for adjustment, no interaction
    ## ----------------------------------------------
    tst = fishpidetect(ToyData$Y, ToyData$Z, X=X, test.stat=SKS.stat.cov.pool, B=B, grid.size=grid.size)
    tst = fishpidetect(ToyData$Y, ToyData$Z, X=X, test.stat=SKS.stat.cov, B=B, grid.size=grid.size)
    tst = fishpidetect(ToyData$Y, ToyData$Z, X=X, test.stat=SKS.stat.cov.rq, B=B, grid.size=grid.size)
    tst = fishpidetect(ToyData$Y, ToyData$Z, X=X, test.stat=rq.stat.cond.cov, B=B, grid.size=grid.size)
    tst = fishpidetect(ToyData$Y, ToyData$Z, X=X, test.stat=rq.stat.uncond.cov, B=B, grid.size=grid.size)

    ## ----------------------------------------------
    ## Test statistics for interaction, no adjustment
    ## ----------------------------------------------
    tst = fishpidetect(ToyData$Y, ToyData$Z, W=W, test.stat=SKS.stat.int.cov.pool, B=B, grid.size=grid.size)
    tst = fishpidetect(ToyData$Y, ToyData$Z, W=W, test.stat=SKS.stat.int.cov, B=B, grid.size=grid.size)
    tst = fishpidetect(ToyData$Y, ToyData$Z, W=W.fact, test.stat=WSKS.t, B=B, grid.size=grid.size)
    tst = fishpidetect(ToyData$Y, ToyData$Z, W=W.fact, test.stat=SKS.pool.t, B=B, grid.size=grid.size)

    ## ----------------------------------------------
    ## Test statistics for interaction and adjustment
    ## ----------------------------------------------
    tst = fishpidetect(ToyData$Y, ToyData$Z, W=W, X=X, test.stat=SKS.stat.int.cov.pool, B=B, grid.size=grid.size)
    tst = fishpidetect(ToyData$Y, ToyData$Z, W=W, X=X, test.stat=SKS.stat.int.cov, B=B, grid.size=grid.size)
    
})
