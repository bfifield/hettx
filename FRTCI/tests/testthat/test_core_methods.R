context("Core methods")

test_that("FRTCI Runs", {
    data(ToyData)
    ## Set parameters
    B <- 20
    grid.size = 11

    ## Test for ideosyncratic treatment effect variation without covariates
    tst = fishpidetect(Y ~ Z, data = ToyData, B=B, grid.size = grid.size, verbose=FALSE)

    expect_false( is.null( tst ) )
    expect_is( tst, "FRTCI.test" )
})



test_that( "Example code from documentation, copied over", {
    data(ToyData)

    B <- 20
    grid.size = 51
    tst = fishpidetect(Y ~ Z, data = ToyData, B=B, grid.size = grid.size)
    tst
    plot( tst )
    
    tst = fishpidetect(Y ~ Z, data = ToyData, plugin = TRUE)

    tst <- fishpidetect(Y ~ Z, data = ToyData, interaction.formula = ~ x1, B=B )
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
    ToyData$W.fact <- as.factor(sample(c("A", "B"), nrow(ToyData), replace = TRUE))
    
    tst = fishpidetect(Y ~ Z, data = ToyData, B=B, grid.size=grid.size)
    tst = fishpidetect(Y ~ Z, data = ToyData, control.formula = ~ x1 + x2, B=B, grid.size=grid.size)
    tst = fishpidetect(Y ~ Z, data = ToyData, interaction.formula = ~ x3 + x4, B=B, grid.size=grid.size)
    tst = fishpidetect(Y ~ Z, data = ToyData, interaction.formula = ~ x3 + x4, control.formula = ~ x1 + x2, B=B, grid.size=grid.size)

    ## ------------------------------------------------
    ## Test statistics for no adjustment or interaction
    ## ------------------------------------------------
    tst = fishpidetect(Y ~ Z, data = ToyData, test.stat=KS.stat, B=B, grid.size=grid.size)
    tst = fishpidetect(Y ~ Z, data = ToyData, test.stat=SKS.stat, B=B, grid.size=grid.size)
    tst = fishpidetect(Y ~ Z, data = ToyData, test.stat=rq.stat, B=B, grid.size=grid.size)

    ## ----------------------------------------------
    ## Test statistics for adjustment, no interaction
    ## ----------------------------------------------
    tst = fishpidetect(Y ~ Z, data = ToyData, control.formula = ~ x1 + x2, test.stat=SKS.stat.cov.pool, B=B, grid.size=grid.size)
    tst = fishpidetect(Y ~ Z, data = ToyData, control.formula = ~ x1 + x2, test.stat=SKS.stat.cov, B=B, grid.size=grid.size)
    tst = fishpidetect(Y ~ Z, data = ToyData, control.formula = ~ x1 + x2, test.stat=SKS.stat.cov.rq, B=B, grid.size=grid.size)
    tst = fishpidetect(Y ~ Z, data = ToyData, control.formula = ~ x1 + x2, test.stat=rq.stat.cond.cov, B=B, grid.size=grid.size)
    tst = fishpidetect(Y ~ Z, data = ToyData, control.formula = ~ x1 + x2, test.stat=rq.stat.uncond.cov, B=B, grid.size=grid.size)

    ## ----------------------------------------------
    ## Test statistics for interaction, no adjustment
    ## ----------------------------------------------
    tst = fishpidetect(Y ~ Z, data = ToyData, interaction.formula = ~ x3 + x4, test.stat=SKS.stat.int.cov.pool, B=B, grid.size=grid.size)
    tst = fishpidetect(Y ~ Z, data = ToyData, interaction.formula = ~ x3 + x4, test.stat=SKS.stat.int.cov, B=B, grid.size=grid.size)
    tst = fishpidetect(Y ~ Z, data = ToyData, interaction.formula = ~ W.fact, test.stat=WSKS.t, B=B, grid.size=grid.size)
    tst = fishpidetect(Y ~ Z, data = ToyData, interaction.formula = ~ W.fact, test.stat=SKS.pool.t, B=B, grid.size=grid.size)

    ## ----------------------------------------------
    ## Test statistics for interaction and adjustment
    ## ----------------------------------------------
    tst = fishpidetect(Y ~ Z, data = ToyData, interaction.formula = ~ x3 + x4, control.formula = ~ x1 + x2, test.stat=SKS.stat.int.cov.pool, B=B, grid.size=grid.size)
    tst = fishpidetect(Y ~ Z, data = ToyData, interaction.formula = ~ x3 + x4, control.formula = ~ x1 + x2, test.stat=SKS.stat.int.cov, B=B, grid.size=grid.size)
    
})
