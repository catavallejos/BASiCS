test_that("RBFMinMax works as expected", {
  set.seed(42)
  d <- makeExampleBASiCS_Data()
  pp <- BASiCS_PriorParam(d, RBFMinMax = FALSE, k = 10)
  c <- run_MCMC(d, N = 10, Thin = 2, Burn = 4, Regression = TRUE, PriorParam = pp)
  l1 <- c@parameters$RBFLocations

  pp <- BASiCS_PriorParam(d, RBFMinMax = TRUE, k = 12)
  c <- run_MCMC(d, N = 10, Thin = 2, Burn = 4, Regression = TRUE, PriorParam = pp)
  l2 <- c@parameters$RBFLocations
  expect_equal(l2[2:9], l1[1:8])
})


test_that("MinGenesPerRBF works as expected", {
  set.seed(42)
  d <- makeExampleBASiCS_Data()
  pp <- BASiCS_PriorParam(d, MinGenesPerRBF = 5)
  c <- run_MCMC(d,
    N = 10,
    Thin = 2,
    Burn = 4,
    Regression = TRUE,
    PriorParam = pp
  )
  l1 <- c@parameters$RBFLocations
  expect_equal(round(l1[, 1], digits = 2), c(1.84, 2.3, 2.76, 3.21))

  d <- makeExampleBASiCS_Data()
  pp <- BASiCS_PriorParam(d, MinGenesPerRBF = 2)
  c <- run_MCMC(d,
    N = 10,
    Thin = 2,
    Burn = 4,
    Regression = TRUE,
    PriorParam = pp
  )
  l1 <- c@parameters$RBFLocations
  expect_equal(round(l1[, 1], digits = 2), c(1.75, 2.18, 2.61, 3.04, 4.34))
})

test_that("MinGenesPerRBF fails with high value (low k)", {
  d <- makeExampleBASiCS_Data()
  pp <- BASiCS_PriorParam(d, MinGenesPerRBF = 20)
  expect_error(
    run_MCMC(d,
      N = 10,
      Thin = 2,
      Burn = 4,
      Regression = TRUE,
      PriorParam = pp
    ),
    "The number of basis functions needs to be >= 4"
  )
})
