context("RBF settings")

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
  pp <- BASiCS_PriorParam(d)
  c <- run_MCMC(d,
    N = 10,
    Thin = 2,
    Burn = 4,
    Regression = TRUE,
    PriorParam = pp,
    MinGenesPerRBF = 5
  )
  l1 <- c@parameters$RBFLocations
  expect_equal(round(l1[, 1], digits = 2), c(0.93, 2.3, 3.67, 5.03))

  d <- makeExampleBASiCS_Data()
  pp <- BASiCS_PriorParam(d)
  c <- run_MCMC(d,
    N = 10,
    Thin = 2,
    Burn = 4,
    Regression = TRUE,
    PriorParam = pp,
    MinGenesPerRBF = 2
  )
  l1 <- c@parameters$RBFLocations
  expect_equal(round(l1[, 1], digits = 2), c(1.32, 2.29, 3.26, 4.23, 5.2))
})

test_that("MinGenesPerRBF fails with high value (low k)", {
  d <- makeExampleBASiCS_Data()
  pp <- BASiCS_PriorParam(d)
  expect_error(
    run_MCMC(d,
      N = 10,
      Thin = 2,
      Burn = 4,
      Regression = TRUE,
      PriorParam = pp,
      MinGenesPerRBF = 20
    ),
    "The number of basis functions needs to be >= 4"
  )
})
