context("ESS")

test_that("ess matches effectiveSize", {
  data("ChainSC")
  expect_equal(
    ess(ChainSC@parameters[["mu"]]),
    effectiveSize(ChainSC@parameters[["mu"]])
  )
})

test_that("ess doesn't fail for zero variance", {
  x <- matrix(rep(1, 10000), ncol = 2)
  expect_error(e <- ess(x), NA)
  expect_equal(e, rep(0, 2))
})

test_that("ess doesn't fail for NA", {
  x <- matrix(rnorm(100), ncol = 2)
  x[1:5, ] <- NA
  expect_error(e <- ess(x), NA)
  expect_equal(e, rep(0, 2))
})
