context("DenoisedCounts/Rates")

test_that("DenoisedCounts and Rates require matching dimensions", {
  data <- makeExampleBASiCS_Data()
  chain <- run_MCMC(
    data,
    N = 10,
    Thin = 2,
    Burn = 4,
    Regression = FALSE,
    WithSpikes = TRUE
  )
  expect_error(
    BASiCS_DenoisedCounts(data, ChainSC),
    "Chain and Data are different dimensions"
  )
  expect_error(
    BASiCS_DenoisedRates(data, ChainSC),
    "Chain and Data are different dimensions"
  )
})
