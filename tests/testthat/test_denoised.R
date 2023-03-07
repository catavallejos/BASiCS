test_that("DenoisedCounts and Rates require matching dimensions", {
  data <- makeExampleBASiCS_Data()
  data(ChainSC)
  expect_error(
    BASiCS_DenoisedCounts(data, ChainSC),
    "Chain and Data are different dimensions"
  )
  expect_error(
    BASiCS_DenoisedRates(data, ChainSC),
    "Chain and Data are different dimensions"
  )
})
