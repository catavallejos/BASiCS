context("Data examples\n")

test_that("Generated data does not match given seed (spikes case)", {
  Data <- makeExampleBASiCS_Data(WithSpikes = TRUE)
  DataCheck0 <- c(0, 0, 0, 1, 1, 3, 17, 18, 1, 0)
  DataCheck <- as.vector(assay(Data)[1:10,1])
  expect_that(all.equal(DataCheck0, DataCheck), is_true())
})