context("Data examples\n")

test_that("Generated data does not match given seed (spikes case)", {
  set.seed(1)
  Data <- makeExampleBASiCS_Data(WithSpikes = TRUE)
  
  # Checks some specific cells
  DataCheck0 <- c(0,  5,  1,  0,  0,  1, 17, 11,  1,  2)
  DataCheck <- as.vector(assay(Data)[1:10,1])
  expect_equal(DataCheck0, DataCheck)
  
  # Checks isSpike info
  TechCheck0 <- c(rep(FALSE, 50), rep(TRUE, 20))
  TechCheck <- c(rep(FALSE, nrow(Data)), rep(TRUE, nrow(altExp(Data))))
  expect_equal(TechCheck0, TechCheck)
  
  # Checks total count for spike in genes
  TechCount0 <- c(11087, 101, 344, 633, 20)
  TechCount <- matrixStats::rowSums2(assay(altExp(Data)))[1:5]
  expect_equal(TechCount0, TechCount)
})

test_that("Generated data does not match given seed (no spikes case)", {
  
  # Checks some specific cells
  set.seed(2)
  Data <- makeExampleBASiCS_Data(WithSpikes = FALSE)
  DataCheck0 <- c(0, 5, 3, 6, 28, 4, 3, 170, 0, 10)
  DataCheck <- as.vector(assay(Data)[1:10,1])
  expect_equal(DataCheck0, DataCheck)

  # Checks isSpike info
  TechCheck0 <- rep(FALSE, 50)
  TechCheck <- rep(FALSE, nrow(Data))
  expect_equal(TechCheck0, TechCheck)
})
