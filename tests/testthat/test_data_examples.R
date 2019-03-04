context("Data examples\n")

test_that("Generated data does not match given seed (spikes case)", {
  Data <- makeExampleBASiCS_Data(WithSpikes = TRUE)
  
  # Checks some specific cells
  DataCheck0 <- c(0, 0, 0, 1, 1, 3, 17, 18, 1, 0)
  DataCheck <- as.vector(assay(Data)[1:10,1])
  expect_equal(DataCheck0, DataCheck)
  
  # Checks isSpike info
  TechCheck0 <- c(rep(FALSE, 50), rep(TRUE, 20))
  TechCheck <- isSpike(Data, "ERCC")
  expect_equal(TechCheck0, TechCheck)
  
  # Checks total count for spike in genes
  TechCount0 <- c( 11213,    87, 355,  707,   18, 2807, 363, 336, 328, 711,
                   44356, 44573,  48, 1449, 2757,   87, 175, 1343,  35, 195)
  TechCount <- matrixStats::rowSums2(assay(Data)[isSpike(Data, "ERCC"),])
  expect_equal(TechCount0, TechCount)
})

test_that("Generated data does not match given seed (no spikes case)", {
  
  # Checks some specific cells
  Data <- makeExampleBASiCS_Data(WithSpikes = FALSE)
  DataCheck0 <- c(3, 1, 0, 6, 7, 15, 100, 84, 2, 3)
  DataCheck <- as.vector(assay(Data)[1:10,1])
  expect_equal(DataCheck0, DataCheck)
  
  # Checks isSpike info
  TechCheck0 <- rep(FALSE, 50)
  TechCheck <- isSpike(Data)
  expect_equal(TechCheck0, TechCheck)
})

