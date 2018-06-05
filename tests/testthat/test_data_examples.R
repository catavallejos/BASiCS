context("Data examples\n")

test_that("Generated data does not match given seed (spikes case)", {
  Data <- makeExampleBASiCS_Data(WithSpikes = TRUE)
  DataCheck0 <- c(0, 0, 0, 1, 1, 3, 17, 18, 1, 0)
  DataCheck <- as.vector(assay(Data)[1:10,1])
  expect_that(all.equal(DataCheck0, DataCheck), is_true())
  
  # Checks isSpike info
  TechCheck0 <- c(rep(FALSE, 50), rep(TRUE, 20))
  TechCheck <- isSpike(Data, "ERCC")
  expect_that(all.equal(TechCheck0, TechCheck), is_true())
  
  # Checks total count for spike in genes
  TechCount0 <- c( 8181,    64, 257,  511,   15, 2013, 279, 234, 244, 507,
                   31945, 32292,  32, 1041, 1986,   65, 124, 957,  27, 136)
  TechCount <- matrixStats::rowSums2(assay(Data)[isSpike(Data, "ERCC"),])
  expect_that(all.equal(TechCount0, TechCount), is_true())
})

test_that("Generated data does not match given seed (no spikes case)", {
  Data <- makeExampleBASiCS_Data(WithSpikes = FALSE)
  DataCheck0 <- c(3, 1, 0, 6, 7, 15, 100, 84, 2, 3)
  DataCheck <- as.vector(assay(Data)[1:10,1])
  expect_that(all.equal(DataCheck0, DataCheck), is_true())
  
  # Checks isSpike info
  TechCheck0 <- rep(FALSE, 50)
  TechCheck <- isSpike(Data)
  expect_that(all.equal(TechCheck0, TechCheck), is_true())
})
