test_that("Generated data does not match given seed (all cases)", {
  
  # Data example
  set.seed(13)
  Data1 <- makeExampleBASiCS_Data(WithSpikes = FALSE, WithBatch = TRUE)
  Data2 <- makeExampleBASiCS_Data(WithSpikes = FALSE, WithBatch = FALSE)
  Data3 <- makeExampleBASiCS_Data(WithSpikes = TRUE, WithBatch = TRUE)
  Data4 <- makeExampleBASiCS_Data(WithSpikes = TRUE, WithBatch = FALSE)
  
  DataCheck01 <- c(2665,   1095,   3680,   1012,   1895)
  DataCheck02 <- c(1800,    823,   3705,   2969,   1132)
  DataCheck03 <- c( 242,    488,    429,    215,    159)
  DataCheck04 <- c( 528,    355,    123,    732,    154)
  
  DataCheck1 <- as.vector(colSums(assay(Data1)[,1:5]))
  DataCheck2 <- as.vector(colSums(assay(Data2)[,1:5]))
  DataCheck3 <- as.vector(colSums(assay(Data3)[,1:5]))
  DataCheck4 <- as.vector(colSums(assay(Data4)[,1:5]))
  
  expect_equal(DataCheck01, DataCheck1)
  expect_equal(DataCheck02, DataCheck2)
  expect_equal(DataCheck03, DataCheck3)
  expect_equal(DataCheck04, DataCheck4)
  
  # Specific checks for WithSpikes = TRUE cases
  # Checks isSpike info
  TechCheck0 <- c(rep(FALSE, 50), rep(TRUE, 20))
  TechCheck3 <- c(rep(FALSE, nrow(Data3)), rep(TRUE, nrow(altExp(Data3))))
  TechCheck4 <- c(rep(FALSE, nrow(Data4)), rep(TRUE, nrow(altExp(Data4))))
  expect_equal(TechCheck0, TechCheck3)
  expect_equal(TechCheck0, TechCheck4)
  
  # Checks total count for spike in genes
  TechCount03 <- c(11243,    80,  354,   688,    25)
  TechCount04 <- c(10280,    74,  354,   663,    26)
  TechCount3 <- as.vector(matrixStats::rowSums2(assay(altExp(Data3)[1:5,])))
  TechCount4 <- as.vector(matrixStats::rowSums2(assay(altExp(Data4)[1:5,])))
  expect_equal(TechCount03, TechCount3)
  expect_equal(TechCount04, TechCount4)
})

