context("Data examples\n")

test_that("Generated data does not match given seed (spikes case)", {
  set.seed(1)
  Data <- makeExampleBASiCS_Data(WithSpikes = TRUE)
  
  # Checks some specific cells
  DataCheck0 <- c(0,  5,  1,  0,  0,  1, 17, 11,  1,  2)
  DataCheck <- as.vector(assay(Data)[1:10,1])
<<<<<<< HEAD
  expect_true(all.equal(DataCheck0, DataCheck))
=======
  expect_equal(DataCheck0, DataCheck)
>>>>>>> 5d23e114161e35201f17adb88647424ac38c9273
  
  # Checks isSpike info
  TechCheck0 <- c(rep(FALSE, 50), rep(TRUE, 20))
  TechCheck <- isSpike(Data, "ERCC")
<<<<<<< HEAD
  expect_true(all.equal(TechCheck0, TechCheck))
  
  # Checks total count for spike in genes
  TechCount0 <- c( 11087,   101,   344,   633,    20)
  TechCount <- matrixStats::rowSums2(assay(Data)[isSpike(Data, "ERCC"),])[1:5]
  expect_true(all.equal(TechCount0, TechCount))
=======
  expect_equal(TechCheck0, TechCheck)
  
  # Checks total count for spike in genes
  TechCount0 <- c( 11213,    87, 355,  707,   18, 2807, 363, 336, 328, 711,
                   44356, 44573,  48, 1449, 2757,   87, 175, 1343,  35, 195)
  TechCount <- matrixStats::rowSums2(assay(Data)[isSpike(Data, "ERCC"),])
  expect_equal(TechCount0, TechCount)
>>>>>>> 5d23e114161e35201f17adb88647424ac38c9273
})

test_that("Generated data does not match given seed (no spikes case)", {
  
  # Checks some specific cells
  set.seed(2)
  Data <- makeExampleBASiCS_Data(WithSpikes = FALSE)
  DataCheck0 <- c(0,   5,   3,   6,  28,   4,   3, 170,   0,  10)
  DataCheck <- as.vector(assay(Data)[1:10,1])
<<<<<<< HEAD
  expect_true(all.equal(DataCheck0, DataCheck))
=======
  expect_equal(DataCheck0, DataCheck)
>>>>>>> 5d23e114161e35201f17adb88647424ac38c9273
  
  # Checks isSpike info
  TechCheck0 <- rep(FALSE, 50)
  TechCheck <- isSpike(Data)
<<<<<<< HEAD
  expect_true(all.equal(TechCheck0, TechCheck))
=======
  expect_equal(TechCheck0, TechCheck)
>>>>>>> 5d23e114161e35201f17adb88647424ac38c9273
})

