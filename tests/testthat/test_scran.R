context("Testing scran::computeSumFactors")

test_that("scran::computeSumFactors changed for shallow data -- talk to Aaron!", {

  # Shallow data
  set.seed(1)
  X <- matrix(rpois(100*1000, 10), ncol = 100, nrow = 1000)
  ColSumX <- colSums(X)
  RowSumX <- rowSums(X)

  Check0 <- c(9776, 10124, 10119,  9986, 10143)
  Check <- as.vector(ColSumX[1:5])
  expect_equal(Check0, Check)

  Check0 <- c(1020, 1033,  975,  992, 1025)
  Check <- as.vector(RowSumX[1:5])
  expect_equal(Check0, Check)

  sce <- SingleCellExperiment(assays=list(counts = X))
  sce <- scran::computeSumFactors(sce)

  Check0 <- c(0.974, 1.010, 1.005, 1.012, 1.018)
  Check <- round(sizeFactors(sce)[1:5], 3)
  expect_equal(Check0, Check)
})

test_that("scran::computeSumFactors changed for deep data -- talk to Aaron!", {

  # Shallow data
  set.seed(1)
  X <- matrix(rpois(100*1000, 200), ncol = 100, nrow = 1000)
  ColSumX <- colSums(X)
  RowSumX <- rowSums(X)

  Check0 <- c(199570, 200028, 200546, 199463, 200187)
  Check <- as.vector(ColSumX[1:5])
  expect_equal(Check0, Check)

  Check0 <- c(20093, 19966, 20058, 20060, 19854)
  Check <- as.vector(RowSumX[1:5])
  expect_equal(Check0, Check)

  sce <- SingleCellExperiment(assays=list(counts = X))
  sce <- scran::computeSumFactors(sce)

  Check0 <- c(0.998, 0.997, 1.001, 0.999, 1.003)
  Check <- round(sizeFactors(sce)[1:5], 3)
  expect_equal(Check0, Check)
}) 
