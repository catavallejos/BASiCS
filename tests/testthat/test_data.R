context("data")

test_that("data preserves colnames", {
  # Expression counts
  set.seed(1)
  Counts <- matrix(rpois(50*10, 2), ncol = 10)
  rownames(Counts) <- paste0('Gene', 1:50)
  colnames(Counts) <- paste('Cell', seq_len(ncol(Counts)))
  BatchInfo <- c(rep(1, 5), rep(2, 5))
  # Creating a BASiCS_Data object (with batch effect)
  expect_warning(
    DataExample <- newBASiCS_Data(Counts, BatchInfo = BatchInfo)
  )
  expect_identical(colnames(DataExample), colnames(Counts))
})
