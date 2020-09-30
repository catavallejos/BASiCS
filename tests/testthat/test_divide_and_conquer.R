context("Divide and conquer")

set.seed(42)
Data <- BASiCS_MockSCE()

test_that(".generateSubsets produces valid output by cell and gene with & w/o spikes", {
  check_subsets <- function(subsetby, nsubsets) {
    for (nsubset in nsubsets) {
      out <- .generateSubsets(Data, NSubsets=2, SubsetBy=subsetby, WithSpikes = TRUE)
      lapply(out, function(x) {
        expect_is(x, "SingleCellExperiment")
      })    
    }
  }
  check_subsets("cell", c(2, 32))
  check_subsets("gene", c(128))
})

test_that(".generateSubsets fails with low N", {
  DataS <- Data[, 1:20]
  expect_error(.generateSubsets(DataS, NSubsets=2, SubsetBy="cell"),
    "Cannot generate subsets with at least 20 cells per subset")
  expect_error(.generateSubsets(Data, NSubsets=8, SubsetBy="cell"),
    "Cannot generate subsets with at least 20 cells per subset")
})

test_that(".generateSubsets partitions well by gene", {
  l <- .generateSubsets(Data, SubsetBy = "gene", NSubsets = 128, WithSpikes = TRUE)
  expect_true(all(sapply(l, function(x) inherits(x, "SingleCellExperiment"))))
})


test_that("BASiCS_DivideAndConquer runs with standard settings", {
  expect_error(
    capture.output(
      m <- BASiCS_DivideAndConquer(
        Data,
        NSubsets = 2,
        SubsetBy = "gene",
        Regression = TRUE,
        PrintProgress = FALSE,
        WithSpikes = TRUE,
        BPParam = BiocParallel::SerialParam(),
        N = 8,
        Thin = 2,
        Burn = 4
      )
    ),
    NA
  )
})

test_that(".consensus_average produces sensible results (cell-wise)", {
  set.seed(42)

  capture.output(
    m <- BASiCS_DivideAndConquer(
      Data,
      NSubsets = 2,
      SubsetBy = "cell",
      Regression = TRUE,
      WithSpikes = TRUE,
      PrintProgress = FALSE,
      N = 50,
      BPParam = BiocParallel::SerialParam(),
      Thin = 5,
      Burn = 10
    )
  )
  capture.output({
    c1 <- BASiCS:::.consensus_average(m, Weighting="naive", SubsetBy = "cell", BPParam = BiocParallel::SerialParam(),)
    c2 <- BASiCS:::.consensus_average(m, Weighting="n_weight", SubsetBy = "cell", BPParam = BiocParallel::SerialParam(),)
  })
  m[[2]] <- BASiCS:::.offset_correct(m[[2]], m[[1]])
  gene <- "Gene 10"
  mean <- mean(
    c(
      m[[1]]@parameters[["mu"]][1, gene], 
      m[[2]]@parameters[["mu"]][1, gene]
    )
  )
  expect_equal(unname(c1@parameters[["mu"]][1, gene]),  mean)
  weighted_mean <- sum(
    c(
      m[[1]]@parameters[["mu"]][1, gene] * ncol(m[[1]]@parameters[["s"]]), 
      m[[2]]@parameters[["mu"]][1, gene] * ncol(m[[2]]@parameters[["s"]])
    )
  ) / ncol(Data)
  expect_equal(unname(c2@parameters[["mu"]][1, gene]),  weighted_mean)
})
