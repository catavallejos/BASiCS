context("Divide and conquer")


set.seed(42)
mu <- rlnorm(100)
nu <- rgamma(100, 5, 5)
mat <- mu %*% t(nu)
counts <- matrix(
  rpois(length(mat), lambda = mat),
  nrow = nrow(mat),
  ncol = ncol(mat),
  dimnames = list(
    paste("Gene", seq_len(nrow(mat))),
    paste("Cell", seq_len(ncol(mat)))
  )
)
spikes <- rlnorm(20)
spike_mat <- spikes %*% t(nu)
spike_counts <- matrix(
  rpois(length(spike_mat), lambda = spike_mat),
  nrow = nrow(spike_mat),
  ncol = ncol(spike_mat),
  dimnames = list(
    paste("Spike", seq_len(nrow(spike_mat))),
    paste("Cell", seq_len(ncol(spike_mat)))
  )
)
Data <- SingleCellExperiment::SingleCellExperiment(
  assays = list(counts = counts)
)
colData(Data)$BatchInfo <- sample(1:2, nrow(sce), replace = TRUE)
Spikes <- SingleCellExperiment::SingleCellExperiment(
  assays = list(counts = spike_counts)
)
rowData(Spikes)[[1]] <- rownames(Spikes)
rowData(Spikes)[[2]] <- round(spikes) + 1
altExp(Data, "spike-ins") <- Spikes


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
      Thin = 5,
      Burn = 10
    )
  )
  capture.output({
    c1 <- BASiCS:::.consensus_average(m, weight_method="naive", subset_by = "cell")
    c2 <- BASiCS:::.consensus_average(m, weight_method="n_weight", subset_by = "cell")
  })
  m[[1]] <- BASiCS:::.offset_correct(m[[1]], m[[2]])
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
