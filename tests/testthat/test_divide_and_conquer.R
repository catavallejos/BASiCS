context("Divide and conquer")

set.seed(42)
Data <- BASiCS_MockSCE()

bp <- BiocParallel::SerialParam()
BiocParallel::register(bp)

test_that("BASiCS:::.generateSubsets produces valid output by cell and gene with & w/o spikes", {
  check_subsets <- function(subsetby, nsubsets) {
    for (nsubset in nsubsets) {
      out <- BASiCS:::.generateSubsets(Data, NSubsets=2, SubsetBy=subsetby, WithSpikes = TRUE)
      lapply(out, function(x) {
        expect_is(x, "SingleCellExperiment")
      })    
    }
  }
  check_subsets("cell", c(2, 32))
  check_subsets("gene", c(128))
})

test_that("BASiCS:::.generateSubsets fails with low N", {
  DataS <- Data[, 1:20]
  expect_error(BASiCS:::.generateSubsets(DataS, NSubsets=2, SubsetBy="cell"),
    "Cannot generate subsets with at least 20 cells per subset")
  expect_error(BASiCS:::.generateSubsets(Data, NSubsets=8, SubsetBy="cell"),
    "Cannot generate subsets with at least 20 cells per subset")
})

test_that("BASiCS:::.generateSubsets partitions well by gene", {
  l <- BASiCS:::.generateSubsets(Data, SubsetBy = "gene", NSubsets = 128, WithSpikes = TRUE)
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
        BPPARAM = bp,
        N = 8,
        Thin = 2,
        Burn = 4
      )
    ),
    NA
  )
})


test_that("BASiCS_MCMC runs with divide and conquer", {
  for (x in c("gene", "cell")) {
    fun <- if (x == "gene") function(x) expect_error(x, NA) else expect_warning
    set.seed(1)
    fun(
      run_MCMC(
        Data,
        NSubsets = 2,
        SubsetBy = x,
        Regression = TRUE,
        PrintProgress = FALSE,
        WithSpikes = TRUE,
        BPPARAM = bp,
        N = 8,
        Thin = 2,
        Burn = 4
      )
    )
    set.seed(2)
    fun(
      run_MCMC(
        Data,
        NSubsets = 2,
        SubsetBy = x,
        Regression = FALSE,
        PrintProgress = FALSE,
        WithSpikes = TRUE,
        BPPARAM = bp,
        N = 8,
        Thin = 2,
        Burn = 4
      )
    )
    set.seed(3)
    fun(
      run_MCMC(
        Data,
        NSubsets = 2,
        SubsetBy = x,
        Regression = TRUE,
        PrintProgress = FALSE,
        WithSpikes = FALSE,
        BPPARAM = bp,
        N = 8,
        Thin = 2,
        Burn = 4
      )
    )
    set.seed(4)
    fun(
      run_MCMC(
        Data,
        NSubsets = 2,
        SubsetBy = x,
        Regression = FALSE,
        PrintProgress = FALSE,
        WithSpikes = FALSE,
        BPPARAM = bp,
        N = 8,
        Thin = 2,
        Burn = 4
      )
    )
  }
})

test_that("NSubsets check is appropriate", {
  expect_error(
    run_MCMC(
      Data,
      NSubsets = 1,
      SubsetBy = "gene",
      Regression = FALSE,
      PrintProgress = FALSE,
      WithSpikes = FALSE,
      BPPARAM = bp,
      N = 8,
      Thin = 2,
      Burn = 4
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
      BPPARAM = bp,
      Thin = 5,
      Burn = 10
    )
  )
  capture.output({
    c1 <- BASiCS:::.consensus_average(m, Weighting="naive", SubsetBy = "cell", BPPARAM = bp)
    c2 <- BASiCS:::.consensus_average(m, Weighting="n_weight", SubsetBy = "cell", BPPARAM = bp)
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


test_that("combine_subposteriors with NA", {
  Chain2 <- Chain1 <- ChainSC
  Chain1@parameters$delta <- Chain1@parameters$delta[, 1:10]
  Chain1@parameters$mu <- Chain1@parameters$mu[, 1:10]

  Chain2@parameters$mu <- Chain2@parameters$mu[, 1:10]
  Chain1@parameters$delta[, 1] <- NA
  expect_error(
    .combine_subposteriors(list(Chain1, Chain2), SubsetBy = "gene"),
    "Too many draws for parameter"
  )
  Chain2@parameters$delta <- Chain2@parameters$delta[, 11:20]
  Chain2@parameters$mu <- ChainSC@parameters$mu[, 11:20]
  expect_error(
    .combine_subposteriors(list(Chain1, Chain2), SubsetBy = "gene"),
    NA
  )
})


test_that("cut doesn't fail if some quantiles are the same", {
  c <- counts(Data)
  for (i in 1:20) c[i, ] <- c[1, ]
  counts(Data) <- c
  expect_message(
    capture.output(
      m <- BASiCS_DivideAndConquer(
        Data,
        NSubsets = 2,
        SubsetBy = "gene",
        Regression = TRUE,
        WithSpikes = TRUE,
        PrintProgress = FALSE,
        N = 50,
        BPPARAM = bp,
        Thin = 5,
        Burn = 10
      )
    ), "Cannot find a balanced split"
  )
})
