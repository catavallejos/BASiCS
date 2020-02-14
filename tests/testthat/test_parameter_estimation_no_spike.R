context("Parameter estimation and denoised data (no-spikes), 
        original data doesn't have spikes\n")

test_that("Estimates match the given seed (no-spikes)", {
  # Data example
  set.seed(10)
  Data <- makeExampleBASiCS_Data(
    WithSpikes = FALSE,
    WithBatch = TRUE
  )
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts(Data)),
    colData = data.frame(
      BatchInfo = SummarizedExperiment::colData(Data)$BatchInfo
    )
  )

  # Fixing starting values
  n <- ncol(Data)
  PriorParam <- BASiCS_PriorParam(Data, k = 12)
  set.seed(2018)
  Start <- BASiCS:::HiddenBASiCS_MCMC_Start(
    Data,
    PriorParam,
    WithSpikes = FALSE,
    Regression = FALSE
  )
  # Running the sampler on Data and sce object
  set.seed(14)
  Chain <- run_MCMC(Data,
    N = 1000,
    Thin = 10,
    Burn = 500,
    Regression = FALSE,
    k = 12,
    PrintProgress = FALSE,
    WithSpikes = FALSE
  )
  set.seed(14)
  ChainSCE <- run_MCMC(sce,
    N = 1000,
    Thin = 10,
    Burn = 500,
    Regression = FALSE,
    k = 12,
    PrintProgress = FALSE,
    WithSpikes = FALSE
  )
  # Calculating a posterior summary
  PostSummary <- Summary(Chain)
  PostSummarySCE <- Summary(ChainSCE)

  # Checking parameter names
  ParamNames <- c("mu", "delta", "s", "nu", "theta", "RefFreq")
  ParamNames1 <- c("mu", "delta", "s", "nu", "theta")
  expect_equal(names(Chain@parameters), ParamNames, tolerance = 1, scale = 1)
  expect_equal(names(PostSummary@parameters), ParamNames1, tolerance = 1, scale = 1)

  # Check if parameter estimates match for the first 5 genes and cells
  Mu <- c(9.625, 14.665, 6.997, 8.724, 31.469)
  MuObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "mu")[1:5, 1], 3))
  MuObsSCE <- as.vector(round(displaySummaryBASiCS(
    PostSummarySCE,
    "mu"
  )[1:5, 1], 3))
  expect_equal(MuObs, Mu, tolerance = 1, scale = 1)
  expect_equal(MuObsSCE, Mu, tolerance = 1, scale = 1)

  Delta <- c(1.234, 0.949, 1.710, 1.414, 0.440)
  DeltaObs <- as.vector(round(displaySummaryBASiCS(
    PostSummary,
    "delta"
  )[1:5, 1], 3))
  DeltaObsSCE <- as.vector(round(displaySummaryBASiCS(
    PostSummarySCE,
    "delta"
  )[1:5, 1], 3))
  expect_equal(DeltaObs, Delta, tolerance = 1, scale = 1)
  expect_equal(DeltaObsSCE, Delta, tolerance = 1, scale = 1)

  S <- c(1.387, 1.552, 0.610, 2.184, 1.457)
  SObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "s")[1:5, 1], 3))
  SObsSCE <- as.vector(round(displaySummaryBASiCS(
    PostSummarySCE,
    "s"
  )[1:5, 1], 3))
  expect_equal(SObs, S, tolerance = 1, scale = 1)
  expect_equal(SObsSCE, S, tolerance = 1, scale = 1)

  Theta <- c(0.120, 0.109)
  ThetaObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "theta")[, 1], 3))
  ThetaObsSCE <- as.vector(round(displaySummaryBASiCS(
    PostSummarySCE,
    "theta"
  )[, 1], 3))
  expect_equal(ThetaObs, Theta, tolerance = 1, scale = 1)
  expect_equal(ThetaObsSCE, Theta, tolerance = 1, scale = 1)

  # Obtaining denoised counts
  set.seed(2018)
  DC <- BASiCS_DenoisedCounts(Data, Chain)
  set.seed(2018)
  DCSCE <- BASiCS_DenoisedCounts(sce, ChainSCE)
  # Checks for an arbitrary set of genes / cells
  DCcheck0 <- c(31.007, 21.633, 7.211, 3.605, 63.456)
  DCcheck <- as.vector(round(DC[1:5, 1], 3))
  DCSCEcheck <- as.vector(round(DCSCE[1:5, 1], 3))
  expect_equal(DCcheck, DCcheck0, tolerance = 1, scale = 5)
  expect_equal(DCSCEcheck, DCcheck0, tolerance = 1, scale = 1)

  # Obtaining denoised rates
  set.seed(2018)
  DR <- BASiCS_DenoisedRates(Data, Chain)
  set.seed(2018)
  DRSCE <- BASiCS_DenoisedRates(sce, ChainSCE)
  # Checks for an arbitrary set of genes / cells
  DRcheck0 <- c(2.193, 2.981, 19.924, 15.005, 5.930)
  DRcheck <- as.vector(round(DR[10, 1:5], 3))
  DRSCEcheck <- as.vector(round(DRSCE[10, 1:5], 3))
  expect_equal(DRcheck, DRcheck0, tolerance = 1, scale = 2)
  expect_equal(DRSCEcheck, DRcheck0, tolerance = 1, scale = 2)
})


test_that("Chain creation works when StoreAdapt=TRUE (no spikes)", {
  # Data example
  set.seed(11)
  Data <- makeExampleBASiCS_Data(
    WithSpikes = FALSE,
    WithBatch = TRUE
  )
  # Fixing starting values
  PriorParam <- BASiCS_PriorParam(Data, k = 12)
  set.seed(2018)
  Start <- BASiCS:::HiddenBASiCS_MCMC_Start(
    Data,
    PriorParam,
    WithSpikes = FALSE,
    Regression = FALSE
  )
  # Running the sampler
  set.seed(14)
  Chain <- run_MCMC(Data,
    N = 8, Thin = 2, Burn = 4,
    Regression = FALSE, StoreAdapt = TRUE,
    PrintProgress = FALSE, WithSpikes = FALSE
  )
  expect_s4_class(Chain, "BASiCS_Chain")
})
