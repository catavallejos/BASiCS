context("Parameter estimation and denoised data (no-spikes),
        original data doesn't have spikes\n")

test_that("Estimates match the given seed (no-spikes)",
{
  # Data example
  Data <- makeExampleBASiCS_Data(WithSpikes = FALSE,
                                 WithBatch = TRUE)
  sce <- SingleCellExperiment(assays = list(counts = counts(Data)),
                      colData = data.frame(BatchInfo = colData(Data)$BatchInfo))

  # Fixing starting values
  n <- ncol(Data)
  PriorParam <- list(s2.mu = 0.5, s2.delta = 0.5, a.delta = 1,
                     b.delta = 1, p.phi = rep(1, times = n),
                     a.s = 1, b.s = 1, a.theta = 1, b.theta = 1)
  set.seed(2018)
  Start <- BASiCS:::HiddenBASiCS_MCMC_Start(Data, PriorParam, WithSpikes = FALSE)
  # Running the sampler on Data and sce object
  set.seed(14)
  Chain <- run_BASiCS_MCMC(Data, N = 1000, Thin = 10, Burn = 500,
                       Regression = FALSE,
                       PrintProgress = FALSE, WithSpikes = FALSE)
  set.seed(14)
  ChainSCE <- run_BASiCS_MCMC(sce, N = 1000, Thin = 10, Burn = 500,
                       Regression = FALSE,
                       PrintProgress = FALSE, WithSpikes = FALSE)
  # Calculating a posterior summary
  PostSummary <- Summary(Chain)
  PostSummarySCE <- Summary(ChainSCE)

  # Checking parameter names
  ParamNames <- c("mu", "delta", "s", "nu", "theta", "RefFreq")
  ParamNames1 <- c("mu", "delta", "s", "nu", "theta")
  expect_equal(names(Chain@parameters), ParamNames)
  expect_equal(names(PostSummary@parameters), ParamNames1)

  # Check if parameter estimates match for the first 5 genes and cells
  Mu <- c(14.031, 11.601, 6.799, 10.882, 30.704)
  MuObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "mu")[1:5,1],3))
  MuObsSCE <- as.vector(round(displaySummaryBASiCS(PostSummarySCE,
                                                   "mu")[1:5,1],3))
  expect_equal(MuObs, Mu)
  expect_equal(MuObsSCE, Mu)

  Delta <- c(1.080, 1.125, 1.589, 1.272, 0.607)
  DeltaObs <- as.vector(round(displaySummaryBASiCS(PostSummary,
                                                   "delta")[1:5,1],3))
  DeltaObsSCE <- as.vector(round(displaySummaryBASiCS(PostSummarySCE,
                                                   "delta")[1:5,1],3))
  expect_equal(DeltaObs, Delta)
  expect_equal(DeltaObsSCE, Delta)

  S <- c(0.641, 1.277, 0.241, 0.631, 1.264)
  SObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "s")[1:5,1],3))
  SObsSCE <- as.vector(round(displaySummaryBASiCS(PostSummarySCE,
                                                  "s")[1:5,1],3))
  expect_equal(SObs, S)
  expect_equal(SObsSCE, S)

  Theta <- c(0.187, 0.313)
  ThetaObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "theta")[,1],3))
  ThetaObsSCE <- as.vector(round(displaySummaryBASiCS(PostSummarySCE,
                                                      "theta")[,1],3))
  expect_equal(ThetaObs, Theta)
  expect_equal(ThetaObsSCE, Theta)

  # Obtaining denoised counts
  set.seed(2018)
  DC <- BASiCS_DenoisedCounts(Data, Chain)
  set.seed(2018)
  DCSCE <- BASiCS_DenoisedCounts(sce, ChainSCE)
  # Checks for an arbitrary set of genes / cells
  DCcheck0 <- c(5.017, 1.672, 0.000, 10.033, 11.706)
  DCcheck <- as.vector(round(DC[1:5,1], 3))
  DCSCEcheck <- as.vector(round(DCSCE[1:5,1], 3))
  expect_equal(DCcheck, DCcheck0)
  expect_equal(DCSCEcheck, DCcheck0)

  # Obtaining denoised rates
  set.seed(2018)
  DR <- BASiCS_DenoisedRates(Data, Chain)
  set.seed(2018)
  DRSCE <- BASiCS_DenoisedRates(sce, ChainSCE)
  # Checks for an arbitrary set of genes / cells
  DRcheck0 <- c(5.786,  3.272, 10.572,  1.535,  4.880)
  DRcheck <- as.vector(round(DR[10,1:5], 3))
  DRSCEcheck <- as.vector(round(DRSCE[10,1:5], 3))
  expect_equal(DRcheck, DRcheck0)
  expect_equal(DRSCEcheck, DRcheck0)
})


test_that("Chain creation works when StoreAdapt=TRUE (no spikes)",
{
  # Data example
  Data <- makeExampleBASiCS_Data(WithSpikes = FALSE,
                                 WithBatch = TRUE)
  # Fixing starting values
  n <- ncol(Data)
  PriorParam <- list(s2.mu = 0.5, s2.delta = 0.5, a.delta = 1,
                     b.delta = 1, p.phi = rep(1, times = n),
                     a.s = 1, b.s = 1, a.theta = 1, b.theta = 1)
  set.seed(2018)
  Start <- BASiCS:::HiddenBASiCS_MCMC_Start(Data, PriorParam, WithSpikes = FALSE)
  # Running the sampler
  set.seed(14)
  Chain <- run_BASiCS_MCMC(Data, N = 50, Thin = 10, Burn = 10,
                       Regression = FALSE, StoreAdapt=TRUE,
                       PrintProgress = FALSE, WithSpikes = FALSE)
  expect_s4_class(Chain, "BASiCS_Chain")
})
