context("Parameter estimation and denoised data (no-spikes+regression)\n")

test_that("Estimates match the given seed (no-spikes+regression)", {
  # Data example
  set.seed(13)
  Data <- makeExampleBASiCS_Data(WithSpikes = FALSE, WithBatch = TRUE)
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts(Data)),
    colData = data.frame(BatchInfo = SingleCellExperiment::colData(Data)$BatchInfo)
  )

  # Fixing starting values
  n <- ncol(Data)
  k <- 12
  PriorParam <- list(
    mu.mu = 0, s2.mu = 0.5, s2.delta = 0.5, a.delta = 1,
    b.delta = 1, p.phi = rep(1, times = n),
    GeneExponent = 1, CellExponent = 1,
    a.s = 1, b.s = 1, a.theta = 1, b.theta = 1
  )
  PriorParam$m <- rep(0, k)
  PriorParam$V <- diag(k)
  PriorParam$a.sigma2 <- 2
  PriorParam$b.sigma2 <- 2
  PriorParam$eta <- 5
  PriorParam$RBFNTile <- FALSE
  PriorParam$FixLocations <- FALSE
  PriorParam$locations <- rep(0, k)
  set.seed(2018)
  Start <- BASiCS:::HiddenBASiCS_MCMC_Start(Data, PriorParam,
                                            Regression = TRUE,
                                            WithSpikes = FALSE)

  # Running the sampler
  set.seed(14)
  Chain <- run_MCMC(Data, N = 1000, Thin = 10, Burn = 500,
                       PrintProgress = FALSE, WithSpikes = FALSE,
                       MinGenesPerRBF = NA,
                       Regression = TRUE)
  set.seed(14)
  ChainSCE <- run_MCMC(sce, N = 1000, Thin = 10, Burn = 500,
                       PrintProgress = FALSE, WithSpikes = FALSE,
                       MinGenesPerRBF = NA,
                       Regression = TRUE)

  # Calculating a posterior summary
  PostSummary <- Summary(Chain)
  PostSummarySCE <- Summary(ChainSCE)

  # Checking parameter names
  ParamNames <- c("mu", "delta", "s", "nu", "theta",
                  "beta", "sigma2", "epsilon", "RefFreq", "locations")
  ParamNames1 <- c("mu", "delta", "s", "nu", "theta",
                  "beta", "sigma2", "epsilon")
  expect_equal(names(Chain@parameters), ParamNames)
  expect_equal(names(PostSummary@parameters), ParamNames1)
            
  # Check if parameter estimates match for the first 5 genes and cells
  Mu <- c(13.927, 17.978, 5.653, 12.183, 37.181)
  MuObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "mu")[1:5, 1], 3))
  MuObsSCE <- as.vector(round(displaySummaryBASiCS(
    PostSummarySCE,
    "mu"
  )[1:5, 1], 3))
  expect_equal(MuObs, Mu, tolerance = 1, scale = 1)
  expect_equal(MuObsSCE, Mu, tolerance = 1, scale = 1)

  Delta <- c(1.487, 1.264, 1.998, 1.638, 0.472)
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

  S <- c(1.421, 0.916, 1.967, 0.812, 1.111)
  SObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "s")[1:5, 1], 3))
  SObsSCE <- as.vector(round(displaySummaryBASiCS(
    PostSummarySCE,
    "s"
  )[1:5, 1], 3))
  expect_equal(SObs, S, tolerance = 1, scale = 1)
  expect_equal(SObsSCE, S, tolerance = 1, scale = 1)

  Theta <- c(0.282, 0.143)
  ThetaObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "theta")[, 1], 3))
  ThetaObsSCE <- as.vector(round(displaySummaryBASiCS(
    PostSummarySCE,
    "theta"
  )[, 1], 3))
  expect_equal(ThetaObs, Theta, tolerance = 1, scale = 1)
  expect_equal(ThetaObsSCE, Theta, tolerance = 1, scale = 1)

  Beta <- c(0.302, -0.309, 0.303, 0.267, 0.063)
  BetaObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "beta")[1:5, 1], 3))
  BetaObsSCE <- as.vector(round(displaySummaryBASiCS(
    PostSummarySCE,
    "beta"
  )[1:5, 1], 3))
  expect_equal(BetaObs, Beta, tolerance = 1, scale = 1)
  expect_equal(BetaObsSCE, Beta, tolerance = 1, scale = 1)

  Sigma2 <- 0.232
  Sigma2Obs <- round(displaySummaryBASiCS(PostSummary, "sigma2")[1], 3)
  Sigma2ObsSCE <- round(displaySummaryBASiCS(PostSummarySCE, "sigma2")[1], 3)
  expect_equal(Sigma2Obs, Sigma2, tolerance = 1, scale = 1)
  expect_equal(Sigma2ObsSCE, Sigma2, tolerance = 1, scale = 1)

  # Obtaining denoised counts
  DC <- BASiCS_DenoisedCounts(Data, Chain)
  DCSCE <- BASiCS_DenoisedCounts(sce, ChainSCE)
  # Checks for an arbitrary set of genes / cells
  DCcheck0 <- c(4.527, 23.203, 7.357, 0.000, 19.241)
  DCcheck <- as.vector(round(DC[1:5, 1], 3))
  DCSCEcheck <- as.vector(round(DCSCE[1:5, 1], 3))
  expect_equal(DCcheck, DCcheck0, tolerance = 1, scale = 1)
  expect_equal(DCSCEcheck, DCcheck0, tolerance = 1, scale = 1)

  # Obtaining denoised rates
  DR <- BASiCS_DenoisedRates(Data, Chain)
  DRSCE <- BASiCS_DenoisedRates(sce, ChainSCE)
  # Checks for an arbitrary set of genes / cells
  DRcheck0 <- c(11.836, 19.464, 31.474, 32.466, 7.247)
  DRcheck <- as.vector(round(DR[10, 1:5], 3))
  DRSCEcheck <- as.vector(round(DRSCE[10, 1:5], 3))
  expect_equal(DRcheck, DRcheck0, tolerance = 1, scale = 1)
  expect_equal(DRSCEcheck, DRcheck0, tolerance = 1, scale = 1)
})

test_that("Chain creation works when regression, no spikes, and StoreAdapt=TRUE", {
  # Data example
  set.seed(14)
  Data <- makeExampleBASiCS_Data(WithSpikes = FALSE, WithBatch = TRUE)
  # Fixing starting values
  n <- ncol(Data)
  k <- 10
  PriorParam <- list(mu.mu = 0, s2.mu = 0.5, s2.delta = 0.5, a.delta = 1,
                     b.delta = 1, p.phi = rep(1, times = n),
                     a.s = 1, b.s = 1, a.theta = 1, b.theta = 1)
  PriorParam$m <- rep(0, k)
  PriorParam$V <- diag(k)
  PriorParam$a.sigma2 <- 2
  PriorParam$b.sigma2 <- 2
  PriorParam$eta <- 5
  PriorParam$RBFNTile <- FALSE
  PriorParam$FixLocations <- FALSE
  PriorParam$locations <- rep(0, k)
  set.seed(2018)
  Start <- BASiCS:::HiddenBASiCS_MCMC_Start(Data, PriorParam,
                                            Regression = TRUE,
                                            WithSpikes = FALSE)
  # Running the sampler
  set.seed(42)
  Chain <- run_MCMC(Data, N = 50, Thin = 10, Burn = 10,
                     PrintProgress = FALSE, WithSpikes = FALSE,
                     Regression = TRUE, StoreAdapt = TRUE,
                     MinGenesPerRBF = NA,
                     Start = Start, PriorParam = PriorParam)
  expect_s4_class(Chain, "BASiCS_Chain")
})
