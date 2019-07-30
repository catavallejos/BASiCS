context("Parameter estimation and denoised data (no-spikes+regression)\n")

test_that("Estimates match the given seed (no-spikes+regression)", 
{
  # Data example
  set.seed(13)
  Data <- makeExampleBASiCS_Data(WithSpikes = FALSE, WithBatch = TRUE)
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts(Data)),
                      colData = data.frame(BatchInfo = SingleCellExperiment::colData(Data)$BatchInfo))
  
  # Fixing starting values
  n <- ncol(Data)
  k <- 12
  PriorParam <- list(s2.mu = 0.5, s2.delta = 0.5, a.delta = 1,
                     b.delta = 1, p.phi = rep(1, times = n),
                     a.s = 1, b.s = 1, a.theta = 1, b.theta = 1)
  PriorParam$m <- rep(0, k); PriorParam$V <- diag(k)
  PriorParam$a.sigma2 <- 2; PriorParam$b.sigma2 <- 2  
  PriorParam$eta <- 5
  set.seed(2018)
  Start <- BASiCS:::HiddenBASiCS_MCMC_Start(Data, PriorParam,
                                            WithSpikes = FALSE)

  # Running the sampler
  set.seed(14)
  Chain <- run_MCMC(Data, N = 1000, Thin = 10, Burn = 500,
                       PrintProgress = FALSE, WithSpikes = FALSE,
                       Regression = TRUE)
  set.seed(14)
  ChainSCE <- run_MCMC(sce, N = 1000, Thin = 10, Burn = 500,
                       PrintProgress = FALSE, WithSpikes = FALSE,
                       Regression = TRUE)

  # Calculating a posterior summary
  PostSummary <- Summary(Chain)
  PostSummarySCE <- Summary(ChainSCE)
  
  # Checking parameter names
  ParamNames <- c("mu", "delta", "s", "nu", "theta",
                  "beta", "sigma2", "epsilon", "RefFreq")
  ParamNames1 <- c("mu", "delta", "s", "nu", "theta",
                  "beta", "sigma2", "epsilon")
  expect_true(all.equal(names(Chain@parameters), ParamNames))
  expect_true(all.equal(names(PostSummary@parameters), ParamNames1))
            
  # Check if parameter estimates match for the first 5 genes and cells
  Mu <- c(14.266, 16.181,  5.343, 11.540, 41.315)
  MuObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "mu")[1:5,1],3))
  MuObsSCE <- as.vector(round(displaySummaryBASiCS(PostSummarySCE, 
                                                "mu")[1:5,1],3))
  expect_true(all.equal(MuObs, Mu))
  expect_true(all.equal(MuObsSCE, Mu))
            
  Delta <- c(1.286, 1.354, 1.783, 1.677, 0.493)
  DeltaObs <- as.vector(round(displaySummaryBASiCS(PostSummary, 
                                                   "delta")[1:5,1],3))
  DeltaObsSCE <- as.vector(round(displaySummaryBASiCS(PostSummarySCE, 
                                                   "delta")[1:5,1],3))
  expect_true(all.equal(DeltaObs, Delta))
  expect_true(all.equal(DeltaObsSCE, Delta))
            
  S <- c(1.653, 0.848, 2.096, 0.687, 1.143)
  SObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "s")[1:5,1],3))
  SObsSCE <- as.vector(round(displaySummaryBASiCS(PostSummarySCE, 
                                                  "s")[1:5,1],3))
  expect_true(all.equal(SObs, S))
  expect_true(all.equal(SObsSCE, S))
  
  Theta <- c(0.088, 0.073)
  ThetaObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "theta")[,1],3))
  ThetaObsSCE <- as.vector(round(displaySummaryBASiCS(PostSummarySCE, 
                                                   "theta")[,1],3))
  expect_true(all.equal(ThetaObs, Theta))
  expect_true(all.equal(ThetaObsSCE, Theta))
  
  Beta <- c(0.240, -0.253,  0.410,  0.281,  0.078)
  BetaObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "beta")[1:5,1],3))
  BetaObsSCE <- as.vector(round(displaySummaryBASiCS(PostSummarySCE, 
                                                  "beta")[1:5,1],3))
  expect_true(all.equal(BetaObs, Beta))
  expect_true(all.equal(BetaObsSCE, Beta))
  
  Sigma2 <- 0.239
  Sigma2Obs <- round(displaySummaryBASiCS(PostSummary, "sigma2")[1],3)
  Sigma2ObsSCE <- round(displaySummaryBASiCS(PostSummarySCE, "sigma2")[1],3)
  expect_true(all.equal(Sigma2Obs, Sigma2))
  expect_true(all.equal(Sigma2ObsSCE, Sigma2))
  
  # Obtaining denoised counts     
  DC <- BASiCS_DenoisedCounts(Data, Chain)
  DCSCE <- BASiCS_DenoisedCounts(sce, ChainSCE)
  # Checks for an arbitrary set of genes / cells
  DCcheck0 <- c(4.867, 24.941,  7.908,  0.000, 20.683)
  DCcheck <- as.vector(round(DC[1:5,1], 3))
  DCSCEcheck <- as.vector(round(DCSCE[1:5,1], 3))
  expect_true(all.equal(DCcheck, DCcheck0))
  expect_true(all.equal(DCSCEcheck, DCcheck0))
  
  # Obtaining denoised rates
  DR <- BASiCS_DenoisedRates(Data, Chain)
  DRSCE <- BASiCS_DenoisedRates(sce, ChainSCE)
  # Checks for an arbitrary set of genes / cells
  DRcheck0 <- c(12.137, 19.305, 34.570, 34.671,  7.422)
  DRcheck <- as.vector(round(DR[10,1:5], 3))
  DRSCEcheck <- as.vector(round(DRSCE[10,1:5], 3))
  expect_true(all.equal(DRcheck, DRcheck0))
  expect_true(all.equal(DRSCEcheck, DRcheck0))
})

test_that("Chain creation works when regression, no spikes, and StoreAdapt=TRUE", {
  # Data example
  set.seed(14)
  Data <- makeExampleBASiCS_Data(WithSpikes = FALSE, WithBatch = TRUE)
  # Fixing starting values
  n <- ncol(Data)
  k <- 12
  PriorParam <- list(s2.mu = 0.5, s2.delta = 0.5, a.delta = 1,
                     b.delta = 1, p.phi = rep(1, times = n),
                     a.s = 1, b.s = 1, a.theta = 1, b.theta = 1)
  PriorParam$m <- rep(0, k); PriorParam$V <- diag(k)
  PriorParam$a.sigma2 <- 2; PriorParam$b.sigma2 <- 2
  PriorParam$eta <- 5
  set.seed(2018)
  Start <- BASiCS:::HiddenBASiCS_MCMC_Start(Data, PriorParam, 
                                            WithSpikes = FALSE)
  # Running the sampler
  set.seed(42)
  Chain <- run_MCMC(Data, N = 50, Thin = 10, Burn = 10,
                     PrintProgress = FALSE, WithSpikes = FALSE,
                     Regression = TRUE, StoreAdapt = TRUE,
                     Start = Start, PriorParam = PriorParam)
  expect_s4_class(Chain, "BASiCS_Chain")
})
