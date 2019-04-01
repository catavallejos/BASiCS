context("Parameter estimation and denoised data (no-spikes), 
        original data doesn't have spikes\n")

test_that("Estimates match the given seed (no-spikes)", 
{
  # Data example
  set.seed(10)
  Data <- makeExampleBASiCS_Data(WithSpikes = FALSE, 
                                 WithBatch = TRUE)
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts(Data)),
                      colData = data.frame(BatchInfo = SummarizedExperiment::colData(Data)$BatchInfo))
  
  # Fixing starting values
  n <- ncol(Data)
  PriorParam <- list(s2.mu = 0.5, s2.delta = 0.5, a.delta = 1, 
                     b.delta = 1, p.phi = rep(1, times = n), 
                     a.s = 1, b.s = 1, a.theta = 1, b.theta = 1)
  set.seed(2018)
  Start <- BASiCS:::HiddenBASiCS_MCMC_Start(Data, PriorParam, WithSpikes = FALSE)
  # Running the sampler on Data and sce object
  set.seed(14)
  Chain <- run_MCMC(Data, N = 1000, Thin = 10, Burn = 500, 
                       Regression = FALSE,
                       PrintProgress = FALSE, WithSpikes = FALSE)
  set.seed(14)
  ChainSCE <- run_MCMC(sce, N = 1000, Thin = 10, Burn = 500, 
                       Regression = FALSE,
                       PrintProgress = FALSE, WithSpikes = FALSE)
  # Calculating a posterior summary
  PostSummary <- Summary(Chain)
  PostSummarySCE <- Summary(ChainSCE)
  
  # Checking parameter names
  ParamNames <- c("mu", "delta", "s", "nu", "theta", "RefFreq")
  ParamNames1 <- c("mu", "delta", "s", "nu", "theta")
  expect_that(all.equal(names(Chain@parameters), ParamNames), is_true())
  expect_that(all.equal(names(PostSummary@parameters), ParamNames1), is_true())
            
  # Check if parameter estimates match for the first 5 genes and cells
  Mu <- c(9.625, 14.665,  6.997,  8.724, 31.469)
  MuObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "mu")[1:5,1],3))
  MuObsSCE <- as.vector(round(displaySummaryBASiCS(PostSummarySCE, 
                                                   "mu")[1:5,1],3))
  expect_that(all.equal(MuObs, Mu), is_true())
  expect_that(all.equal(MuObsSCE, Mu), is_true())
            
  Delta <- c(1.234, 0.949, 1.710, 1.414, 0.440)
  DeltaObs <- as.vector(round(displaySummaryBASiCS(PostSummary, 
                                                   "delta")[1:5,1],3))
  DeltaObsSCE <- as.vector(round(displaySummaryBASiCS(PostSummarySCE, 
                                                   "delta")[1:5,1],3))
  expect_that(all.equal(DeltaObs, Delta), is_true())
  expect_that(all.equal(DeltaObsSCE, Delta), is_true())
  
  S <- c(1.387, 1.552, 0.610, 2.184, 1.457)
  SObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "s")[1:5,1],3))
  SObsSCE <- as.vector(round(displaySummaryBASiCS(PostSummarySCE, 
                                                  "s")[1:5,1],3))
  expect_that(all.equal(SObs, S), is_true())
  expect_that(all.equal(SObsSCE, S), is_true())
  
  Theta <- c(0.120, 0.109)
  ThetaObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "theta")[,1],3))
  ThetaObsSCE <- as.vector(round(displaySummaryBASiCS(PostSummarySCE, 
                                                      "theta")[,1],3))
  expect_that(all.equal(ThetaObs, Theta), is_true())
  expect_that(all.equal(ThetaObsSCE, Theta), is_true())
  
  # Obtaining denoised counts   
  set.seed(2018)
  DC <- BASiCS_DenoisedCounts(Data, Chain)
  set.seed(2018)
  DCSCE <- BASiCS_DenoisedCounts(sce, ChainSCE)
  # Checks for an arbitrary set of genes / cells
  DCcheck0 <- c(31.007, 21.633,  7.211,  3.605, 63.456)
  DCcheck <- as.vector(round(DC[1:5,1], 3))
  DCSCEcheck <- as.vector(round(DCSCE[1:5,1], 3))
  expect_that(all.equal(DCcheck, DCcheck0), is_true())
  expect_that(all.equal(DCSCEcheck, DCcheck0), is_true())
  
  # Obtaining denoised rates
  set.seed(2018)
  DR <- BASiCS_DenoisedRates(Data, Chain)
  set.seed(2018)
  DRSCE <- BASiCS_DenoisedRates(sce, ChainSCE)
  # Checks for an arbitrary set of genes / cells
  DRcheck0 <- c(2.193,  2.981, 19.924, 15.005,  5.930)
  DRcheck <- as.vector(round(DR[10,1:5], 3))
  DRSCEcheck <- as.vector(round(DRSCE[10,1:5], 3))
  expect_that(all.equal(DRcheck, DRcheck0), is_true())
  expect_that(all.equal(DRSCEcheck, DRcheck0), is_true())
})


test_that("Chain creation works when StoreAdapt=TRUE (no spikes)", 
{
  # Data example
  set.seed(11)
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
  Chain <- run_MCMC(Data, N = 50, Thin = 10, Burn = 10,
                       Regression = FALSE, StoreAdapt=TRUE,
                       PrintProgress = FALSE, WithSpikes = FALSE)
  expect_s4_class(Chain, "BASiCS_Chain")
})
