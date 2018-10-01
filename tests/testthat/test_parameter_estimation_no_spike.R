context("Parameter estimation and denoised data (no-spikes), 
        original data doesn't have spikes\n")

test_that("Estimates match the given seed (no-spikes)", 
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
  Chain <- BASiCS_MCMC(Data, N = 1000, Thin = 10, Burn = 500, 
                       Regression = FALSE,
                       PrintProgress = FALSE, WithSpikes = FALSE)
  # Calculating a posterior summary
  PostSummary <- Summary(Chain)
  
  # Checking parameter names
  ParamNames <- c("mu", "delta", "s", "nu", "theta", "RefFreq")
  ParamNames1 <- c("mu", "delta", "s", "nu", "theta")
  expect_that(all.equal(names(Chain@parameters), ParamNames), is_true())
  expect_that(all.equal(names(PostSummary@parameters), ParamNames1), is_true())
            
  # Check if parameter estimates match for the first 5 genes and cells
  Mu <- c(12.714, 9.361, 5.666, 10.267, 28.430)
  MuObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "mu")[1:5,1],3))
  expect_that(all.equal(MuObs, Mu), is_true())
            
  Delta <- c(1.369, 0.986, 1.650, 1.411, 0.660)
  DeltaObs <- as.vector(round(displaySummaryBASiCS(PostSummary, 
                                                   "delta")[1:5,1],3))
  expect_that(all.equal(DeltaObs, Delta), is_true())
  
  S <- c(0.901, 1.628, 0.319, 0.768, 1.504)
  SObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "s")[1:5,1],3))
  expect_that(all.equal(SObs, S), is_true())
  
  Theta <- c(0.320, 0.115)
  ThetaObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "theta")[,1],3))
  expect_that(all.equal(ThetaObs, Theta), is_true())
  
  # Obtaining denoised counts     
  DC <- BASiCS_DenoisedCounts(Data, Chain)
  # Checks for an arbitrary set of genes / cells
  DCcheck0 <- c(3.424, 1.141, 0.000, 6.847, 7.988)
  DCcheck <- as.vector(round(DC[1:5,1], 3))
  expect_that(all.equal(DCcheck, DCcheck0), is_true())
  
  # Obtaining denoised rates
  DR <- BASiCS_DenoisedRates(Data, Chain)
  # Checks for an arbitrary set of genes / cells
  DRcheck0 <- c(4.178,  2.461, 8.513,  1.517,  3.840)
  DRcheck <- as.vector(round(DR[10,1:5], 3))
  expect_that(all.equal(DRcheck, DRcheck0), is_true())
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
  Chain <- BASiCS_MCMC(Data, N = 50, Thin = 10, Burn = 10,
                       Regression = FALSE, StoreAdapt=TRUE,
                       PrintProgress = FALSE, WithSpikes = FALSE)
  expect_s4_class(Chain, "BASiCS_Chain")
})
