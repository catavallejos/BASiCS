context("Parameter estimation and denoised data (spikes) \n")

test_that("Estimates match the given seed (spikes)", {
  # Data example
  set.seed(7)
  Data <- makeExampleBASiCS_Data(WithSpikes = TRUE)
  # Fixing starting values
  n <- ncol(Data)
  PriorParam <- list(mu.mu = rep(0, times = nrow(Data)), s2.mu = 0.5, 
                     s2.delta = 0.5, a.delta = 1, 
                     b.delta = 1, p.phi = rep(1, times = n), 
                     a.s = 1, b.s = 1, a.theta = 1, b.theta = 1)
  set.seed(2018)
  Start <- BASiCS:::HiddenBASiCS_MCMC_Start(Data, PriorParam, WithSpikes = TRUE)
  # Running the samples
  set.seed(18)
  Chain <- run_MCMC(Data, N = 1000, Thin = 10, Burn = 500, 
                       Regression = FALSE, PrintProgress = FALSE, 
                       Start = Start, PriorParam = PriorParam)
  # Calculating a posterior summary
  PostSummary <- Summary(Chain)
  
  # Checking parameter names
  ParamNames <- c("mu", "delta", "phi", "s", "nu", "theta")
  expect_equal(names(Chain@parameters), ParamNames)
  expect_equal(names(PostSummary@parameters), ParamNames)
            
  # Check if parameter estimates match for the first 5 genes and cells
  Mu <- c(9.983,  6.903,  3.242,  5.589, 23.492)
  MuObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "mu")[1:5,1],3))
  expect_equal(MuObs, Mu, tolerance = 1, scale = 1)
            
  Delta <- c(0.870, 0.731, 1.614, 1.496, 0.472)
  DeltaObs <- as.vector(round(displaySummaryBASiCS(PostSummary, 
                                                   "delta")[1:5,1],3))
  expect_equal(DeltaObs, Delta, tolerance = 1, scale = 1)

  Phi <- c(0.998, 0.682, 1.131, 1.146, 0.859)
  PhiObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "phi")[1:5,1],3))
  expect_equal(PhiObs, Phi, tolerance = 1, scale = 1)
            
  S <- c(1.017, 0.114, 0.606, 1.095, 0.289)
  SObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "s")[1:5,1],3))
  expect_equal(SObs, S, tolerance = 1, scale = 1)

  Theta <- 0.251
  ThetaObs <- round(displaySummaryBASiCS(PostSummary, "theta")[1],3)
  expect_equal(ThetaObs, Theta, tolerance = 1, scale = 1)
  
  # Obtaining denoised counts     
  DC <- BASiCS_DenoisedCounts(Data, Chain)
  # Checks for an arbitrary set of genes / cells
  DCcheck0 <- c(22.559,  0.940,  0.000,  1.880, 27.259)
  DCcheck <- as.vector(round(DC[1:5,1], 3))
  expect_equal(DCcheck, DCcheck0, tolerance = 1.5, scale = 1)
  
  # Obtaining denoised rates
  DR <- BASiCS_DenoisedRates(Data, Chain)
  # Checks for an arbitrary set of genes / cells
  DRcheck0 <- c(0.503, 2.591, 7.458, 4.614, 1.592)
  DRcheck <- as.vector(round(DR[10,1:5], 3))
  expect_equal(DRcheck, DRcheck0, tolerance = 1.5, scale = 1)
})

test_that("Chain creation works when StoreAdapt=TRUE (spikes)", 
{
  # Data example
  set.seed(8)
  Data <- makeExampleBASiCS_Data(WithSpikes = TRUE)
  # Fixing starting values
  n <- ncol(Data)
  PriorParam <- list(mu.mu = rep(0, times = nrow(Data)), s2.mu = 0.5, 
                     s2.delta = 0.5, a.delta = 1, 
                     b.delta = 1, p.phi = rep(1, times = n), 
                     a.s = 1, b.s = 1, a.theta = 1, b.theta = 1)
  set.seed(2018)
  Start <- BASiCS:::HiddenBASiCS_MCMC_Start(Data, PriorParam, WithSpikes = TRUE)
  # Running the samples
  set.seed(18)
  Chain <- run_MCMC(Data, N = 50, Thin = 10, Burn = 10,
                       Regression = FALSE, PrintProgress = FALSE, 
                       StoreAdapt = TRUE,
                       Start = Start, PriorParam = PriorParam)
  expect_s4_class(Chain, "BASiCS_Chain")
})
