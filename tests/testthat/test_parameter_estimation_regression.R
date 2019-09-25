context("Parameter estimation and denoised data (spikes+regression)\n")

test_that("Estimates match the given seed (spikes+regression)", 
{
  # Data example
  set.seed(15)
  Data <- makeExampleBASiCS_Data(WithSpikes = TRUE, WithBatch = TRUE)
  # Fixing starting values
  n <- ncol(Data); k <- 12
  PriorParam <- list(s2.mu = 0.5, s2.delta = 0.5, a.delta = 1, 
                     b.delta = 1, p.phi = rep(1, times = n), 
                     a.s = 1, b.s = 1, a.theta = 1, b.theta = 1)
  PriorParam$m <- rep(0, k); PriorParam$V <- diag(k) 
  PriorParam$a.sigma2 <- 2; PriorParam$b.sigma2 <- 2  
  PriorParam$eta <- 5
  set.seed(2018)
  Start <- BASiCS:::HiddenBASiCS_MCMC_Start(Data, PriorParam, WithSpikes = TRUE)
  # Running the sampler
  set.seed(12)
  Chain <- run_MCMC(Data, N = 1000, Thin = 10, Burn = 500, 
                       PrintProgress = FALSE, Regression = TRUE,
                       Start = Start, PriorParam = PriorParam)
  # Calculating a posterior summary
  PostSummary <- Summary(Chain)
  
  # Checking parameter names
  ParamNames <- c("mu", "delta", "phi", "s", "nu", "theta",
                  "beta", "sigma2", "epsilon")
  expect_true(all.equal(names(Chain@parameters), ParamNames))
  expect_true(all.equal(names(PostSummary@parameters), ParamNames))
            
  # Check if parameter estimates match for the first 5 genes and cells
  Mu <- c(6.661, 11.292,  4.271,  3.818, 25.961)
  MuObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "mu")[1:5,1],3))
  expect_true(all.equal(MuObs, Mu))
            
  Delta <- c(1.343, 0.497, 1.932, 1.755, 0.426)
  DeltaObs <- as.vector(round(displaySummaryBASiCS(PostSummary, 
                                                   "delta")[1:5,1],3))
  expect_true(all.equal(DeltaObs, Delta))
            
  Phi <- c(0.902, 1.266, 0.891, 1.056, 0.780)
  PhiObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "phi")[1:5,1],3))
  expect_true(all.equal(PhiObs, Phi))
            
  S <- c(0.481, 1.031, 0.356, 0.242, 0.145)
  SObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "s")[1:5,1],3))
  expect_true(all.equal(SObs, S))
            
  Theta <- c(0.629, 0.366)
  ThetaObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "theta")[,1],3))
  expect_true(all.equal(ThetaObs, Theta))
  
  Beta <- c(0.429, -0.355,  0.485,  0.448,  0.230)
  BetaObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "beta")[1:5,1],3))
  expect_true(all.equal(BetaObs, Beta))
  
  Sigma2 <- 0.291
  Sigma2Obs <- round(displaySummaryBASiCS(PostSummary, "sigma2")[1],3)
  expect_true(all.equal(Sigma2Obs, Sigma2))
  
  # Obtaining denoised counts     
  DC <- BASiCS_DenoisedCounts(Data, Chain)
  
  # Checks for an arbitrary set of genes / cells
  DCcheck0 <- c(0.000,  9.489,  0.000, 22.140,  3.163)
  DCcheck <- as.vector(round(DC[1:5,1], 3))
  expect_true(all.equal(DCcheck, DCcheck0))
  
  # Obtaining denoised rates
  DR <- BASiCS_DenoisedRates(Data, Chain)
  
  # Checks for an arbitrary set of genes / cells
  DRcheck0 <- c( 28.273,  0.541,  2.229,  2.438,  4.218)
  DRcheck <- as.vector(round(DR[10,1:5], 3))
  expect_true(all.equal(DRcheck, DRcheck0))
})
test_that("Chain creation works when StoreAdapt=TRUE (spikes+regression)", 
{
  # Data example
  set.seed(18)
  Data <- makeExampleBASiCS_Data(WithSpikes = TRUE, WithBatch = TRUE)
  # Fixing starting values
  n <- ncol(Data); k <- 12
  PriorParam <- list(s2.mu = 0.5, s2.delta = 0.5, a.delta = 1, 
                     b.delta = 1, p.phi = rep(1, times = n), 
                     a.s = 1, b.s = 1, a.theta = 1, b.theta = 1)
  PriorParam$m <- rep(0, k); PriorParam$V <- diag(k) 
  PriorParam$a.sigma2 <- 2; PriorParam$b.sigma2 <- 2  
  PriorParam$eta <- 5
  set.seed(2018)
  Start <- BASiCS:::HiddenBASiCS_MCMC_Start(Data, PriorParam, WithSpikes = TRUE)
  # Running the sampler
  set.seed(12)
  Chain <- run_MCMC(Data, 
                    N = 50, Thin = 10, Burn = 10,
                    PrintProgress = FALSE, Regression = TRUE,
                    StoreAdapt = TRUE,Start = Start, PriorParam = PriorParam)
  expect_s4_class(Chain, "BASiCS_Chain")
})
