context("Parameter estimation and denoised data (no-spikes+regression)\n")

test_that("Estimates match the given seed (no-spikes+regression)", 
{
  # Data example
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
  set.seed(14)
  Chain <- BASiCS_MCMC(Data, N = 1000, Thin = 10, Burn = 500,
                       PrintProgress = FALSE, WithSpikes = FALSE,
                       Regression = TRUE,
                       Start = Start, PriorParam = PriorParam)

  # Calculating a posterior summary
  PostSummary <- Summary(Chain)
  
  # Checking parameter names
  ParamNames <- c("mu", "delta", "s", "nu", "theta",
                  "beta", "sigma2", "epsilon", "RefFreq")
  ParamNames1 <- c("mu", "delta", "s", "nu", "theta",
                  "beta", "sigma2", "epsilon")
  expect_that(all.equal(names(Chain@parameters), ParamNames), is_true())
  expect_that(all.equal(names(PostSummary@parameters), ParamNames1), is_true())
            
  # Check if parameter estimates match for the first 5 genes and cells
  Mu <- c(16.482, 11.091,  8.765, 12.476, 33.748)
  MuObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "mu")[1:5,1],3))
  expect_that(all.equal(MuObs, Mu), is_true())
            
  Delta <- c(1.294, 1.175, 1.783, 1.209, 0.442)
  DeltaObs <- as.vector(round(displaySummaryBASiCS(PostSummary, 
                                                   "delta")[1:5,1],3))
  expect_that(all.equal(DeltaObs, Delta), is_true())
            
  S <- c(0.657, 1.407, 0.208, 0.614, 1.348)
  SObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "s")[1:5,1],3))
  expect_that(all.equal(SObs, S), is_true())
  
  Theta <- c(0.072, 0.235)
  ThetaObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "theta")[,1],3))
  expect_that(all.equal(ThetaObs, Theta), is_true())
  
  Beta <- c(0.151, -0.264,  0.225,  0.293,  0.476)
  BetaObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "beta")[1:5,1],3))
  expect_that(all.equal(BetaObs, Beta), is_true())
  
  Sigma2 <- 0.358
  Sigma2Obs <- round(displaySummaryBASiCS(PostSummary, "sigma2")[1],3)
  expect_that(all.equal(Sigma2Obs, Sigma2), is_true())
  
  # Obtaining denoised counts     
  DC <- BASiCS_DenoisedCounts(Data, Chain)
  # Checks for an arbitrary set of genes / cells
  DCcheck0 <- c(4.535,  1.512,  0.000,  9.070, 10.582)
  DCcheck <- as.vector(round(DC[1:5,1], 3))
  expect_that(all.equal(DCcheck, DCcheck0), is_true())
  
  # Obtaining denoised rates
  DR <- BASiCS_DenoisedRates(Data, Chain)
  # Checks for an arbitrary set of genes / cells
  DRcheck0 <- c(5.418,  3.275, 11.200,  1.642,  5.057)
  DRcheck <- as.vector(round(DR[10,1:5], 3))
  expect_that(all.equal(DRcheck, DRcheck0), is_true())
})

test_that("Chain creation works when regression, no spikes, and StoreAdapt=TRUE", {
  # Data example
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
  Chain <- BASiCS_MCMC(Data, N = 100, Thin = 5, Burn = 5,
                     PrintProgress = FALSE, WithSpikes = FALSE,
                     Regression = TRUE, StoreAdapt = TRUE,
                     Start = Start, PriorParam = PriorParam)
})
