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
  Mu <- c(12.766, 10.717,  7.324, 10.877, 30.084)
  MuObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "mu")[1:5,1],3))
  expect_that(all.equal(MuObs, Mu), is_true())
            
  Delta <- c(1.002, 0.993, 1.875, 1.355, 0.605)
  DeltaObs <- as.vector(round(displaySummaryBASiCS(PostSummary, 
                                                   "delta")[1:5,1],3))
  expect_that(all.equal(DeltaObs, Delta), is_true())
            
  S <- c(0.859, 1.452, 0.225, 0.611, 1.362)
  SObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "s")[1:5,1],3))
  expect_that(all.equal(SObs, S), is_true())
  
  Theta <- c(0.172, 0.128)
  ThetaObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "theta")[,1],3))
  expect_that(all.equal(ThetaObs, Theta), is_true())
  
  Beta <- c(0.385, -0.338,  0.536,  0.294,  0.211)
  BetaObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "beta")[1:5,1],3))
  expect_that(all.equal(BetaObs, Beta), is_true())
  
  Sigma2 <- 0.35
  Sigma2Obs <- round(displaySummaryBASiCS(PostSummary, "sigma2")[1],3)
  expect_that(all.equal(Sigma2Obs, Sigma2), is_true())
  
  # Obtaining denoised counts     
  DC <- BASiCS_DenoisedCounts(Data, Chain)
  # Checks for an arbitrary set of genes / cells
  DCcheck0 <- c(3.533,  1.178,  0.000,  7.066, 8.244)
  DCcheck <- as.vector(round(DC[1:5,1], 3))
  expect_that(all.equal(DCcheck, DCcheck0), is_true())
  
  # Obtaining denoised rates
  DR <- BASiCS_DenoisedRates(Data, Chain)
  # Checks for an arbitrary set of genes / cells
  DRcheck0 <- c(4.234,  2.923, 9.706,  1.517,  4.359)
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
  Chain <- BASiCS_MCMC(Data, N = 50, Thin = 10, Burn = 10,
                     PrintProgress = FALSE, WithSpikes = FALSE,
                     Regression = TRUE, StoreAdapt = TRUE,
                     Start = Start, PriorParam = PriorParam)
  expect_s4_class(Chain, "BASiCS_Chain")
})
