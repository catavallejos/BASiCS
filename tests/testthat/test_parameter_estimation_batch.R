context("Parameter estimation and denoised data (spikes+batch)\n")

test_that("Estimates match the given seed (spikes+batch)", {

  Data <- makeExampleBASiCS_Data(WithSpikes = TRUE, WithBatch = TRUE)
  # Fixing starting values
  n <- ncol(Data)
  PriorParam <- list(s2.mu = 0.5, s2.delta = 0.5, a.delta = 1, 
                     b.delta = 1, p.phi = rep(1, times = n), 
                     a.s = 1, b.s = 1, a.theta = 1, b.theta = 1)
  set.seed(2018)
  Start <- BASiCS:::HiddenBASiCS_MCMC_Start(Data, PriorParam, WithSpikes = TRUE)
  # Running the samples
  set.seed(18)
  Chain <- BASiCS_MCMC(Data, N = 1000, Thin = 10, Burn = 500, 
                       PrintProgress = FALSE, 
                       Start = Start, PriorParam = PriorParam)
  # Calculating a posterior summary
  PostSummary <- Summary(Chain)
            
  # Check if parameter estimates match for the first 5 genes and cells
  Mu <- c(6.875,  4.750,  4.015,  5.352, 18.519)
  MuObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "mu")[1:5,1],3))
  expect_that(all.equal(MuObs, Mu), is_true())
            
  Delta <- c(1.161, 1.814, 0.654, 1.588, 0.699)
  DeltaObs <- as.vector(round(displaySummaryBASiCS(PostSummary, 
                                                   "delta")[1:5,1],3))
  expect_that(all.equal(DeltaObs, Delta), is_true())
            
  Phi <- c(1.161, 1.017, 0.762, 1.012, 0.964)
  PhiObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "phi")[1:5,1],3))
  expect_that(all.equal(PhiObs, Phi), is_true())
            
  S <- c(0.304, 0.525, 0.123, 0.269, 0.560)
  SObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "s")[1:5,1],3))
  expect_that(all.equal(SObs, S), is_true())
            
  Theta <- c(0.705, 0.344)
  ThetaObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "theta")[,1],3))
  expect_that(all.equal(ThetaObs, Theta), is_true())
            
  # Obtaining denoised counts     
  DC <- BASiCS_DenoisedCounts(Data, Chain)
            
  # Checks for 2 arbitrary sets of genes / cells
  DCcheck0 <- c(0.000, 0.000, 0.000, 4.537, 4.537)
  DCcheck <- as.vector(round(DC[1:5,1], 3))
  expect_that(all.equal(DCcheck, DCcheck0), is_true())
            
  DCcheck0 <- c(0.000, 2.228, 0.000, 0.000, 2.478)
  DCcheck <- as.vector(round(DC[10,1:5], 3))
  expect_that(all.equal(DCcheck, DCcheck0), is_true())
            
  # Obtaining denoised rates
  DR <- BASiCS_DenoisedRates(Data, Chain)
            
  # Checks for 2 arbitrary sets of genes / cells
  DRcheck0 <- c(2.490, 1.640, 2.476, 4.905, 8.260)
  DRcheck <- as.vector(round(DR[1:5,1], 3))
  expect_that(all.equal(DRcheck, DRcheck0), is_true())
            
  DRcheck0 <- c(2.118, 3.136, 3.848, 2.609, 3.388)
  DRcheck <- as.vector(round(DR[10,1:5], 3))
  expect_that(all.equal(DRcheck, DRcheck0), is_true())
})