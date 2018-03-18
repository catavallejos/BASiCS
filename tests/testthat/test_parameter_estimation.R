context("Parameter estimation and denoised data (spikes) \n")

test_that("Estimates match the given seed (spikes)", {
  # Data example
  Data <- makeExampleBASiCS_Data(WithSpikes = TRUE)
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
  Mu <- c(7.271,  4.842,  4.322,  4.865, 18.003)
  MuObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "mu")[1:5,1],3))
  expect_that(all.equal(MuObs, Mu), is_true())
            
  Delta <- c(1.151, 2.024, 0.681, 1.441, 0.607)
  DeltaObs <- as.vector(round(displaySummaryBASiCS(PostSummary, 
                                                   "delta")[1:5,1],3))
  expect_that(all.equal(DeltaObs, Delta), is_true())
            
  Phi <- c(0.976, 0.986, 1.003, 1.012, 0.880)
  PhiObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "phi")[1:5,1],3))
  expect_that(all.equal(PhiObs, Phi), is_true())
            
  S <- c(0.360, 0.692, 0.122, 0.225, 0.593)
  SObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "s")[1:5,1],3))
  expect_that(all.equal(SObs, S), is_true())
            
  # Obtaining denoised counts     
  DC <- BASiCS_DenoisedCounts(Data, Chain)
            
  # Checks for 2 arbitrary sets of genes / cells
  DCcheck0 <- c(0.000, 0.000, 0.000, 5.363, 5.363)
  DCcheck <- as.vector(round(DC[1:5,1], 3))
  expect_that(all.equal(DCcheck, DCcheck0), is_true())
            
  DCcheck0 <- c(0.000, 2.305, 0.000, 0.000, 2.711)
  DCcheck <- as.vector(round(DC[10,1:5], 3))
  expect_that(all.equal(DCcheck, DCcheck0), is_true())
            
  # Obtaining denoised rates
  DR <- BASiCS_DenoisedRates(Data, Chain)
            
  # Checks for 2 arbitrary sets of genes / cells
  DRcheck0 <- c(2.870, 1.781, 2.711, 5.060, 9.569)
  DRcheck <- as.vector(round(DR[1:5,1], 3))
  expect_that(all.equal(DRcheck, DRcheck0), is_true())
            
  DRcheck0 <- c(2.393, 3.195, 3.713, 2.698, 3.577)
  DRcheck <- as.vector(round(DR[10,1:5], 3))
  expect_that(all.equal(DRcheck, DRcheck0), is_true())
})