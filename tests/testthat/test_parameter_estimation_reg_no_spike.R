context("Parameter estimation no spikes BASiCS sampler (regression case)\n")

test_that("paramater estimations match the given seed", 
{
  Data <- makeExampleBASiCS_Data(WithSpikes = FALSE, WithBatch = TRUE)
  set.seed(14)
  Chain <- BASiCS_MCMC(Data, N = 10000, Thin = 10, Burn = 5000, 
                       PrintProgress = FALSE, WithSpikes = FALSE,
                       Regression = TRUE, k = 12, Var = 1.2)
  PostSummary <- Summary(Chain)
            
  # Check if parameter estimates match for the first 5 genes and cells
  Mu <- c(17.266, 10.645,  9.045, 11.344, 32.515)
  MuObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "mu")[1:5,1],3))
  expect_that(all.equal(MuObs, Mu), is_true())
            
  Delta <- c(1.191, 1.352, 1.999, 1.083, 0.442)
  DeltaObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "delta")[1:5,1],3))
  expect_that(all.equal(DeltaObs, Delta), is_true())

  Phi <- rep(1, 5)
  PhiObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "phi")[1:5,1],3))
  expect_that(all.equal(PhiObs, Phi), is_true())
            
  S <- c(0.782, 1.357, 0.225, 0.616, 1.397)
  SObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "s")[1:5,1],3))
  expect_that(all.equal(SObs, S), is_true())
  
  Theta <- c(0.256, 0.226)
  ThetaObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "theta")[,1],3))
  expect_that(all.equal(ThetaObs, Theta), is_true())
})

