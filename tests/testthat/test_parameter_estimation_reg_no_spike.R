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
  Mu <- c(15.870, 10.173,  8.935, 11.256, 31.471)
  MuObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "mu")[1:5,1],3))
  expect_that(all.equal(MuObs, Mu), is_true())
            
  Delta <- c(1.309, 1.384, 1.983, 1.131, 0.412)
  DeltaObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "delta")[1:5,1],3))
  expect_that(all.equal(DeltaObs, Delta), is_true())

  Phi <- rep(1, 5)
  PhiObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "phi")[1:5,1],3))
  expect_that(all.equal(PhiObs, Phi), is_true())
            
  S <- c(0.767, 1.463, 0.227, 0.607, 1.431)
  SObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "s")[1:5,1],3))
  expect_that(all.equal(SObs, S), is_true())
  
  Theta <- c(0.225, 0.236)
  ThetaObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "theta")[,1],3))
  expect_that(all.equal(ThetaObs, Theta), is_true())
})

