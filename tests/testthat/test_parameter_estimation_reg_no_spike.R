context("Parameter estimation no spikes BASiCS sampler (regression case)\n")

test_that("paramater estimations match the given seed", 
{
  Data <- makeExampleBASiCS_Data(WithSpikes = FALSE, WithBatch = TRUE)
  set.seed(14)
  Chain <- BASiCS_MCMC(Data, N = 2000, Thin = 10, Burn = 1000, 
                       PrintProgress = FALSE, WithSpikes = FALSE,
                       Regression = TRUE)
  PostSummary <- Summary(Chain)
            
  # Check if parameter estimates match for the first 5 genes and cells
  Mu <- c(16.956, 10.735,  8.572, 11.141, 33.688)
  MuObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "mu")[1:5,1],3))
  expect_that(all.equal(MuObs, Mu), is_true())
            
  Delta <- c(1.307, 1.494, 1.709, 1.007, 0.384)
  DeltaObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "delta")[1:5,1],3))
  expect_that(all.equal(DeltaObs, Delta), is_true())

  Phi <- rep(1, 5)
  PhiObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "phi")[1:5,1],3))
  expect_that(all.equal(PhiObs, Phi), is_true())
            
  S <- c(0.734, 1.353, 0.205, 0.558, 1.500)
  SObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "s")[1:5,1],3))
  expect_that(all.equal(SObs, S), is_true())
  
  Theta <- c(0.208, 0.218)
  ThetaObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "theta")[,1],3))
  expect_that(all.equal(ThetaObs, Theta), is_true())
})

