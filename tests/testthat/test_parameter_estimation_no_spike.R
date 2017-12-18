context("Parameter estimation (no spike case), original data doesn't have spikes\n")

test_that("paramater estimations match the given seed", 
{
  Data <- makeExampleBASiCS_Data(WithSpikes = FALSE, 
                                 WithBatch = TRUE)
  set.seed(14)
  Chain <- BASiCS_MCMC(Data, N = 2000, Thin = 10, Burn = 1000, 
                       PrintProgress = FALSE, WithSpikes = FALSE)
  PostSummary <- Summary(Chain)
            
  # Check if parameter estimates match for the first 5 genes and cells
  Mu <- c(16.141, 11.314,  9.405, 11.561, 32.808)
  MuObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "mu")[1:5,1],3))
  expect_that(all.equal(MuObs, Mu), is_true())
            
  Delta <- c(1.183, 1.143, 1.817, 1.020, 0.479)
  DeltaObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "delta")[1:5,1],3))
  expect_that(all.equal(DeltaObs, Delta), is_true())

  Phi <- rep(1, 5)
  PhiObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "phi")[1:5,1],3))
  expect_that(all.equal(PhiObs, Phi), is_true())
            
  S <- c(0.686, 1.434, 0.237, 0.632, 1.391)
  SObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "s")[1:5,1],3))
  expect_that(all.equal(SObs, S), is_true())
  
  Theta <- c(0.283, 0.338)
  ThetaObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "theta")[,1],3))
  expect_that(all.equal(ThetaObs, Theta), is_true())
})

