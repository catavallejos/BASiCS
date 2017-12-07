context("Parameter estimation (no spike case), original data doesn't have spikes\n")

test_that("paramater estimations match the given seed", 
{
  Data <- makeExampleBASiCS_Data(WithSpikes = FALSE, 
                                 WithBatch = TRUE)
  set.seed(14)
  Chain <- BASiCS_MCMC(Data, N = 10000, Thin = 10, Burn = 5000, 
                       PrintProgress = FALSE, WithSpikes = FALSE)
  PostSummary <- Summary(Chain)
            
  # Check if parameter estimates match for the first 5 genes and cells
  Mu <- c(15.775, 10.513,  8.975, 11.168, 32.375)
  MuObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "mu")[1:5,1],3))
  expect_that(all.equal(MuObs, Mu), is_true())
            
  Delta <- c( 1.253, 1.186, 1.741, 1.059, 0.476)
  DeltaObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "delta")[1:5,1],3))
  expect_that(all.equal(DeltaObs, Delta), is_true())

  Phi <- rep(1, 5)
  PhiObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "phi")[1:5,1],3))
  expect_that(all.equal(PhiObs, Phi), is_true())
            
  S <- c(0.738, 1.403, 0.226, 0.639, 1.430)
  SObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "s")[1:5,1],3))
  expect_that(all.equal(SObs, S), is_true())
  
  Theta <- c(0.259, 0.259)
  ThetaObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "theta")[,1],3))
  expect_that(all.equal(ThetaObs, Theta), is_true())
})

