context("Parameter estimation (regression case)\n")

test_that("paramater estimations match the given seed", 
{
  Data <- makeExampleBASiCS_Data(WithSpikes = TRUE)
  set.seed(12)
  Chain <- BASiCS_MCMC(Data, N = 2000, Thin = 10, Burn = 1000, 
                       PrintProgress = FALSE, Regression = TRUE)
  PostSummary <- Summary(Chain)
            
  # Check if parameter estimates match for the first 5 genes and cells
  Mu <- c( 7.874,  4.796,  4.091,  5.103, 18.225)
  MuObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "mu")[1:5,1],3))
  expect_that(all.equal(MuObs, Mu), is_true())
            
  Delta <- c(1.129, 2.111, 1.196, 1.520, 0.533)
  DeltaObs <- as.vector(round(displaySummaryBASiCS(PostSummary, 
                                                   "delta")[1:5,1],3))
  expect_that(all.equal(DeltaObs, Delta), is_true())
            
  Phi <- c(1.258, 1.093, 0.853, 1.037, 0.920)
  PhiObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "phi")[1:5,1],3))
  expect_that(all.equal(PhiObs, Phi), is_true())
            
  S <- c(0.254, 0.545, 0.080, 0.219, 0.558)
  SObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "s")[1:5,1],3))
  expect_that(all.equal(SObs, S), is_true())
            
  Theta <- 0.404
  ThetaObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "theta")[,1],3))
  expect_that(all.equal(ThetaObs, Theta), is_true())
})

