context("Parameter estimation (batch case)\n")

test_that("paramater estimations match the given seed", 
{
  Data <- makeExampleBASiCS_Data(WithSpikes = TRUE, 
                                 WithBatch = TRUE)
  set.seed(18)
  Chain <- BASiCS_MCMC(Data, N = 10000, Thin = 10, Burn = 5000, 
                       PrintProgress = FALSE)
  PostSummary <- Summary(Chain)
            
  # Check if parameter estimates match for the first 5 genes and cells
<<<<<<< HEAD
  Mu <- c(7.002,  5.029,  4.272,  5.444, 18.058)
  MuObs <- as.vector(round(PostSummary@parameters$mu[1:5,1],3))
  expect_that(all.equal(MuObs, Mu), is_true())
            
  Delta <- c(1.217, 2.027, 0.735, 1.408, 0.633)
  DeltaObs <- as.vector(round(PostSummary@parameters$delta[1:5,1],3))
  expect_that(all.equal(DeltaObs, Delta), is_true())

  Phi <- c(1.117, 1.037, 0.938, 0.992, 0.925)
  PhiObs <- as.vector(round(PostSummary@parameters$phi[1:5,1],3))
  expect_that(all.equal(PhiObs, Phi), is_true())
            
  S <- c(0.339, 0.578, 0.101, 0.249, 0.554)
  SObs <- as.vector(round(PostSummary@parameters$s[1:5,1],3))
  expect_that(all.equal(SObs, S), is_true())
  
  Theta <- c(0.567, 0.489)
  ThetaObs <- as.vector(round(PostSummary@parameters$theta[,1],3))
=======
  Mu <- c(6.942,  5.007,  4.252,  5.378, 17.889)
  MuObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "mu")[1:5,1],3))
  expect_that(all.equal(MuObs, Mu), is_true())
            
  Delta <- c(1.193, 1.965, 0.721, 1.427, 0.677)
  DeltaObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "delta")[1:5,1],3))
  expect_that(all.equal(DeltaObs, Delta), is_true())

  Phi <- c(1.099, 1.042, 0.913, 0.952, 0.932)
  PhiObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "phi")[1:5,1],3))
  expect_that(all.equal(PhiObs, Phi), is_true())
            
  S <- c(0.284, 0.553, 0.093, 0.241, 0.512)
  SObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "s")[1:5,1],3))
  expect_that(all.equal(SObs, S), is_true())
  
  Theta <- c(0.501, 0.434)
  ThetaObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "theta")[,1],3))
>>>>>>> master
  expect_that(all.equal(ThetaObs, Theta), is_true())
})

