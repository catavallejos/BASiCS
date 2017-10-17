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
  expect_that(all.equal(ThetaObs, Theta), is_true())
})

