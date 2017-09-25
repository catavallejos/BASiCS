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
  Mu <- c(7.071,  5.032,  4.379,  5.256, 18.358)
  MuObs <- as.vector(round(PostSummary@mu[1:5,1],3))
  expect_that(all.equal(MuObs, Mu), is_true())
            
  Delta <- c(1.226, 2.085, 0.597, 1.475, 0.598)
  DeltaObs <- as.vector(round(PostSummary@delta[1:5,1],3))
  expect_that(all.equal(DeltaObs, Delta), is_true())

  Phi <- c(1.084, 1.032, 0.930, 0.940, 0.917)
  PhiObs <- as.vector(round(PostSummary@phi[1:5,1],3))
  expect_that(all.equal(PhiObs, Phi), is_true())
            
  S <- c(0.304, 0.546, 0.092, 0.229, 0.527)
  SObs <- as.vector(round(PostSummary@s[1:5,1],3))
  expect_that(all.equal(SObs, S), is_true())
  
  Theta <- c(0.533, 0.354)
  ThetaObs <- as.vector(round(PostSummary@theta[,1],3))
  expect_that(all.equal(ThetaObs, Theta), is_true())
})

