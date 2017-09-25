context("Parameter estimation (batch case)")

test_that("paramater estimations match the given seed", 
{
  Data <- makeExampleBASiCS_Data(WithSpikes = TRUE, Example = 1,
                                 WithBatch = TRUE)
  set.seed(18)
  Chain <- BASiCS_MCMC(Data, N = 10000, Thin = 10, Burn = 5000, 
                       PrintProgress = FALSE)
  PostSummary <- Summary(Chain)
            
  # Check if parameter estimates match for the first 5 genes and cells
  Mu <- c(7.132,  4.956,  4.207,  5.388, 18.447)
  MuObs <- as.vector(round(PostSummary@mu[1:5,1],3))
  expect_that(all.equal(MuObs, Mu), is_true())
            
  Delta <- c(1.217, 2.109, 0.625, 1.417, 0.595)
  DeltaObs <- as.vector(round(PostSummary@delta[1:5,1],3))
  expect_that(all.equal(DeltaObs, Delta), is_true())

  Phi <- c(1.097, 1.016, 0.855, 0.985, 0.939)
  PhiObs <- as.vector(round(PostSummary@phi[1:5,1],3))
  expect_that(all.equal(PhiObs, Phi), is_true())
            
  S <- c(0.319, 0.545, 0.097, 0.248, 0.538)
  SObs <- as.vector(round(PostSummary@s[1:5,1],3))
  expect_that(all.equal(SObs, S), is_true())
  
  Theta <- c(0.594, 0.435)
  ThetaObs <- as.vector(round(PostSummary@theta[,1],3))
  expect_that(all.equal(ThetaObs, Theta), is_true())
})

