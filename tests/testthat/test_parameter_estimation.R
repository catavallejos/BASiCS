context("Parameter estimation\n")

test_that("paramater estimations match the given seed", 
{
  Data <- makeExampleBASiCS_Data(WithSpikes = TRUE)
  set.seed(18)
  Chain <- BASiCS_MCMC(Data, N = 10000, Thin = 10, Burn = 5000, 
                       PrintProgress = FALSE)
  PostSummary <- Summary(Chain)
            
  # Check if parameter estimates match for the first 5 genes and cells
  Mu <- c(7.627, 4.786,  4.024,  5.290, 18.464)
  MuObs <- as.vector(round(PostSummary@parameters$mu[1:5,1],3))
  expect_that(all.equal(MuObs, Mu), is_true())
            
  Delta <- c(1.161, 2.118, 0.721, 1.496, 0.646)
  DeltaObs <- as.vector(round(PostSummary@parameters$delta[1:5,1],3))
  expect_that(all.equal(DeltaObs, Delta), is_true())

  Phi <- c(1.103, 1.027, 0.911, 1.018, 0.892)
  PhiObs <- as.vector(round(PostSummary@parameters$phi[1:5,1],3))
  expect_that(all.equal(PhiObs, Phi), is_true())
            
  S <- c(0.272, 0.526, 0.073, 0.221, 0.502)
  SObs <- as.vector(round(PostSummary@parameters$s[1:5,1],3))
  expect_that(all.equal(SObs, S), is_true())
})

