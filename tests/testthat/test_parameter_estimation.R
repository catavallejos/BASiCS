context("Parameter estimation\n")

test_that("paramater estimations match the given seed", 
{
  Data <- makeExampleBASiCS_Data(WithSpikes = TRUE)
  set.seed(18)
  Chain <- BASiCS_MCMC(Data, N = 10000, Thin = 10, Burn = 5000, 
                       PrintProgress = FALSE)
  PostSummary <- Summary(Chain)
            
  # Check if parameter estimates match for the first 5 genes and cells
  Mu <- c(7.827, 4.818,  4.079,  5.101, 18.603)
  MuObs <- as.vector(round(PostSummary@parameters$mu[1:5,1],3))
  expect_that(all.equal(MuObs, Mu), is_true())
            
  Delta <- c(1.140, 2.030, 0.727, 1.562, 0.646)
  DeltaObs <- as.vector(round(PostSummary@parameters$delta[1:5,1],3))
  expect_that(all.equal(DeltaObs, Delta), is_true())

  Phi <- c(1.102, 1.034, 0.872, 1.005, 0.906)
  PhiObs <- as.vector(round(PostSummary@parameters$phi[1:5,1],3))
  expect_that(all.equal(PhiObs, Phi), is_true())
            
  S <- c(0.290, 0.548, 0.082, 0.225, 0.505)
  SObs <- as.vector(round(PostSummary@parameters$s[1:5,1],3))
  expect_that(all.equal(SObs, S), is_true())
})

