context("Parameter estimation\n")

test_that("paramater estimations match the given seed", 
{
  Data <- makeExampleBASiCS_Data(WithSpikes = TRUE)
  set.seed(18)
  Chain <- BASiCS_MCMC(Data, N = 10000, Thin = 10, Burn = 5000, 
                       PrintProgress = FALSE)
  PostSummary <- Summary(Chain)
            
  # Check if parameter estimates match for the first 5 genes and cells
  Mu <- c(7.403,  4.795,  4.098,  5.278, 18.812)
  MuObs <- as.vector(round(PostSummary@mu[1:5,1],3))
  expect_that(all.equal(MuObs, Mu), is_true())
            
  Delta <- c(1.169, 2.308, 0.620, 1.594, 0.577)
  DeltaObs <- as.vector(round(PostSummary@delta[1:5,1],3))
  expect_that(all.equal(DeltaObs, Delta), is_true())

  Phi <- c(1.084, 1.014, 0.918, 0.989, 0.908)
  PhiObs <- as.vector(round(PostSummary@phi[1:5,1],3))
  expect_that(all.equal(PhiObs, Phi), is_true())
            
  S <- c(0.282, 0.532, 0.081, 0.213, 0.528)
  SObs <- as.vector(round(PostSummary@s[1:5,1],3))
  expect_that(all.equal(SObs, S), is_true())
})

