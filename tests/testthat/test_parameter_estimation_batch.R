context("Parameter estimation and denoised data (spikes+batch)")

test_that("Estimates match the given seed (spikes+batch)", {
  set.seed(9)
  Data <- makeExampleBASiCS_Data(
    WithSpikes = TRUE,
    WithBatch = TRUE
  )
  # Fixing starting values
  PriorParam <- BASiCS_PriorParam(Data, k = 12)
  set.seed(2018)
  Start <- BASiCS:::.BASiCS_MCMC_Start(
    Data,
    PriorParam,
    WithSpikes = TRUE,
    Regression = FALSE
  )
  # Running the samples
  set.seed(18)
  Chain <- run_MCMC(
    Data,
    N = 1000,
    Thin = 10,
    Burn = 500,
    Regression = FALSE,
    PrintProgress = FALSE,
    Start = Start,
    PriorParam = PriorParam
  )

  # Calculating a posterior summary
  PostSummary <- Summary(Chain)

  # Checking parameter names
  ParamNames <- c("mu", "delta", "phi", "s", "nu", "theta")
  expect_equal(names(Chain@parameters), ParamNames)
  expect_equal(names(PostSummary@parameters), ParamNames)

  # Check if parameter estimates match for the first 5 genes and cells
  Mu <- c(6.296, 9.44, 2.732, 5.504, 18.885)
  MuObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "mu")[1:5, 1],3))
  expect_equal(MuObs, Mu)

  Delta <- c(0.993, 0.766, 1.134, 1.413, 0.515)
  DeltaObs <- as.vector(
    round(displaySummaryBASiCS(PostSummary, "delta")[1:5, 1], 3)
  )
  expect_equal(DeltaObs, Delta)

  Phi <- c(1.118, 1, 0.913, 1.473, 0.62)
  PhiObs <- as.vector(
    round(displaySummaryBASiCS(PostSummary, "phi")[1:5, 1], 3)
  )
  expect_equal(PhiObs, Phi)

  S <- c(0.171, 0.215, 0.289, 1.03, 0.448)
  SObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "s")[1:5, 1], 3))
  expect_equal(SObs, S)

  Theta <- c(0.385, 0.808)
  ThetaObs <- as.vector(
    round(displaySummaryBASiCS(PostSummary, "theta")[,1], 3)
  )
  expect_equal(ThetaObs, Theta)

  # Obtaining denoised counts
  DC <- BASiCS_DenoisedCounts(Data, Chain)

  # Checks for an arbitrary set of genes / cells
  DCcheck0 <- c(0, 6.675, 20.025, 20.025, 40.049)
  DCcheck <- as.vector(round(DC[1:5, 1], 3))
  expect_equal(DCcheck, DCcheck0)

  # Obtaining denoised rates
  DR <- BASiCS_DenoisedRates(Data, Chain)

  # Checks for an arbitrary set of genes / cells
  DRcheck0 <- c(8.62, 15.017, 2.102, 11.07, 8.518)
  DRcheck <- as.vector(round(DR[10, 1:5], 3))
  expect_equal(DRcheck, DRcheck0)
})
