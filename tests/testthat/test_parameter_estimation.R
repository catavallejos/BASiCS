context("Parameter estimation and denoised data (spikes) \n")

test_that("Estimates match the given seed (spikes)",
{
  # Data example
  Data <- makeExampleBASiCS_Data(WithSpikes = TRUE)
  # Fixing starting values
  n <- ncol(Data)
  PriorParam <- list(s2.mu = 0.5, s2.delta = 0.5, a.delta = 1,
                     b.delta = 1, p.phi = rep(1, times = n),
                     a.s = 1, b.s = 1, a.theta = 1, b.theta = 1)
  set.seed(2018)
  Start <- BASiCS:::HiddenBASiCS_MCMC_Start(Data, PriorParam, WithSpikes = TRUE)
  # Running the samples
  set.seed(18)
  Chain <- BASiCS_MCMC(Data, N = 1000, Thin = 10, Burn = 500,
                       Regression = FALSE, PrintProgress = FALSE,
                       Start = Start, PriorParam = PriorParam)
  # Calculating a posterior summary
  PostSummary <- Summary(Chain)

  # Checking parameter names
  ParamNames <- c("mu", "delta", "phi", "s", "nu", "theta")
  expect_that(all.equal(names(Chain@parameters), ParamNames), is_true())
  expect_that(all.equal(names(PostSummary@parameters), ParamNames), is_true())

  # Check if parameter estimates match for the first 5 genes and cells
  Mu <- c(7.828,  7.290,  4.166,  5.286, 20.882)
  MuObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "mu")[1:5,1],3))
  expect_that(all.equal(MuObs, Mu, tolerance = 1, scale = 1), is_true())

  Delta <- c(1.183, 1.941, 0.691, 1.519, 0.656)
  DeltaObs <- as.vector(round(displaySummaryBASiCS(PostSummary,
                                                   "delta")[1:5,1],3))
  expect_that(all.equal(DeltaObs, Delta, tolerance = 1, scale = 1), is_true())

  Phi <- c(1.064, 0.986, 0.581, 0.961, 0.830)
  PhiObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "phi")[1:5,1],3))
  expect_that(all.equal(PhiObs, Phi, tolerance = 1, scale = 1), is_true())

  S <- c(0.390, 0.742, 0.110, 0.232, 0.593)
  SObs <- as.vector(round(displaySummaryBASiCS(PostSummary, "s")[1:5,1],3))
  expect_that(all.equal(SObs, S, tolerance = 1, scale = 1), is_true())

  Theta <- 0.541
  ThetaObs <- round(displaySummaryBASiCS(PostSummary, "theta")[1],3)
  expect_that(all.equal(ThetaObs, Theta, tolerance = 1, scale = 1), is_true())

  # Obtaining denoised counts
  DC <- BASiCS_DenoisedCounts(Data, Chain)
  # Checks for an arbitrary set of genes / cells
  DCcheck0 <- c(0.000, 0.000, 0.000, 4.935, 4.935)
  DCcheck <- as.vector(round(DC[1:5,1], 3))
  expect_that(all.equal(DCcheck, DCcheck0, tolerance = 1.5, scale = 1), is_true())

  # Obtaining denoised rates
  DR <- BASiCS_DenoisedRates(Data, Chain)
  # Checks for an arbitrary set of genes / cells
  DRcheck0 <- c(2.107, 2.918, 3.661, 2.517, 3.406)
  DRcheck <- as.vector(round(DR[10,1:5], 3))
  expect_that(all.equal(DRcheck, DRcheck0, tolerance = 1.5, scale = 1), is_true())
})

test_that("Chain creation works when StoreAdapt=TRUE (spikes)",
{
  # Data example
  Data <- makeExampleBASiCS_Data(WithSpikes = TRUE)
  # Fixing starting values
  n <- ncol(Data)
  PriorParam <- list(s2.mu = 0.5, s2.delta = 0.5, a.delta = 1,
                     b.delta = 1, p.phi = rep(1, times = n),
                     a.s = 1, b.s = 1, a.theta = 1, b.theta = 1)
  set.seed(2018)
  Start <- BASiCS:::HiddenBASiCS_MCMC_Start(Data, PriorParam, WithSpikes = TRUE)
  # Running the samples
  set.seed(18)
  Chain <- run_BASiCS_MCMC(Data, N = 50, Thin = 10, Burn = 10,
                       Regression = FALSE, PrintProgress = FALSE,
                       StoreAdapt = TRUE,
                       Start = Start, PriorParam = PriorParam)
  expect_s4_class(Chain, "BASiCS_Chain")
})
