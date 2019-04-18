context("plots")

test_that("plot of BASiCS_Summary works without spikes", {
  # Data example
  set.seed(19)
  Data <- makeExampleBASiCS_Data(WithSpikes = FALSE, WithBatch = TRUE)
  set.seed(42)
  n <- ncol(Data)
  k <- 12

  PriorParam <- list(s2.mu = 0.5, s2.delta = 0.5, a.delta = 1,
                   b.delta = 1, p.phi = rep(1, times = n),
                   a.s = 1, b.s = 1, a.theta = 1, b.theta = 1)

  Start <- BASiCS:::HiddenBASiCS_MCMC_Start(
    Data,
    PriorParam,
    WithSpikes = FALSE)

  Chain <- run_MCMC(
    Data = Data,
    N = 100,
    Thin = 10,
    Burn = 10,
    PrintProgress = FALSE,
    Regression = FALSE,
    WithSpikes = FALSE)
  S <- Summary(Chain)
  pdf(NULL)
  expect_error(plot(S), NA)
  dev.off()
})

test_that("plot works for BASiCS_Chain (non-regression, spikes)", {
  Data <- makeExampleBASiCS_Data(WithSpikes = TRUE)
  set.seed(42)
  Chain <- run_MCMC(
    Data = Data,
    N = 100,
    Thin = 10,
    Burn = 10,
    PrintProgress = FALSE,
    Regression = FALSE,
    WithSpikes = TRUE)
  pdf(NULL)
  plot(Chain, Param = "mu", Gene = 1)
  plot(Chain, Param = "delta", Gene = 1)
  plot(Chain, Param = "phi", Cell = 1)
  plot(Chain, Param = "s", Cell = 1)
  plot(Chain, Param = "nu", Cell = 1)
  plot(Chain, Param = "theta")
  dev.off()
})

test_that("plot works for BASiCS_Chain (regression, no spikes)", {

  Data <- makeExampleBASiCS_Data(WithSpikes = FALSE, WithBatch = TRUE)
  set.seed(42)
  Chain <- run_MCMC(
    Data = Data,
    N = 100,
    Thin = 10,
    Burn = 10,
    PrintProgress = FALSE,
    Regression = TRUE,
    WithSpikes = FALSE)

  pdf(NULL)
  expect_error({
    plot(Chain, Param = "mu", Gene = 1)
    plot(Chain, Param = "delta", Gene = 1)
    plot(Chain, Param = "epsilon", Gene = 1)
    plot(Chain, Param = "s", Cell = 1)
    plot(Chain, Param = "nu", Cell = 1)
    plot(Chain, Param = "theta")
    plot(Chain, Param = "beta", RegressionTerm = 1)
  }, NA)
  dev.off()
})

test_that("plot works for BASiCS_Summary with all valid combinations", {
  Data <- makeExampleBASiCS_Data(WithSpikes = TRUE)
  set.seed(42)
  Chain <- run_MCMC(
    Data = Data,
    N = 100,
    Thin = 10,
    Burn = 10,
    PrintProgress = FALSE,
    Regression = TRUE,
    WithSpikes = TRUE)
  
  SChain <- Summary(Chain)
  pdf(NULL)
  expect_error({
    plot(SChain, Param = "mu", Param2 = "delta")
    plot(SChain, Param = "mu", Param2 = "epsilon")
    plot(SChain, Param = "delta", Param2 = "epsilon")
    plot(SChain, Param = "s", Param2 = "phi")
    plot(SChain, Param = "s", Param2 = "nu")
    plot(SChain, Param = "phi", Param2 = "nu")
    plot(SChain, Param = "sigma2")
    plot(SChain, Param = "beta")
  }, NA)
  dev.off()
})


test_that("BASiCS_showFit works", {
  Data <- makeExampleBASiCS_Data(WithSpikes = TRUE)
  set.seed(42)
  Chain <- run_MCMC(
    Data = Data,
    N = 100,
    Thin = 10,
    Burn = 10,
    PrintProgress = FALSE,
    Regression = TRUE,
    WithSpikes = TRUE)

  pdf(NULL)
  expect_error({
    BASiCS_showFit(Chain, smooth = TRUE)
    BASiCS_showFit(Chain, smooth = FALSE)
  }, NA)
  dev.off()

})


test_that("Diagnostic plot works", {
  Data <- makeExampleBASiCS_Data(WithSpikes = TRUE)
  set.seed(42)
  Chain <- run_MCMC(
    Data = Data,
    N = 100,
    Thin = 10,
    Burn = 10,
    PrintProgress = FALSE,
    Regression = TRUE,
    WithSpikes = TRUE)
  g <- BASiCS_diagPlot(Chain)
  expect_is(g, "ggplot")
  g <- BASiCS_diagPlot(Chain, Param = "delta")
  expect_is(g, "ggplot")
  g <- BASiCS_diagPlot(Chain, Param = "epsilon")
  expect_is(g, "ggplot")
  g <- BASiCS_diagPlot(Chain, Param = "delta", x = "mu", y = "epsilon")
  expect_is(g, "ggplot")
  expect_error(
    BASiCS_diagPlot(Chain, Param = "nu", x = "mu", y = "epsilon"),
    "Invalid combination of parameters"
  )
  expect_error(
    BASiCS_diagPlot(Chain, Param = "delta", x = "delta", y = "nu"),
    "Invalid combination of parameters:"
  )
  g <- BASiCS_diagPlot(Chain, Param = "delta", x = "mu", y = "epsilon")
  expect_is(g, "ggplot")

  Chain@parameters[["epsilon"]][, 1] <- NA
  g <- BASiCS_diagPlot(Chain, Param = "epsilon", x = "mu", y = "epsilon")
  expect_is(g, "ggplot")
  g <- BASiCS_diagPlot(Chain, Param = "epsilon")
  expect_is(g, "ggplot")
})


test_that("Diagnostic hist work", {
  Data <- makeExampleBASiCS_Data(WithSpikes = TRUE)
  set.seed(42)
  Chain <- run_MCMC(
    Data,
    N = 100,
    Thin = 10,
    Burn = 10,
    PrintProgress = FALSE,
    Regression = TRUE,
    WithSpikes = TRUE)
  expect_warning(g <- BASiCS_diagHist(Chain),
    "coda::effectiveSize failed for some parameters: theta, sigma2"
  )
  expect_is(g, "ggplot")
  g <- BASiCS_diagHist(Chain, Param = "delta")
  expect_is(g, "ggplot")
  g <- BASiCS_diagHist(Chain, Param = "nu")
  expect_is(g, "ggplot")
  g <- BASiCS_diagHist(Chain, Param = "nu")
  expect_is(g, "ggplot")
})
