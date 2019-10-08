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

  data(ChainSC)
  S <- Summary(ChainSC)
  pdf(NULL)
  expect_error(plot(S), NA)
  dev.off()
})

test_that("plot works for BASiCS_Chain (non-regression, spikes)", {
  data(ChainSC)
  pdf(NULL)
  expect_error({
    plot(ChainSC, Param = "mu", Gene = 1)
    plot(ChainSC, Param = "delta", Gene = 1)
    plot(ChainSC, Param = "phi", Cell = 1)
    plot(ChainSC, Param = "s", Cell = 1)
    plot(ChainSC, Param = "nu", Cell = 1)
    plot(ChainSC, Param = "theta")
  }, NA)
  dev.off()
})

test_that("plot works for BASiCS_Chain (regression, no spikes)", {
  data(ChainSCReg)
  pdf(NULL)
  expect_error({
    plot(ChainSCReg, Param = "mu", Gene = 1)
    plot(ChainSCReg, Param = "delta", Gene = 1)
    plot(ChainSCReg, Param = "epsilon", Gene = 1)
    plot(ChainSCReg, Param = "s", Cell = 1)
    plot(ChainSCReg, Param = "nu", Cell = 1)
    plot(ChainSCReg, Param = "theta")
    plot(ChainSCReg, Param = "beta", RegressionTerm = 1)
  }, NA)
  dev.off()
})

test_that("plot works for BASiCS_Summary with all valid combinations", {
  data(ChainSCReg)
  SChain <- Summary(ChainSCReg)
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


test_that("BASiCS_ShowFit works", {
  data(ChainSCReg)
  expect_warning(BASiCS_showFit(ChainSCReg, smooth = FALSE), "deprecated")
  expect_error({
    BASiCS_ShowFit(ChainSCReg, smooth = TRUE)
    BASiCS_ShowFit(ChainSCReg, smooth = FALSE)
  }, NA)
})


test_that("Diagnostic plot works", {
  data(ChainSCReg)
  expect_warning(BASiCS_diagPlot(ChainSCReg), "deprecated")
  g <- BASiCS_DiagPlot(ChainSCReg)
  expect_is(g, "ggplot")
  g <- BASiCS_DiagPlot(ChainSCReg, Param = "delta")
  expect_is(g, "ggplot")
  g <- BASiCS_DiagPlot(ChainSCReg, Param = "epsilon")
  expect_is(g, "ggplot")
  g <- BASiCS_DiagPlot(ChainSCReg, Param = "delta", x = "mu", y = "epsilon")
  expect_is(g, "ggplot")
  expect_error(
    BASiCS_DiagPlot(ChainSCReg, Param = "nu", x = "mu", y = "epsilon"),
    "Invalid combination of parameters"
  )
  expect_error(
    BASiCS_DiagPlot(ChainSCReg, Param = "delta", x = "delta", y = "nu"),
    "Invalid combination of parameters:"
  )
  g <- BASiCS_DiagPlot(ChainSCReg, Param = "delta", x = "mu", y = "epsilon")
  expect_is(g, "ggplot")

  ChainSCReg@parameters[["epsilon"]][, 1] <- NA
  g <- BASiCS_DiagPlot(ChainSCReg, Param = "epsilon", x = "mu", y = "epsilon")
  expect_is(g, "ggplot")
  g <- BASiCS_DiagPlot(ChainSCReg, Param = "epsilon")
  expect_is(g, "ggplot")
})


test_that("Diagnostic hist work", {
  data(ChainSC)
  expect_warning(g <- BASiCS_diagHist(ChainSCReg), "deprecated")
  expect_is(g, "ggplot")
  g <- BASiCS_DiagHist(ChainSC, Param = "delta")
  expect_is(g, "ggplot")
  g <- BASiCS_DiagHist(ChainSC, Param = "nu")
  expect_is(g, "ggplot")
  g <- BASiCS_DiagHist(ChainSC, Param = "nu")
  expect_is(g, "ggplot")
})


test_that("DE plots work (non-regression)", {
  data(ChainSC)
  data(ChainRNA)
  Test <- BASiCS_TestDE(ChainSC, ChainRNA, Plot = FALSE, PlotOffset = FALSE)
  BASiCS_PlotDE(Test)
  BASiCS_PlotDE(Test@Results[[1]])
  BASiCS_PlotDE(Test@Results[[2]])

  for (type in c("MAPlot", "GridPlot", "VolcanoPlot")) {
    BASiCS_PlotDE(Test@Results[[1]], Which = type)
    BASiCS_PlotDE(Test@Results[[2]], Which = type)
  }
})

test_that("DE plots work (regression)", {
  data(ChainSCReg)
  data(ChainRNAReg)
  Test <- BASiCS_TestDE(ChainSCReg, ChainRNAReg, Plot = FALSE, PlotOffset = FALSE)
  BASiCS_PlotDE(Test)
  BASiCS_PlotDE(Test@Results[[1]])
  BASiCS_PlotDE(Test@Results[[2]])
  BASiCS_PlotDE(Test@Results[[3]])
  for (type in c("MAPlot", "GridPlot", "VolcanoPlot")) {
    BASiCS_PlotDE(Test@Results[[1]], Which = type)
    BASiCS_PlotDE(Test@Results[[2]], Which = type)
    BASiCS_PlotDE(Test@Results[[3]], Which = type)
  }
})
