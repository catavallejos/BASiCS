context("plots")

test_that("plot of BASiCS_Summary works without spikes", {
  # Data example
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

  capture.output(
    Chain <- BASiCS_MCMC(
      Data = Data, 
      N = 100, 
      Thin = 10, 
      Burn = 10, 
      PrintProgress = FALSE, 
      Regression = FALSE,
      WithSpikes = FALSE)
  )
  print(plot)
  S <- Summary(Chain)
  pdf(NULL)
  expect_error(plot(S), NA)
  dev.off()
})
