context("Starting values for MCMC\n")

test_that("Generated starting values do not match given seed (spikes case)", {
  Data <- makeExampleBASiCS_Data(WithSpikes = TRUE)
  n <- ncol(Data)
  PriorParam <- list(s2.mu = 0.5, s2.delta = 0.5, a.delta = 1, 
                     b.delta = 1, p.phi = rep(1, times = n), 
                     a.s = 1, b.s = 1, a.theta = 1, b.theta = 1)
  set.seed(2018)
  Start <- BASiCS:::HiddenBASiCS_MCMC_Start(Data, PriorParam, WithSpikes = TRUE)
  
  Check0 <- c(11.763,  7.504,  5.791,  7.638, 22.870)
  Check <- as.vector(round(Start$mu0[1:5], 3))
  expect_that(all.equal(Check0, Check), is_true())
  
  Check0 <- c(1.552, 2.567, 2.010, 1.877, 0.682)
  Check <- as.vector(round(Start$delta0[1:5], 3))
  expect_that(all.equal(Check0, Check), is_true())
  
  Check0 <- c(1.068, 0.998, 1.058, 1.084, 0.982)
  Check <- as.vector(round(Start$phi0[1:5], 3))
  expect_that(all.equal(Check0, Check), is_true())
  
  Check0 <- c(0.189, 0.438, 0.043, 0.137, 0.417)
  Check <- as.vector(round(Start$s0[1:5], 3))
  expect_that(all.equal(Check0, Check), is_true())
  
  Check0 <- c(0.189, 0.438, 0.043, 0.137, 0.417)
  Check <- as.vector(round(Start$nu0[1:5], 3))
  expect_that(all.equal(Check0, Check), is_true())
  
  Check0 <- 0.817
  Check <- round(Start$theta0, 3)
  expect_that(all.equal(Check0, Check), is_true())
  
  Check0 <- c(-4, -4)
  Check <- as.vector(round(Start$ls.mu0[1:2], 3))
  expect_that(all.equal(Check0, Check), is_true())
  
  Check0 <- c(-2, -2)
  Check <- as.vector(round(Start$ls.delta0[1:2], 3))
  expect_that(all.equal(Check0, Check), is_true())
  
  Check0 <- 11
  Check <- round(Start$ls.phi0, 3)
  expect_that(all.equal(Check0, Check), is_true())
  
  Check0 <- c(-6.804, -8.210)
  Check <- as.vector(round(Start$ls.nu0[1:2], 3))
  expect_that(all.equal(Check0, Check), is_true())
  
  Check0 <- -4
  Check <- round(Start$ls.theta0, 3)
  expect_that(all.equal(Check0, Check), is_true())  
})