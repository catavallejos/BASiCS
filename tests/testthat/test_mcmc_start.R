context("Starting values for MCMC\n")

test_that("Generated starting values do not match given seed (spikes case)", {
  Data <- makeExampleBASiCS_Data(WithSpikes = TRUE)
  n <- ncol(Data)
  PriorParam <- list(s2.mu = 0.5, s2.delta = 0.5, a.delta = 1, 
                     b.delta = 1, p.phi = rep(1, times = n), 
                     a.s = 1, b.s = 1, a.theta = 1, b.theta = 1)
  set.seed(2018)
  Start <- BASiCS:::HiddenBASiCS_MCMC_Start(Data, PriorParam, WithSpikes = TRUE)

  Check0 <- c(11.566, 10.377,  6.010,  8.245, 25.929)
  Check <- as.vector(round(Start$mu0[1:5], 3))
  expect_that(all.equal(Check0, Check), is_true())
  
  Check0 <- c(1.769, 2.282, 2.876, 2.243, 0.858)
  Check <- as.vector(round(Start$delta0[1:5], 3))
  expect_that(all.equal(Check0, Check), is_true())
  
  Check0 <- c(1.109, 1.102, 1.202, 0.758, 0.651)
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

test_that("Generated starting values do not match given seed (no spikes case)", {
  Data <- makeExampleBASiCS_Data(WithSpikes = FALSE)
  n <- ncol(Data)
  PriorParam <- list(s2.mu = 0.5, s2.delta = 0.5, a.delta = 1, 
                     b.delta = 1, p.phi = rep(1, times = n), 
                     a.s = 1, b.s = 1, a.theta = 1, b.theta = 1)
  set.seed(2018)
  Start <- BASiCS:::HiddenBASiCS_MCMC_Start(Data, PriorParam, WithSpikes = FALSE)
  
  Check0 <- c(14.378, 11.894,  6.418, 10.961, 36.917)
  Check <- as.vector(round(Start$mu0[1:5], 3))
  expect_that(all.equal(Check0, Check), is_true())
  
  Check0 <- c(1.204, 1.281, 1.958, 1.448, 0.654)
  Check <- as.vector(round(Start$delta0[1:5], 3))
  expect_that(all.equal(Check0, Check), is_true())
  
  Check0 <- NULL
  Check <- Start$phi0
  expect_that(all.equal(Check0, Check), is_true())
  
  Check0 <- c(0.943, 1.295, 0.181, 0.410, 1.781)
  Check <- as.vector(round(Start$s0[1:5], 3))
  expect_that(all.equal(Check0, Check), is_true())
  
  Check0 <- c(0.943, 1.295, 0.181, 0.410, 1.781)
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
  
  Check0 <- c(-10, -10)
  Check <- as.vector(round(Start$ls.nu0[1:2], 3))
  expect_that(all.equal(Check0, Check), is_true())
  
  Check0 <- -4
  Check <- round(Start$ls.theta0, 3)
  expect_that(all.equal(Check0, Check), is_true())  
})

test_that("Generated starting values do not match given seed (regression+spike)", {
  Data <- makeExampleBASiCS_Data(WithSpikes = TRUE)
  n <- ncol(Data); k <- 12
  PriorParam <- list(s2.mu = 0.5, s2.delta = 0.5, a.delta = 1, 
                     b.delta = 1, p.phi = rep(1, times = n), 
                     a.s = 1, b.s = 1, a.theta = 1, b.theta = 1)
  PriorParam$m <- rep(0, k); PriorParam$V <- diag(k) 
  PriorParam$a.sigma2 <- 2; PriorParam$b.sigma2 <- 2  
  PriorParam$eta <- 5

  set.seed(2018)
  Start <- BASiCS:::HiddenBASiCS_MCMC_Start(Data, PriorParam, WithSpikes = TRUE)
  
  Check0 <- c(0.800,  1.939, -0.906,  0.554,  0.897)
  Check <- as.vector(round(Start$beta0[1:5], 3))
  expect_that(all.equal(Check0, Check), is_true())
  
  Check0 <- 2.117
  Check <- round(Start$sigma20, 3)
  expect_that(all.equal(Check0, Check), is_true())
  
  Check0 <- c(0.794, 0.914, 0.891, 0.237, 1.015)
  Check <- as.vector(round(Start$lambda0[1:5], 3))
  expect_that(all.equal(Check0, Check), is_true())
})