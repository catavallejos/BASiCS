context("Starting values for MCMC\n")

test_that("Generated starting values do not match given seed (spikes case)", {
  set.seed(4)
  Data <- makeExampleBASiCS_Data(WithSpikes = TRUE)
  n <- ncol(Data)
  PriorParam <- list(s2.mu = 0.5, s2.delta = 0.5, a.delta = 1,
                     b.delta = 1, p.phi = rep(1, times = n),
                     a.s = 1, b.s = 1, a.theta = 1, b.theta = 1)
  set.seed(2018)
  Start <- BASiCS:::HiddenBASiCS_MCMC_Start(Data, PriorParam, WithSpikes = TRUE,
    Regression = FALSE)

  Check0 <- c(7.262, 10.828,  6.816,  8.258, 21.141)
  Check <- as.vector(round(Start$mu0[1:5], 3))
  expect_equal(Check0, Check)
  
  Check0 <- c(1.542, 1.229, 2.055, 1.688, 0.662)
  Check <- as.vector(round(Start$delta0[1:5], 3))
  expect_equal(Check0, Check)
  
  Check0 <- c(0.781, 1.265, 1.586, 1.139, 1.036)
  Check <- as.vector(round(Start$phi0[1:5], 3))
  expect_equal(Check0, Check)
  
  Check0 <- c(0.331, 0.185, 0.152, 0.812, 0.411)
  Check <- as.vector(round(Start$s0[1:5], 3))
  expect_equal(Check0, Check)
  
  Check0 <- c(0.331, 0.185, 0.152, 0.812, 0.411)
  Check <- as.vector(round(Start$nu0[1:5], 3))
  expect_equal(Check0, Check)
  
  Check0 <- 0.817
  Check <- round(Start$theta0, 3)
  expect_equal(Check0, Check)
  
  Check0 <- c(-4, -4)
  Check <- as.vector(round(Start$ls.mu0[1:2], 3))
  expect_equal(Check0, Check)
  
  Check0 <- c(-2, -2)
  Check <- as.vector(round(Start$ls.delta0[1:2], 3))
  expect_equal(Check0, Check)
  
  Check0 <- 11
  Check <- round(Start$ls.phi0, 3)
  expect_equal(Check0, Check)
  
  Check0 <- c(-7.624, -6.779)
  Check <- as.vector(round(Start$ls.nu0[1:2], 3))
  expect_equal(Check0, Check)
  
  Check0 <- -4
  Check <- round(Start$ls.theta0, 3)
  expect_equal(Check0, Check)
})

test_that("Generated starting values do not match given seed (no spikes case)", {

  set.seed(5)
  Data <- makeExampleBASiCS_Data(WithSpikes = FALSE)
  n <- ncol(Data)
  PriorParam <- list(s2.mu = 0.5, s2.delta = 0.5, a.delta = 1,
                     b.delta = 1, p.phi = rep(1, times = n),
                     a.s = 1, b.s = 1, a.theta = 1, b.theta = 1)
  set.seed(2018)
  Start <- BASiCS:::HiddenBASiCS_MCMC_Start(Data, PriorParam, WithSpikes = FALSE,
    Regression = FALSE)
  
  Check0 <- c(9.019, 20.049,  7.819, 13.087, 35.027)
  Check <- as.vector(round(Start$mu0[1:5], 3))
  expect_equal(Check0, Check)
  
  Check0 <- c(0.868, 1.964, 1.711, 1.275, 0.768)
  Check <- as.vector(round(Start$delta0[1:5], 3))
  expect_equal(Check0, Check)
  
  Check0 <- NULL
  Check <- Start$phi0
  expect_equal(Check0, Check)
  
  Check0 <- c(0.446, 0.725, 1.309, 2.570, 0.511)
  Check <- as.vector(round(Start$s0[1:5], 3))
  expect_equal(Check0, Check)
  
  Check0 <- c(0.446, 0.725, 1.309, 2.570, 0.511)
  Check <- as.vector(round(Start$nu0[1:5], 3))
  expect_equal(Check0, Check)
  
  Check0 <- 0.817
  Check <- round(Start$theta0, 3)
  expect_equal(Check0, Check)
  
  Check0 <- c(-4, -4)
  Check <- as.vector(round(Start$ls.mu0[1:2], 3))
  expect_equal(Check0, Check)
  
  Check0 <- c(-2, -2)
  Check <- as.vector(round(Start$ls.delta0[1:2], 3))
  expect_equal(Check0, Check)
  
  Check0 <- 11
  Check <- round(Start$ls.phi0, 3)
  expect_equal(Check0, Check)
  
  Check0 <- c(-8.254, -10.000)
  Check <- as.vector(round(Start$ls.nu0[1:2], 3))
  expect_equal(Check0, Check)
  
  Check0 <- -4
  Check <- round(Start$ls.theta0, 3)
  expect_equal(Check0, Check)
})

test_that("Generated starting values do not match given seed (regression+spike)", {

  set.seed(6)
  Data <- makeExampleBASiCS_Data(WithSpikes = TRUE)
  n <- ncol(Data); k <- 12
  PriorParam <- list(s2.mu = 0.5, s2.delta = 0.5, a.delta = 1, 
                     b.delta = 1, p.phi = rep(1, times = n),
                     a.s = 1, b.s = 1, a.theta = 1, b.theta = 1)
  PriorParam$m <- rep(0, k); PriorParam$V <- diag(k)
  PriorParam$a.sigma2 <- 2; PriorParam$b.sigma2 <- 2
  PriorParam$eta <- 5

  set.seed(2018)
  Start <- BASiCS:::HiddenBASiCS_MCMC_Start(Data, PriorParam, WithSpikes = TRUE,
    Regression = TRUE)
  
  Check0 <- c(0.800,  1.939, -0.906,  0.554,  0.897)
  Check <- as.vector(round(Start$beta0[1:5], 3))
  expect_equal(Check0, Check)
  
  Check0 <- 2.117
  Check <- round(Start$sigma20, 3)
  expect_equal(Check0, Check)

  Check0 <- c(0.794, 0.914, 0.891, 0.237, 1.015)
  Check <- as.vector(round(Start$lambda0[1:5], 3))
  expect_equal(Check0, Check)
})
