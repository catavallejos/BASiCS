test_that("Generated starting values do not match given seed (spikes case)", {
  set.seed(4)
  Data <- makeExampleBASiCS_Data(WithSpikes = TRUE)
  n <- ncol(Data)
  PriorParam <- BASiCS_PriorParam(Data, k = 12)
  set.seed(2018)
  Start <- BASiCS:::.BASiCS_MCMC_Start(
    Data,
    PriorParam,
    WithSpikes = TRUE,
    Regression = FALSE
  )

  muCheck0 <- c(7.262, 10.828,  6.816,  8.258, 21.141)
  muCheck <- as.vector(round(Start$mu0[1:5], 3))
  expect_equal(muCheck0, muCheck)
  
  deltaCheck0 <- c(1.542, 1.229, 2.055, 1.688, 0.662)
  deltaCheck <- as.vector(round(Start$delta0[1:5], 3))
  expect_equal(deltaCheck0, deltaCheck)
  
  phiCheck0 <- c(0.781, 1.265, 1.586, 1.139, 1.036)
  phiCheck <- as.vector(round(Start$phi0[1:5], 3))
  expect_equal(phiCheck0, phiCheck)
  
  sCheck0 <- c(0.331, 0.185, 0.152, 0.812, 0.411)
  sCheck <- as.vector(round(Start$s0[1:5], 3))
  expect_equal(sCheck0, sCheck)
  
  nuCheck0 <- c(0.331, 0.185, 0.152, 0.812, 0.411)
  nuCheck <- as.vector(round(Start$nu0[1:5], 3))
  expect_equal(nuCheck0, nuCheck)
  
  thetaCheck0 <- 0.817
  thetaCheck <- round(Start$theta0, 3)
  expect_equal(thetaCheck0, thetaCheck)
  
  lsmuCheck0 <- c(-4, -4)
  lsmuCheck <- as.vector(round(Start$ls.mu0[1:2], 3))
  expect_equal(lsmuCheck0, lsmuCheck)
  
  lsdeltaCheck0 <- c(-2, -2)
  lsdeltaCheck <- as.vector(round(Start$ls.delta0[1:2], 3))
  expect_equal(lsdeltaCheck0, lsdeltaCheck)
  
  lsphiCheck0 <- 11
  lsphiCheck <- round(Start$ls.phi0, 3)
  expect_equal(lsphiCheck0, lsphiCheck)
  
  lsnuCheck0 <- c(-7.624, -6.779)
  lsnuCheck <- as.vector(round(Start$ls.nu0[1:2], 3))
  expect_equal(lsnuCheck0, lsnuCheck)
  
  lsthetaCheck0 <- -4
  lsthetaCheck <- round(Start$ls.theta0, 3)
  expect_equal(lsthetaCheck0, lsthetaCheck)
})

test_that("Generated starting values do not match given seed (no spikes case)", {

  set.seed(5)
  Data <- makeExampleBASiCS_Data(WithSpikes = FALSE)
  n <- ncol(Data)
  PriorParam <- BASiCS_PriorParam(Data, k = 12)

  set.seed(2018)
  Start <- BASiCS:::.BASiCS_MCMC_Start(
    Data,
    PriorParam,
    WithSpikes = FALSE,
    Regression = FALSE
  )
  
  muCheck0 <- c(9.019, 20.049,  7.819, 13.087, 35.027)
  muCheck <- as.vector(round(Start$mu0[1:5], 3))
  expect_equal(muCheck0, muCheck)
  
  deltaCheck0 <- c(0.868, 1.964, 1.711, 1.275, 0.768)
  deltaCheck <- as.vector(round(Start$delta0[1:5], 3))
  expect_equal(deltaCheck0, deltaCheck)
  
  phiCheck0 <- NULL
  phiCheck <- Start$phi0
  expect_equal(phiCheck0, phiCheck)
  
  sCheck0 <- c(0.446, 0.725, 1.309, 2.570, 0.511)
  sCheck <- as.vector(round(Start$s0[1:5], 3))
  expect_equal(sCheck0, sCheck)
  
  nuCheck0 <- c(0.446, 0.725, 1.309, 2.570, 0.511)
  nuCheck <- as.vector(round(Start$nu0[1:5], 3))
  expect_equal(nuCheck0, nuCheck)
  
  thetaCheck0 <- 0.817
  thetaCheck <- round(Start$theta0, 3)
  expect_equal(thetaCheck0, thetaCheck)
  
  lsmuCheck0 <- c(-4, -4)
  lsmuCheck <- as.vector(round(Start$ls.mu0[1:2], 3))
  expect_equal(lsmuCheck0, lsmuCheck)
  
  lsdeltaCheck0 <- c(-2, -2)
  lsdeltaCheck <- as.vector(round(Start$ls.delta0[1:2], 3))
  expect_equal(lsdeltaCheck0, lsdeltaCheck)
  
  lsphiCheck0 <- 11
  lsphiCheck <- round(Start$ls.phi0, 3)
  expect_equal(lsphiCheck0, lsphiCheck)
  
  lsnuCheck0 <- c(-8.254, -10.000)
  lsnuCheck <- as.vector(round(Start$ls.nu0[1:2], 3))
  expect_equal(lsnuCheck0, lsnuCheck)
  
  lsthetaCheck0 <- -4
  lsthetaCheck <- round(Start$ls.theta0, 3)
  expect_equal(lsthetaCheck0, lsthetaCheck)
})

test_that("Generated starting values do not match given seed (regression+spike)", {

  set.seed(6)
  Data <- makeExampleBASiCS_Data(WithSpikes = TRUE)
  PriorParam <- BASiCS_PriorParam(Data, k = 12)
  
  set.seed(2018)
  Start <- BASiCS:::.BASiCS_MCMC_Start(
    Data,
    PriorParam,
    WithSpikes = TRUE,
    Regression = TRUE
  )
  
  betaCheck0 <- c(0.800,  1.939, -0.906,  0.554,  0.897)
  betaCheck <- as.vector(round(Start$beta0[1:5], 3))
  expect_equal(betaCheck0, betaCheck)
  
  s2Check0 <- 2.117
  s2Check <- round(Start$sigma20, 3)
  expect_equal(s2Check0, s2Check)

  lambdaCheck0 <- c(0.794, 0.914, 0.891, 0.237, 1.015)
  lambdaCheck <- as.vector(round(Start$lambda0[1:5], 3))
  expect_equal(lambdaCheck0, lambdaCheck)
})
