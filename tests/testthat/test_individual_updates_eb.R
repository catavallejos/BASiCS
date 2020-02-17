context("Individual MCMC updates (cpp code) - empirical bayes\n")


test_that("Spikes + no regression", {
  set.seed(3)
  Data <- makeExampleBASiCS_Data(WithSpikes = TRUE)
  CountsBio <- counts(Data)
  CountsAll <- rbind(CountsBio, assay(altExp(Data, "spike-ins")))
  q0 <- nrow(CountsBio)
  n <- ncol(CountsBio)
  PriorParam <- list(
    mu.mu = rep(0, q0), s2.mu = 0.5, s2.delta = 0.5, a.delta = 1,
    b.delta = 1, p.phi = rep(1, times = n),
    a.s = 1, b.s = 1, a.theta = 1, b.theta = 1
  )
  PriorParam <- BASiCS:::.EmpiricalBayesMu(Data, PriorParam)
  
  mu.mu0 <- c(0.43,  1.32,  0.24,  1.05,  1.95,  1.15)
  mu.mu1 <- as.vector(round(PriorParam$mu.mu[1:6], 2))
  expect_equal(mu.mu0, mu.mu1)
    
  set.seed(2018)
  Start <- BASiCS:::HiddenBASiCS_MCMC_Start(Data, PriorParam, WithSpikes = TRUE)
  uGene <- rep(0, times = q0)
  indGene <- rbinom(q0, size = 1, prob = 0.5)
  BatchDesign <- model.matrix(~ 0 + Data$BatchInfo)
  ThetaBatch <- BatchDesign %*% Start$theta0
  uCell <- rep(0, times = n)
  indCell <- rbinom(n, size = 1, prob = 0.5)

  mu1 <- pmax(0, Start$mu0[seq_len(q0)] + rnorm(q0, sd = 0.005))
  mu <- BASiCS:::.muUpdate(
    mu0 = Start$mu0,
    prop_var = exp(Start$ls.mu0),
    Counts = CountsBio,
    invdelta = 1 / Start$delta0,
    phinu = Start$phi0 * Start$nu0,
    sum_bycell_bio = rowSums(CountsBio),
    mu_mu = PriorParam$mu.mu,
    s2_mu = PriorParam$s2.mu,
    q0 = q0,
    n = n,
    mu1 = mu1,
    u = uGene,
    ind = indGene,
    exponent = 1,
    mintol = 1e-3
  )

  mu1 <- c(6.41, 15.67, 5.57, 11.80, 31.68, 12.08)
  mu1_obs <- round(mu[1:6, 1], 2)
  expect_equal(mu1, mu1_obs)

  ind <- c(1, 0, 1, 1, 1, 1)
  ind_obs <- mu[1:6, 2]
  expect_equal(ind, ind_obs)

})


test_that("Spikes + regression", {
  set.seed(3)
  Data <- makeExampleBASiCS_Data(WithSpikes = TRUE)
  CountsBio <- counts(Data)
  CountsAll <- rbind(CountsBio, assay(altExp(Data, "spike-ins")))
  q0 <- nrow(CountsBio)
  n <- ncol(CountsBio)
  k <- 12
  var <- 1.2
  PriorParam <- list(
    mu.mu = rep(0, q0), s2.mu = 0.5, s2.delta = 0.5, a.delta = 1,
    eta = 5, m = rep(0, k), V = diag(k),
    a.sigma2 = 1, b.sigma2 = 1,
    b.delta = 1, p.phi = rep(1, times = n),
    a.s = 1, b.s = 1, a.theta = 1, b.theta = 1
  )
  
  PriorParam <- BASiCS:::.EmpiricalBayesMu(Data, PriorParam)
  mu.mu0 <- c(0.43,  1.32,  0.24,  1.05,  1.95,  1.15)
  mu.mu1 <- as.vector(round(PriorParam$mu.mu[1:6], 2))
  expect_equal(mu.mu0, mu.mu1)
  
  set.seed(2020)
  Start <- BASiCS:::HiddenBASiCS_MCMC_Start(Data, PriorParam, WithSpikes = TRUE)
  uGene <- rep(0, times = q0)
  indGene <- rbinom(q0, size = 1, prob = 0.5)
  BatchDesign <- model.matrix(~ 0 + Data$BatchInfo)
  ThetaBatch <- BatchDesign %*% Start$theta0
  uCell <- rep(0, times = n)
  indCell <- rbinom(n, size = 1, prob = 0.5)
  X <- BASiCS:::.designMatrix(k, Start$mu0, var)

  # Hidden_muUpdate
  mu1 <- pmax(0, Start$mu0[seq_len(q0)] + rnorm(q0, sd = 0.005))
  mu <- BASiCS:::.muUpdateReg(
    mu0 = Start$mu0,
    prop_var = exp(Start$ls.mu0),
    Counts = CountsBio,
    delta = Start$delta0,
    phinu = Start$phi0 * Start$nu0,
    sum_bycell_bio = rowSums(CountsBio),
    mu_mu = PriorParam$mu.mu,
    s2_mu = PriorParam$s2.mu,
    q0 = q0,
    n = n,
    mu1 = mu1,
    u = uGene,
    ind = indGene,
    k = k,
    lambda = Start$lambda0,
    beta = Start$beta0,
    X = X,
    sigma2 = Start$sigma20,
    variance = var,
    exponent = 1,
    mintol = 1e-3
  )

  mu1 <- c(5.92, 15.67, 6.23, 12.26, 28.44)
  mu1_obs <- round(mu[1:5, 1], 2)
  expect_equal(mu1, mu1_obs)

  ind <- c(1, 0, 1, 0, 0)
  ind_obs <- mu[1:5, 2]
  expect_equal(ind, ind_obs)

})


test_that("No Spikes + no regression", {
  set.seed(5)
  Data <- makeExampleBASiCS_Data(WithSpikes = FALSE)
  CountsBio <- counts(Data)
  q0 <- nrow(CountsBio)
  n <- ncol(CountsBio)
  PriorParam <- list(
    mu.mu = rep(0, q0), s2.mu = 0.5, s2.delta = 0.5, a.delta = 1,
    b.delta = 1, p.phi = rep(1, times = n),
    a.s = 1, b.s = 1, a.theta = 1, b.theta = 1
  )
  PriorParam <- BASiCS:::.EmpiricalBayesMu(Data, PriorParam)
  mu.mu0 <- c(1.95, 2.75, 1.81, 2.32, 3.31, 2.77)
  mu.mu1 <- as.vector(round(PriorParam$mu.mu[1:6], 2))
  expect_equal(mu.mu0, mu.mu1)
  
  set.seed(2020)
  Start <- BASiCS:::HiddenBASiCS_MCMC_Start(Data, PriorParam, WithSpikes = FALSE)
  uGene <- rep(0, times = q0)
  indGene <- rbinom(q0, size = 1, prob = 0.5)
  BatchDesign <- model.matrix(~ 0 + Data$BatchInfo)
  ThetaBatch <- BatchDesign %*% Start$theta0
  uCell <- rep(0, times = n)
  indCell <- rbinom(n, size = 1, prob = 0.5)

  ## Components for no-spikes
  means <- rowMeans(CountsBio)
  RefGene <- which.min(means[means >= median(means)])

  # Hidden_muUpdate
  mu1 <- pmax(0, Start$mu0[seq_len(q0)] + rnorm(q0, sd = 0.005))
  mu <- BASiCS:::.muUpdateNoSpikes(
    mu0 = Start$mu0,
    prop_var = exp(Start$ls.mu0),
    Counts = CountsBio,
    invdelta = 1 / Start$delta0,
    nu = Start$nu0,
    sum_bycell_all = rowSums(CountsBio),
    mu_mu = PriorParam$mu.mu,
    s2_mu = PriorParam$s2.mu,
    q0 = q0,
    n = n,
    mu1 = mu1,
    u = uGene,
    ind = indGene,
    Constrain = mean(log(Start$mu0)),
    RefGene = RefGene,
    ConstrainGene = seq_len(q0) - 1,
    NotConstrainGene = 0,
    ConstrainType = 1,
    exponent = 1,
    mintol = 1e-3
  )


  mu1 <- c(6.98, 24.26,  8.12, 13.85, 35.03)
  mu1_obs <- round(mu[1:5, 1], 2)
  expect_equal(mu1, mu1_obs)

  ind <- c(1, 1, 1, 1, 0)
  ind_obs <- mu[1:5, 2]
  expect_equal(ind, ind_obs)

})



test_that("No Spikes + regression", {
  set.seed(12)
  Data <- makeExampleBASiCS_Data(WithSpikes = FALSE)
  CountsBio <- counts(Data)
  q0 <- nrow(CountsBio)
  n <- ncol(CountsBio)
  k <- 12
  var <- 1.2
  PriorParam <- list(
    mu.mu = rep(0, q0), s2.mu = 0.5, s2.delta = 0.5, a.delta = 1,
    eta = 5, m = rep(0, k), V = diag(k),
    a.sigma2 = 1, b.sigma2 = 1,
    b.delta = 1, p.phi = rep(1, times = n),
    a.s = 1, b.s = 1, a.theta = 1, b.theta = 1
  )
  PriorParam <- BASiCS:::.EmpiricalBayesMu(Data, PriorParam)
  mu.mu0 <- c(2.00, 2.73, 1.81, 1.69, 3.49, 2.75)
  mu.mu1 <- as.vector(round(PriorParam$mu.mu[1:6], 2))
  expect_equal(mu.mu0, mu.mu1)

  set.seed(2044)
  Start <- BASiCS:::HiddenBASiCS_MCMC_Start(Data, PriorParam, WithSpikes = FALSE)
  uGene <- rep(0, times = q0)
  indGene <- rbinom(q0, size = 1, prob = 0.5)
  BatchDesign <- model.matrix(~ 0 + Data$BatchInfo)
  ThetaBatch <- BatchDesign %*% Start$theta0
  uCell <- rep(0, times = n)
  indCell <- rbinom(n, size = 1, prob = 0.5)
  X <- BASiCS:::.designMatrix(k, Start$mu0, var)

  ## Components for no-spikes
  means <- rowMeans(CountsBio)
  RefGene <- which.min(means[means >= median(means)])

  mu1 <- pmax(0, Start$mu0[seq_len(q0)] + rnorm(q0, sd = 0.005))
  mu <- BASiCS:::.muUpdateRegNoSpikes(
    mu0 = Start$mu0,
    prop_var = exp(Start$ls.mu0),
    Counts = CountsBio,
    delta = Start$delta0,
    invdelta = 1 / Start$delta0,
    nu = Start$nu0,
    sum_bycell_all = rowSums(CountsBio),
    mu_mu = PriorParam$mu.mu,
    s2_mu = PriorParam$s2.mu,
    q0 = q0,
    n = n,
    mu1 = mu1,
    u = uGene,
    ind = indGene,
    Constrain = mean(log(Start$mu0)),
    RefGene = RefGene,
    ConstrainGene = seq_len(q0) - 1,
    NotConstrainGene = 0,
    ConstrainType = 1,
    k = k,
    lambda = Start$lambda0,
    beta = Start$beta0,
    X = X,
    sigma2 = Start$sigma20,
    variance = var,
    exponent = 1,
    mintol = 1e-3
  )

  mu1 <- c(9.49, 17.42,  6.38, 10.30, 36.60)
  mu1_obs <- round(mu[1:5, 1], 2)
  expect_equal(mu1, mu1_obs)

  ind <- c(0, 1, 1, 1, 1)
  ind_obs <- mu[1:5, 2]
  expect_equal(ind, ind_obs)

})
