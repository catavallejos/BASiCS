context("Individual MCMC updates (cpp code) - empirical bayes")

test_that("Spikes + no regression", {
  set.seed(3)
  Data <- makeExampleBASiCS_Data(WithSpikes = TRUE)
  CountsBio <- counts(Data)
  CountsAll <- rbind(CountsBio, assay(altExp(Data, "spike-ins")))
  q0 <- nrow(CountsBio)
  n <- ncol(CountsBio)
  PriorParam <- BASiCS_PriorParam(Data, PriorMu = "EmpiricalBayes")
  PriorParam$mu.mu <- BASiCS:::.EmpiricalBayesMu(Data, PriorParam$s2.mu, with_spikes=TRUE)
  
  mu.mu0 <- c(2.19, 3.07, 2, 2.81, 3.7, 2.91)
  mu.mu1 <- as.vector(round(PriorParam$mu.mu[1:6], 2))
  expect_equal(mu.mu0, mu.mu1)

  set.seed(2018)
  Start <- BASiCS:::.BASiCS_MCMC_Start(Data, PriorParam, Regression = FALSE, WithSpikes = TRUE)
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

  mu1 <- c(8.12, 21.64, 6.85, 15.99, 40.49, 16.49)
  mu1_obs <- round(mu[1:6, 1], 2)
  expect_equal(mu1, mu1_obs)

  ind <- c(1, 0, 1, 1, 0, 1)
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
  PriorParam <- BASiCS_PriorParam(Data, PriorMu = "EmpiricalBayes")
  PriorParam$mu.mu <- BASiCS:::.EmpiricalBayesMu(Data, PriorParam$s2.mu, with_spikes=TRUE)
  mu.mu0 <- c(2.19, 3.07, 2, 2.81, 3.7, 2.91)
  mu.mu1 <- as.vector(round(PriorParam$mu.mu[1:6], 2))
  expect_equal(mu.mu0, mu.mu1)
  
  set.seed(2020)
  Start <- BASiCS:::.BASiCS_MCMC_Start(Data, PriorParam, Regression = TRUE, WithSpikes = TRUE)
  uGene <- rep(0, times = q0)
  indGene <- rbinom(q0, size = 1, prob = 0.5)
  BatchDesign <- model.matrix(~ 0 + Data$BatchInfo)
  ThetaBatch <- BatchDesign %*% Start$theta0
  uCell <- rep(0, times = n)
  indCell <- rbinom(n, size = 1, prob = 0.5)
  RBFLocations <- BASiCS:::.estimateRBFLocations(log(Start$mu0), k, TRUE)
  X <- BASiCS:::.designMatrix(k, RBFLocations, Start$mu0, var)

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
    FixLocations = FALSE,
    RBFMinMax = TRUE,
    RBFLocations = RBFLocations,
    variance = var,
    exponent = 1,
    mintol = 1e-3
  )

  mu1 <- c(7.49, 21.64, 7.66, 16.68, 40.49)
  mu1_obs <- round(mu[1:5, 1], 2)
  expect_equal(mu1, mu1_obs)

  ind <- c(1, 0, 1, 1, 0)
  ind_obs <- mu[1:5, 2]
  expect_equal(ind, ind_obs)
})

test_that("No Spikes + no regression", {
  set.seed(5)
  Data <- makeExampleBASiCS_Data(WithSpikes = FALSE)
  CountsBio <- counts(Data)
  q0 <- nrow(CountsBio)
  n <- ncol(CountsBio)
  PriorParam <- BASiCS_PriorParam(Data, PriorMu = "EmpiricalBayes")
  PriorParam$mu.mu <- BASiCS:::.EmpiricalBayesMu(Data, PriorParam$s2.mu, with_spikes=FALSE)
  mu.mu0 <- c(1.95, 2.75, 1.81, 2.32, 3.31, 2.77)
  mu.mu1 <- as.vector(round(PriorParam$mu.mu[1:6], 2))
  expect_equal(mu.mu0, mu.mu1)
  
  set.seed(2020)
  Start <- BASiCS:::.BASiCS_MCMC_Start(
    Data,
    PriorParam,
    WithSpikes = FALSE,
    Regression = FALSE
  )
  uGene <- rep(0, times = q0)
  indGene <- rbinom(q0, size = 1, prob = 0.5)
  BatchDesign <- model.matrix(~ 0 + Data$BatchInfo)
  ThetaBatch <- BatchDesign %*% Start$theta0
  uCell <- rep(0, times = n)
  indCell <- rbinom(n, size = 1, prob = 0.5)

  ## Components for no-spikes
  means <- rowMeans(CountsBio)
  RefGene <- which.min(means[means >= median(means)])

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
    exponent = 1,
    mintol = 1e-3
  )

  mu1 <- c(7.02, 18.89,  6.32, 10.78, 29.61)
  mu1_obs <- round(mu[1:5, 1], 2)
  expect_equal(mu1, mu1_obs)

  ind <- c(0, 1, 1, 1, 1)
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
  PriorParam <- BASiCS_PriorParam(Data, a.sigma2 = 1, b.sigma2 = 1, PriorMu = "EmpiricalBayes")
  PriorParam$mu.mu <- BASiCS:::.EmpiricalBayesMu(Data, PriorParam$s2.mu, with_spikes=FALSE)
  mu.mu0 <- c(2.00, 2.73, 1.81, 1.69, 3.49, 2.75)
  mu.mu1 <- as.vector(round(PriorParam$mu.mu[1:6], 2))
  expect_equal(mu.mu0, mu.mu1)

  set.seed(2044)
  Start <- BASiCS:::.BASiCS_MCMC_Start(
    Data,
    PriorParam,
    WithSpikes = FALSE,
    Regression = TRUE
  )
  uGene <- rep(0, times = q0)
  indGene <- rbinom(q0, size = 1, prob = 0.5)
  BatchDesign <- model.matrix(~ 0 + Data$BatchInfo)
  ThetaBatch <- BatchDesign %*% Start$theta0
  uCell <- rep(0, times = n)
  indCell <- rbinom(n, size = 1, prob = 0.5)
  RBFLocations <- BASiCS:::.estimateRBFLocations(log(Start$mu0), k, TRUE)
  X <- BASiCS:::.designMatrix(k, RBFLocations, Start$mu0, var)

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
    k = k,
    lambda = Start$lambda0,
    beta = Start$beta0,
    X = X,
    sigma2 = Start$sigma20,
    variance = var,
    FixLocations = FALSE,
    RBFMinMax = TRUE,
    RBFLocations = RBFLocations,
    exponent = 1,
    mintol = 1e-3
  )

  mu1 <- c(7.39, 13.57,  6.13,  1.79, 32.77)
  mu1_obs <- round(mu[1:5, 1], 2)
  expect_equal(mu1, mu1_obs)

  ind <- c(0, 1, 0, 1, 0)
  ind_obs <- mu[1:5, 2]
  expect_equal(ind, ind_obs)
})
