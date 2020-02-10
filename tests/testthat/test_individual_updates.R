context("Individual MCMC updates (cpp code)\n")

test_that("Dirichlet sampler", {
  # Generating arbitrary hyper-param
  phi0 <- 1:10
  phi0 <- phi0 / sum(phi0)
  set.seed(2018)
  x <- as.vector(BASiCS:::Hidden_rDirichlet(phi0))

  set.seed(2018)
  x0 <- rgamma(length(phi0), shape = phi0, scale = 1)
  x0 <- x0 / sum(x0)

  expect_equal(x, x0)
})

test_that("Spikes + no regression", {
  set.seed(3)
  Data <- makeExampleBASiCS_Data(WithSpikes = TRUE)
  CountsBio <- counts(Data)
  CountsAll <- rbind(CountsBio, assay(altExp(Data, "spike-ins")))
  q0 <- nrow(CountsBio)
  n <- ncol(CountsBio)
  PriorParam <- list(
    mu.mu = 0, s2.mu = 0.5, s2.delta = 0.5, a.delta = 1,
    b.delta = 1, p.phi = rep(1, times = n),
    a.s = 1, b.s = 1, a.theta = 1, b.theta = 1
  )
  set.seed(2018)
  Start <- BASiCS:::HiddenBASiCS_MCMC_Start(Data, PriorParam, WithSpikes = TRUE,
    Regression = FALSE)
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

  mu1 <- c(6.41, 15.67, 5.57, 11.80, 28.44)
  mu1_obs <- round(mu[1:5, 1], 2)
  expect_equal(mu1, mu1_obs)

  ind <- c(1, 0, 1, 1, 0)
  ind_obs <- mu[1:5, 2]
  expect_equal(ind, ind_obs)

  delta1 <- pmax(0, Start$delta0[seq_len(q0)] + rnorm(q0, sd = 0.005))
  delta_prior1 <- BASiCS:::.deltaUpdate(
    delta0 = Start$delta0,
    prop_var = exp(Start$ls.delta0),
    Counts = CountsBio,
    mu = Start$mu0,
    phinu = Start$phi0 * Start$nu0,
    a_delta = PriorParam$a.delta,
    b_delta = PriorParam$b.delta,
    s2delta = PriorParam$s2.delta,
    prior_delta = 1,
    q0 = q0,
    n = n,
    delta1 = delta1,
    u = uGene,
    ind = indGene,
    exponent = 1,
    mintol = 1e-3
  )

  delta1_prior1 <- c(0.73, 0.6, 3.52, 1.2, 0.69)
  delta1_prior1_obs <- round(delta_prior1[1:5, 1], 2)
  expect_equal(delta1_prior1, delta1_prior1_obs)

  ind <- c(1, 1, 1, 1, 0)
  ind_obs <- delta_prior1[1:5, 2]
  expect_equal(ind, ind_obs)

  delta_prior2 <- BASiCS:::.deltaUpdate(
    delta0 = Start$delta0,
    prop_var = exp(Start$ls.delta0),
    Counts = CountsBio,
    mu = Start$mu0,
    phinu = Start$phi0 * Start$nu0,
    a_delta = PriorParam$a.delta,
    b_delta = PriorParam$b.delta,
    s2delta = PriorParam$s2.delta,
    prior_delta = 2,
    q0 = q0,
    n = n,
    delta1 = delta1,
    u = uGene,
    ind = indGene,
    exponent = 1,
    mintol = 1e-3
  )

  delta1_prior2 <- c(1.05, 0.75, 1.71, 2.56, 0.65)
  delta1_prior2_obs <- round(delta_prior2[1:5, 1], 2)
  expect_equal(delta1_prior2, delta1_prior2_obs)

  ind <- c(1, 0, 0, 1, 1)
  ind_obs <- delta_prior2[1:5, 2]
  expect_equal(ind, ind_obs)

  phi1 <- rgamma(n, shape = 0.1, scale = 1)
  phi1 <- phi1 / sum(phi1) * length(phi1)
  phi <- BASiCS:::.phiUpdate(
    phi0 = rep(1, n),
    prop_var = 200,
    Counts = CountsBio,
    mu = Start$mu0,
    invdelta = 1 / Start$delta0,
    nu = Start$nu0,
    aphi = PriorParam$p.phi,
    sum_bygene_bio = colSums(CountsBio),
    q0 = q0,
    n = n,
    phi1 = phi1,
    exponent = 1
  )
  phi1 <- c(0.84, 1.01, 0.94, 1.05, 1.07)
  phi1_obs <- round(phi[[1]][1:5], 2)
  expect_equal(phi1, phi1_obs)

  ind <- 1
  ind_obs <- phi[[2]]
  expect_equal(ind, ind_obs)

  s1 <- pmax(0, Start$s0[seq_len(n)] + rnorm(n, sd = 0.005))
  s <- BASiCS:::.sUpdateBatch(
    s0 = Start$s0,
    nu = Start$nu0,
    thetaBatch = ThetaBatch,
    as = PriorParam$a.s,
    bs = PriorParam$b.s,
    BatchDesign = BatchDesign,
    n = n,
    s1 = s1,
    exponent = 1
  )

  s1 <- c(0.22, 0.12, 0.25, 0.32, 0.41)
  s1Obs <- round(s[1:5], digits = 2)
  expect_equal(s1, s1Obs)


  nu1 <- pmax(0, Start$nu0[seq_len(n)] + rnorm(n, sd = 0.005))
  nu <- BASiCS:::.nuUpdateBatch(
    nu0 = Start$nu0,
    prop_var = exp(Start$ls.nu0),
    Counts = CountsAll,
    SumSpikeInput = sum(metadata(Data)$SpikeInput[, 2]),
    BatchDesign = BatchDesign,
    mu = Start$mu0,
    invdelta = 1 / Start$delta0,
    phi = Start$phi0,
    s = Start$s0,
    thetaBatch = ThetaBatch,
    sum_bygene_all = colSums(CountsAll),
    q0 = q0,
    n = n,
    nu1 = nu1,
    u = uCell,
    ind = indCell,
    exponent = 1,
    mintol = 1e-3
  )

  nu1 <- c(0.11, 0.2, 0.07, 0.38, 0.25)
  nu1_obs <- round(nu[1:5, 1], 2)
  expect_equal(nu1, nu1_obs)

  ind <- c(1, 0, 1, 0, 1)
  ind_obs <- nu[1:5, 2]
  expect_equal(ind, ind_obs)

  theta <- BASiCS:::.thetaUpdateBatch(
    theta0 = Start$theta0,
    prop_var = exp(Start$ls.theta0),
    BatchDesign = BatchDesign,
    BatchSizes = colSums(BatchDesign),
    s = Start$s0,
    nu = Start$nu0,
    a_theta = PriorParam$a.theta,
    b_theta = PriorParam$b.theta,
    n = n,
    nBatch = ncol(BatchDesign),
    exponent = 1,
    mintol = 1e-3
  )
  theta1 <- 0.84
  theta1_obs <- round(theta[1, 1], 2)
  expect_equal(theta1, theta1_obs)

  ind <- 1
  ind_obs <- theta[1, 2]
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
    mu.mu = 0, s2.mu = 0.5, s2.delta = 0.5, a.delta = 1,
    eta = 5, m = rep(0, k), V = diag(k),
    a.sigma2 = 1, b.sigma2 = 1,
    b.delta = 1, p.phi = rep(1, times = n),
    a.s = 1, b.s = 1, a.theta = 1, b.theta = 1
  )

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

  delta1 <- pmax(0, Start$delta0[seq_len(q0)] + rnorm(q0, sd = 0.005))
  delta <- BASiCS:::.deltaUpdateReg(
    delta0 = Start$delta0,
    prop_var = exp(Start$ls.delta0),
    Counts = CountsBio,
    mu = Start$mu0,
    phinu = Start$phi0 * Start$nu0,
    q0 = q0,
    n = n,
    delta1 = delta1,
    u = uGene,
    ind = indGene,
    lambda = Start$lambda0,
    X = X,
    sigma2 = Start$sigma20,
    beta = Start$beta0,
    exponent = 1,
    mintol = 1e-3
  )

  delta1 <- c(1.39, 0.86, 1.23, 1.5, 0.8)
  delta1_obs <- round(delta[1:5, 1], 2)
  expect_equal(delta1, delta1_obs)

  ind <- c(1, 1, 0, 1, 1)
  ind_obs <- delta[1:5, 2]
  expect_equal(ind, ind_obs)


  beta <- BASiCS:::.betaUpdateReg(
    sigma2 = Start$sigma2,
    VAux = PriorParam$V,
    mAux = PriorParam$m
  )
  beta1 <- c(0.05, 1.02, 0.48, -0.14, 1.4)
  beta1_obs <- round(beta[1:5, 1], 2)
  expect_equal(beta1, beta1_obs)


  sigma2 <- BASiCS:::.sigma2UpdateReg(
    delta = Start$delta0,
    beta = Start$beta0,
    lambda = Start$lambda0,
    V1 = PriorParam$V,
    mInvVm0 = 0,
    m = PriorParam$m,
    sigma2_a0 = PriorParam$a.sigma2,
    sigma2_b0 = PriorParam$b.sigma2,
    q0 = q0,
    exponent = 1
  )
  expect_equal(round(sigma2, digits = 3), 0.489)

  lambda <- BASiCS:::.lambdaUpdateReg(
    delta = Start$delta0,
    X = X,
    beta = Start$beta0,
    sigma2 = Start$sigma20,
    eta = PriorParam$eta,
    q0 = q0,
    lambda1 = Start$lambda0,
    exponent = 1
  )
  lambda1 <- c(0.14, 0.13, 0.57, 0.28, 0.06)
  lambda1_obs <- round(lambda[1:5, 1], 2)
  expect_equal(lambda1, lambda1_obs)
})




test_that("No Spikes + no regression", {
  set.seed(5)
  Data <- makeExampleBASiCS_Data(WithSpikes = FALSE)
  CountsBio <- counts(Data)
  q0 <- nrow(CountsBio)
  n <- ncol(CountsBio)
  PriorParam <- list(
    mu.mu = 0, s2.mu = 0.5, s2.delta = 0.5, a.delta = 1,
    b.delta = 1, p.phi = rep(1, times = n),
    a.s = 1, b.s = 1, a.theta = 1, b.theta = 1
  )
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


  mu1 <- c(9.02, 24.26, 8.12, 13.85, 38.03)
  mu1_obs <- round(mu[1:5, 1], 2)
  expect_equal(mu1, mu1_obs)

  ind <- c(0, 1, 1, 1, 1)
  ind_obs <- mu[1:5, 2]
  expect_equal(ind, ind_obs)

  nu1 <- pmax(0, Start$nu0[seq_len(n)] + rnorm(n, sd = 0.005))
  nu <- BASiCS:::.nuUpdateBatchNoSpikes(
    nu0 = Start$nu0,
    prop_var = exp(Start$ls.nu0),
    Counts = CountsBio,
    BatchDesign = BatchDesign,
    mu = Start$mu0,
    invdelta = 1 / Start$delta0,
    s = Start$s0,
    thetaBatch = ThetaBatch,
    sum_bygene_all = colSums(CountsBio),
    q0 = q0,
    n = n,
    nu1 = nu1,
    u = uCell,
    ind = indCell,
    exponent = 1,
    mintol = 1e-3
  )

  nu1 <- c(0.43, 0.73, 1.31, 2.63, 0.52)
  nu1_obs <- round(nu[1:5, 1], 2)
  expect_equal(nu1, nu1_obs)

  ind <- c(1, 1, 0, 1, 1)
  ind_obs <- nu[1:5, 2]
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
    mu.mu = 0, s2.mu = 0.5, s2.delta = 0.5, a.delta = 1,
    eta = 5, m = rep(0, k), V = diag(k),
    a.sigma2 = 1, b.sigma2 = 1,
    b.delta = 1, p.phi = rep(1, times = n),
    a.s = 1, b.s = 1, a.theta = 1, b.theta = 1
  )

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

  mu1 <- c(7.1, 17.42, 6.38, 22.63, 36.6)
  mu1_obs <- round(mu[1:5, 1], 2)
  expect_equal(mu1, mu1_obs)

  ind <- c(1, 1, 1, 1, 1)
  ind_obs <- mu[1:5, 2]
  expect_equal(ind, ind_obs)

  delta1 <- pmax(0, Start$delta0[seq_len(q0)] + rnorm(q0, sd = 0.005))
  delta <- BASiCS:::.deltaUpdateRegNoSpikes(
    delta0 = Start$delta0,
    prop_var = exp(Start$ls.delta0),
    Counts = CountsBio,
    mu = Start$mu0,
    nu = Start$nu0,
    q0 = q0,
    n = n,
    delta1 = delta1,
    u = uGene,
    ind = indGene,
    lambda = Start$lambda0,
    X = X,
    sigma2 = Start$sigma20,
    beta = Start$beta0,
    exponent = 1,
    mintol = 1e-3
  )

  delta1 <- c(0.67, 0.51, 1.27, 2.28, 0.74)
  delta1_obs <- round(delta[1:5, 1], 2)
  expect_equal(delta1, delta1_obs)

  ind <- c(1, 1, 1, 1, 1)
  ind_obs <- delta[1:5, 2]
  expect_equal(ind, ind_obs)
})
