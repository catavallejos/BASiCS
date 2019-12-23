context("Individual MCMC updates (cpp code)\n")

test_that("Dirichlet sampler", {
  # Generating arbitrary hyper-param
  phi0 <- 1:10; phi0 <- phi0/sum(phi0)
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
  q0 <- nrow(CountsBio); n <- ncol(CountsBio)
  PriorParam <- list(mu.mu = 0, s2.mu = 0.5, s2.delta = 0.5, a.delta = 1, 
                     b.delta = 1, p.phi = rep(1, times = n), 
                     a.s = 1, b.s = 1, a.theta = 1, b.theta = 1)
  set.seed(2018)
  Start <- BASiCS:::HiddenBASiCS_MCMC_Start(Data, 
                                            WithSpikes = TRUE)
  uGene <- rep(0, times = q0)
  indGene <- rbinom(q0, size = 1, prob = 0.5)
  
  # Hidden_muUpdate
  mu1 <- pmax(0, Start$mu0[seq_len(q0)] + rnorm(q0, sd = 0.005))
  Aux <- BASiCS:::Hidden_muUpdate(mu0 = Start$mu0, 
                                  prop_var = exp(Start$ls.mu0),

                                  Counts = CountsBio, 
                                  invdelta = 1/Start$delta0, 
                                  phinu = Start$phi0 * Start$nu0,
                                  sum_bycell_bio = rowSums(CountsBio),
                                  mu_mu = PriorParam$mu.mu,
                                  s2_mu = PriorParam$s2.mu, 
                                  q0 = q0, n = n,
                                  mu1 =  mu1, 
                                  u = uGene, 
                                  ind = indGene,
                                  exponent = 1,
                                  mintol = 1e-3)


  mu1 <- c(6.78, 15.67,  5.43, 13.05, 24.20)
  mu1Obs <- round(Aux[1:5,1],2) 
  expect_equal(mu1, mu1Obs)

  ind <- c(1, 0, 1, 1, 1)
  indObs <- Aux[1:5,2]
  expect_equal(ind, indObs)
})
