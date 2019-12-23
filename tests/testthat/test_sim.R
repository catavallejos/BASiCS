context("simulated data")

test_that("BASiCS_Sim works", {
  # Simulated parameter values for 10 genes
  # (7 biogical and 3 spike-in) measured in 5 cells
  Mu <- c(8.36, 10.65, 4.88, 6.29, 21.72, 12.93, 30.19)
  Mu_spikes <-  c(1010.72, 7.90, 31.59)
  Delta <- c(1.29, 0.88, 1.51, 1.49, 0.54, 0.40, 0.85)
  Phi <- c(1.00, 1.06, 1.09, 1.05, 0.80)
  S <- c(0.38, 0.40, 0.38, 0.39, 0.34)
  Theta <- 0.39
  
  # Data with spike-ins, single batch
  set.seed(1)
  Data <- BASiCS_Sim(Mu, Mu_spikes, Delta, Phi, S, Theta)
  
  # Check if values are reproducible given fixed seed
  Aux <- as.vector(SingleCellExperiment::counts(Data)[1:5, 1])
  Aux0 <- c(6, 0, 0, 1, 11)
  expect_equal(Aux, Aux0)
  Aux <- sum(SingleCellExperiment::counts(Data))
  
  Aux0 <- 201
  expect_equal(Aux, Aux0)
  
  # Data with spike-ins, multiple batches
  BatchInfo <- c(rep(1, 3), rep(2, 2))
  Theta2 <- rep(Theta, times = 2)
  set.seed(2)
  Data <- BASiCS_Sim(Mu, Mu_spikes, Delta, Phi, S, Theta2, BatchInfo)
  expect_is(Data, "SingleCellExperiment")
  # Check if values are reproducible given fixed seed
  Aux <- as.vector(SingleCellExperiment::counts(Data)[1:5, 1])
  Aux0 <- c(0, 2, 2, 0, 0)
  expect_equal(Aux, Aux0)
  Aux <- sum(SingleCellExperiment::counts(Data))
  Aux0 <- 74
  expect_equal(Aux, Aux0)
  
  # Data without spike-ins, multiple batches
  set.seed(3)
  Data <- BASiCS_Sim(Mu, Mu_spikes = NULL, Delta, Phi = NULL, S, 
                     Theta2, BatchInfo)
  expect_is(Data, "SingleCellExperiment")
  # Check if values are reproducible given fixed seed
  Aux <- as.vector(SingleCellExperiment::counts(Data)[1:5, 1])
  Aux0 <- c(0, 1, 0, 0, 4)
  expect_equal(Aux, Aux0)
  Aux <- sum(SingleCellExperiment::counts(Data))
  Aux0 <- 80
  expect_equal(Aux, Aux0)
  # When the parameter input is not right
  expect_error(
    BASiCS_Sim(Mu, Mu_spikes = NULL, Delta, Phi = NULL, S, Theta),
    "When spike-ins are not included, 'BatchInfo' is required."
  )
})

test_that("BASiCS_Draw works", {
  data <- makeExampleBASiCS_Data(WithBatch = TRUE)
  chain <- run_MCMC(data,
                    N = 10,
                    Thin = 2,
                    Burn = 4,
                    Regression = FALSE,
                    WithSpikes = FALSE
  )
  draw <- BASiCS_Draw(chain)
  expect_is(draw, "SingleCellExperiment")
})