context("Variance decomposition")

test_that("Variance decomposition is correct", {
  # Spikes case
  data(ChainSC)
  
  Decomp <- BASiCS_VarianceDecomp(ChainSC, Plot = FALSE)
  
  BioVarGlobal0 <- c(0.770, 0.745, 0.666, 0.663, 0.642)
  expect_equal(round(Decomp$BioVarGlobal[1:5],3), BioVarGlobal0)
  TechVarGlobal0 <- c(0.202, 0.210, 0.247, 0.263, 0.275)
  expect_equal(round(Decomp$TechVarGlobal[1:5],3), TechVarGlobal0)
  BioVarBatch10 <- c(0.754, 0.738, 0.659, 0.627, 0.632)
  expect_equal(round(Decomp$BioVarBatch1[1:5],3), BioVarBatch10)
  TechVarBatch10 <- c(0.202, 0.218, 0.243, 0.258, 0.270)
  expect_equal(round(Decomp$TechBatch1[1:5],3), TechVarBatch10)
  
  # To emulate no-spikes case
  ChainSC1 <- ChainSC
  ChainSC1@parameters$phi <- NULL
  
  Decomp <- BASiCS_VarianceDecomp(ChainSC1, Plot = FALSE)
  
  BioVarGlobal0 <- c(0.772, 0.747, 0.675, 0.666, 0.653)
  expect_equal(round(Decomp$BioVarGlobal[1:5],3), BioVarGlobal0)
  TechVarGlobal0 <- c(0.203, 0.210, 0.246, 0.263, 0.275)
  expect_equal(round(Decomp$TechVarGlobal[1:5],3), TechVarGlobal0)
  BioVarBatch10 <- c(0.760, 0.742, 0.669, 0.636, 0.642)
  expect_equal(round(Decomp$BioVarBatch1[1:5],3), BioVarBatch10)
  TechVarBatch10 <- c(0.206, 0.219, 0.245, 0.259, 0.272)
  expect_equal(round(Decomp$TechBatch1[1:5],3), TechVarBatch10)
  Decomp <- BASiCS_VarianceDecomp(ChainSC, Plot = TRUE)
})

test_that("Variance decomposition works with 1 batch", {
  Chain <- run_MCMC(makeExampleBASiCS_Data(), 
    N = 10, Burn = 4, Thin = 2, WithSpikes = TRUE, Regression = FALSE
  )
  expect_error(BASiCS_VarianceDecomp(Chain, Plot = FALSE), NA)
  expect_error(BASiCS_VarianceDecomp(Chain, Plot = FALSE, OrderVariable = "GeneName"), NA)
  expect_error(BASiCS_VarianceDecomp(Chain, Plot = FALSE, OrderVariable = "TechVarGlobal"), NA)
  expect_error(BASiCS_VarianceDecomp(Chain, Plot = FALSE, OrderVariable = "ShotNoiseGlobal"), NA)
})
