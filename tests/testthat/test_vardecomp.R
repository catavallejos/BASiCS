context("Variance decomposition\n")

test_that("Variance decomposition is correct", 
{
  # Spikes case
  data(ChainSC)
  
  VD <- BASiCS_VarianceDecomp(ChainSC, Plot = FALSE)
  
  BioVarGlobal0 <- c(0.771, 0.746, 0.675, 0.670, 0.665)
  expect_that(all.equal(round(VD$BioVarGlobal[1:5],3), 
                        BioVarGlobal0), is_true())
  TechVarGlobal0 <- c(0.203, 0.210, 0.263, 0.250, 0.263)
  expect_that(all.equal(round(VD$TechVarGlobal[1:5],3), 
                        TechVarGlobal0), is_true())
  BioVarBatch10 <- c(0.755, 0.739, 0.653, 0.662, 0.630)
  expect_that(all.equal(round(VD$BioVarBatch1[1:5],3), 
                        BioVarBatch10), is_true())
  TechVarBatch10 <- c(0.205, 0.219, 0.257, 0.244, 0.258)
  expect_that(all.equal(round(VD$TechBatch1[1:5],3), 
                        TechVarBatch10), is_true()) 
  
  # To emulate no-spikes case
  ChainSC1 <- ChainSC
  ChainSC1@parameters$phi <- NULL
  
  VD <- BASiCS_VarianceDecomp(ChainSC1, Plot = FALSE)
  
  BioVarGlobal0 <- c(0.772, 0.747, 0.675, 0.667, 0.666)
  expect_that(all.equal(round(VD$BioVarGlobal[1:5],3), 
                        BioVarGlobal0), is_true())
  TechVarGlobal0 <- c(0.203, 0.210, 0.246, 0.264, 0.263)
  expect_that(all.equal(round(VD$TechVarGlobal[1:5],3), 
                        TechVarGlobal0), is_true())
  BioVarBatch10 <- c(0.760, 0.742, 0.669, 0.663, 0.636)
  expect_that(all.equal(round(VD$BioVarBatch1[1:5],3), 
                        BioVarBatch10), is_true())
  TechVarBatch10 <- c(0.206, 0.219, 0.245, 0.261, 0.259)
  expect_that(all.equal(round(VD$TechBatch1[1:5],3), 
                        TechVarBatch10), is_true()) 
})




