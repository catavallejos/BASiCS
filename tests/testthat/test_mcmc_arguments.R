context("MCMC arguments\n")

test_that("MCMC fails for one or multiple arguments", {
  DataSpikes <- makeExampleBASiCS_Data(WithSpikes = TRUE, 
                                       WithBatch = TRUE)
  DataSpikesNoBatch <- makeExampleBASiCS_Data(WithSpikes = TRUE, 
                                              WithBatch = FALSE)
  DataNoSpikes <- makeExampleBASiCS_Data(WithSpikes = FALSE)
  
  # Check standard MCMC setting - missing arguments
  expect_error(BASiCS_MCMC(Data = DataSpikes))
  expect_error(BASiCS_MCMC(Data = DataSpikes, N = 50))
  expect_error(BASiCS_MCMC(Data = DataSpikes, N = 50, 
                           Thin = 5))
  expect_error(BASiCS_MCMC(Data = DataSpikes, N = 50, 
                           Thin = 5, Burn = 25))
  
  # Check standard MCMC setting - additional parameters
  expect_error(BASiCS_MCMC(Data = DataSpikes, N = 50, 
                           Thin = 5, Burn = 25, Regression = FALSE), NA)
  expect_error(BASiCS_MCMC(Data = DataSpikes, N = 50, 
                           Thin = 5, Burn = 25, Regression = FALSE, 
                           WithSpikes = TRUE), NA)
  expect_error(BASiCS_MCMC(Data = DataNoSpikes, N = 50, 
                           Thin = 5, Burn = 25, Regression = FALSE, 
                           WithSpikes = TRUE))
  Data2 <- DataSpikes
  metadata(Data2)$SpikeInput <- NULL 
  expect_error(BASiCS_MCMC(Data = Data2, N = 50, 
                           Thin = 5, Burn = 25, Regression = FALSE, 
                           WithSpikes = TRUE))
  
  Data2 <- DataSpikes
  Data2 <- clearSpikes(Data2)  
  expect_error(BASiCS_MCMC(Data = Data2, N = 50, 
                           Thin = 5, Burn = 25, Regression = FALSE, 
                           WithSpikes = TRUE))
               
  # Check standard MCMC setting - without spikes
  expect_error(BASiCS_MCMC(Data = DataSpikes, N = 50, 
                           Thin = 5, Burn = 25, Regression = FALSE, 
                           WithSpikes = FALSE), NA)
  expect_error(BASiCS_MCMC(Data = DataSpikesNoBatch, N = 50, 
                           Thin = 5, Burn = 25, Regression = FALSE, 
                           WithSpikes = FALSE))
  
  Data2 <- DataNoSpikes
  metadata(Data2)$SpikeInput <- NULL 
  expect_error(BASiCS_MCMC(Data = Data2, N = 50, 
                           Thin = 5, Burn = 25, Regression = FALSE, 
                           WithSpikes = FALSE), NA)
  
  Data2 <- DataNoSpikes
  Data2 <- clearSpikes(Data2)  
  expect_error(BASiCS_MCMC(Data = Data2, N = 50, 
                           Thin = 5, Burn = 25, Regression = FALSE, 
                           WithSpikes = FALSE))
  
  # Regression BASiCS
  
  # Self-generated SCE object
})



