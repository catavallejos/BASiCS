context("MCMC arguments\n")

test_that("MCMC fails for one or multiple arguments", {
  DataSpikes <- makeExampleBASiCS_Data(WithSpikes = TRUE, 
                                       WithBatch = TRUE)
  DataSpikesNoBatch <- makeExampleBASiCS_Data(WithSpikes = TRUE, 
                                              WithBatch = FALSE)
  DataNoSpikes <- makeExampleBASiCS_Data(WithSpikes = FALSE)
  
  # Check standard MCMC setting - missing arguments
  expect_error(BASiCS_MCMC(Data = DataSpikes),
               regexp = "argument \"N\" is missing, with no default")
  expect_error(BASiCS_MCMC(Data = DataSpikes, N = 50),
               regexp = "argument \"Thin\" is missing, with no default")
  expect_error(BASiCS_MCMC(Data = DataSpikes, N = 50, 
                           Thin = 5),
               regexp = "argument \"Burn\" is missing, with no default")
  expect_error(BASiCS_MCMC(Data = DataSpikes, N = 50, 
                           Thin = 5, Burn = 25),
               regexp = "argument \"Regression\" is missing, with no default")
  
  # Check standard MCMC setting - additional parameters
  expect_error(BASiCS_MCMC(Data = DataSpikes, N = 50, 
                           Thin = 5, Burn = 25, Regression = FALSE), NA)
  expect_error(BASiCS_MCMC(Data = DataSpikes, N = 50, 
                           Thin = 5, Burn = 25, Regression = FALSE, 
                           WithSpikes = TRUE), NA)
  expect_error(BASiCS_MCMC(Data = DataNoSpikes, N = 50, 
                           Thin = 5, Burn = 25, Regression = FALSE, 
                           WithSpikes = TRUE),
               regexp = ".*does not contain TRUE values.*")
  
  Data2 <- DataSpikes
  metadata(Data2)$SpikeInput <- NULL 
  expect_error(BASiCS_MCMC(Data = Data2, N = 50, 
                           Thin = 5, Burn = 25, Regression = FALSE, 
                           WithSpikes = TRUE),
               regexp = ".*does not contain the \'SpikeInput\' slot.*")
  
  Data2 <- DataSpikes
  Data2 <- clearSpikes(Data2)  
  expect_error(BASiCS_MCMC(Data = Data2, N = 50, 
                           Thin = 5, Burn = 25, Regression = FALSE, 
                           WithSpikes = TRUE),
               regexp = ".*does not contain a logical vector to indicate.*")
               
  # Check standard MCMC setting - without spikes
  expect_error(BASiCS_MCMC(Data = DataSpikes, N = 50, 
                           Thin = 5, Burn = 25, Regression = FALSE, 
                           WithSpikes = FALSE), NA)
  expect_error(BASiCS_MCMC(Data = DataSpikesNoBatch, N = 50, 
                           Thin = 5, Burn = 25, Regression = FALSE, 
                           WithSpikes = FALSE),
               regexp = ".*requires the data to contain at least 2 batches*")
  
  Data2 <- DataNoSpikes
  metadata(Data2)$SpikeInput <- NULL 
  expect_error(BASiCS_MCMC(Data = Data2, N = 50, 
                           Thin = 5, Burn = 25, Regression = FALSE, 
                           WithSpikes = FALSE), NA)
  
  Data2 <- DataNoSpikes
  Data2 <- clearSpikes(Data2)  
  expect_error(BASiCS_MCMC(Data = Data2, N = 50, 
                           Thin = 5, Burn = 25, Regression = FALSE, 
                           WithSpikes = FALSE), NA)
  
  Data2 <- DataNoSpikes
  colData(Data2)$BatchInfo <- rep(1, length(colData(Data2)$BatchInfo))  
  expect_error(BASiCS_MCMC(Data = Data2, N = 50, 
                           Thin = 5, Burn = 25, Regression = FALSE, 
                           WithSpikes = FALSE),
               regexp = ".*requires the data to contain at least 2 batches*")
  
  Data2 <- DataNoSpikes
  colData(Data2)$BatchInfo <- NULL
  expect_error(BASiCS_MCMC(Data = Data2, N = 50, 
                           Thin = 5, Burn = 25, Regression = FALSE, 
                           WithSpikes = FALSE),
               regexp = ".*does not contain a BatchInfo vector*")
  
  # Regression BASiCS
  expect_error(BASiCS_MCMC(Data = DataSpikes, N = 50, 
                           Thin = 5, Burn = 25, Regression = TRUE), NA)
  
  expect_error(BASiCS_MCMC(Data = DataSpikes, N = 50, 
                           Thin = 5, Burn = 25, Regression = TRUE,
                           k = 2), 
               regexp = "The number of basis functions needs to be >= 4.")
  
  expect_error(BASiCS_MCMC(Data = DataSpikes, N = 50, 
                           Thin = 5, Burn = 25, Regression = TRUE, 
                           WithSpikes = TRUE), NA)
  expect_error(BASiCS_MCMC(Data = DataNoSpikes, N = 50, 
                           Thin = 5, Burn = 25, Regression = TRUE, 
                           WithSpikes = TRUE),
               regexp = ".*does not contain TRUE values.*")
  
  Data2 <- DataSpikes
  metadata(Data2)$SpikeInput <- NULL 
  expect_error(BASiCS_MCMC(Data = Data2, N = 50, 
                           Thin = 5, Burn = 25, Regression = TRUE, 
                           WithSpikes = TRUE),
               regexp = ".*does not contain the \'SpikeInput\' slot.*")
  
  Data2 <- DataSpikes
  Data2 <- clearSpikes(Data2)  
  expect_error(BASiCS_MCMC(Data = Data2, N = 50, 
                           Thin = 5, Burn = 25, Regression = TRUE, 
                           WithSpikes = TRUE),
               regexp = ".*does not contain a logical vector to indicate.*")
  
  # Check standard MCMC setting - without spikes
  expect_error(BASiCS_MCMC(Data = DataSpikes, N = 50, 
                           Thin = 5, Burn = 25, Regression = TRUE, 
                           WithSpikes = FALSE), NA)
  expect_error(BASiCS_MCMC(Data = DataSpikesNoBatch, N = 50, 
                           Thin = 5, Burn = 25, Regression = TRUE, 
                           WithSpikes = FALSE),
               regexp = ".*requires the data to contain at least 2 batches*")
  
  Data2 <- DataNoSpikes
  metadata(Data2)$SpikeInput <- NULL 
  expect_error(BASiCS_MCMC(Data = Data2, N = 50, 
                           Thin = 5, Burn = 25, Regression = TRUE, 
                           WithSpikes = FALSE), NA)
  
  Data2 <- DataNoSpikes
  Data2 <- clearSpikes(Data2)  
  expect_error(BASiCS_MCMC(Data = Data2, N = 50, 
                           Thin = 5, Burn = 25, Regression = TRUE, 
                           WithSpikes = FALSE), NA)
  
  Data2 <- DataNoSpikes
  colData(Data2)$BatchInfo <- rep(1, length(colData(Data2)$BatchInfo))  
  expect_error(BASiCS_MCMC(Data = Data2, N = 50, 
                           Thin = 5, Burn = 25, Regression = TRUE, 
                           WithSpikes = FALSE),
               regexp = ".*requires the data to contain at least 2 batches*")
  
  Data2 <- DataNoSpikes
  colData(Data2)$BatchInfo <- NULL
  expect_error(BASiCS_MCMC(Data = Data2, N = 50, 
                           Thin = 5, Burn = 25, Regression = TRUE, 
                           WithSpikes = FALSE),
               regexp = ".*does not contain a BatchInfo vector*")
  
  # Self-generated SCE object
  
  ### Without spikes
  
  # Test if it contains a 'counts' slot
  sce <- SingleCellExperiment(assays = list(test = counts(DataSpikes)),
                               colData = colData(DataSpikes))
  expect_error(BASiCS_MCMC(Data = sce, N = 50, 
                           Thin = 5, Burn = 25, Regression = FALSE, 
                           WithSpikes = FALSE),
               regexp = ".*does not contain a \'counts\' slot*")
  
  # Test if it is a SingleCellExperimentObject
  sce <- SummarizedExperiment(assays = list(counts = counts(DataSpikes)),
                               colData = colData(DataSpikes))
  expect_error(BASiCS_MCMC(Data = sce, N = 50, 
                           Thin = 5, Burn = 25, Regression = FALSE, 
                           WithSpikes = FALSE),
               regexp = ".*is not a SingleCellExperiment class object.*")
  
  # Test if it contains a batch vector
  sce <- SingleCellExperiment(assays = list(counts = counts(DataSpikes)))
  expect_error(BASiCS_MCMC(Data = sce, N = 50, 
                           Thin = 5, Burn = 25, Regression = FALSE, 
                           WithSpikes = FALSE),
               regexp = ".*does not contain a BatchInfo vector.*")

  # incorrect batch vector
  sce <- SingleCellExperiment(assays = list(counts = counts(DataSpikes)),
              colData = data.frame(BatchInfo = colData(DataSpikesNoBatch)$BatchInfo))
  expect_error(BASiCS_MCMC(Data = sce, N = 50, 
                           Thin = 5, Burn = 25, Regression = FALSE, 
                           WithSpikes = FALSE), 
               regexp = ".*requires the data to contain at least 2 batches.*")
    
  # Incorporate a batch vector
  sce <- SingleCellExperiment(assays = list(counts = counts(DataSpikes)),
              colData = data.frame(BatchInfo = colData(DataSpikes)$BatchInfo))
  expect_error(BASiCS_MCMC(Data = sce, N = 50, 
                           Thin = 5, Burn = 25, Regression = FALSE, 
                           WithSpikes = FALSE), NA)
  
  ### With spikes
  sce <- SingleCellExperiment(assays = list(counts = counts(DataSpikes)))
  expect_error(BASiCS_MCMC(Data = sce, N = 50, 
                           Thin = 5, Burn = 25, Regression = FALSE, 
                           WithSpikes = TRUE),
               regexp = ".*does not contain a logical vector to indicate.*")
  
  isSpike(sce, "ERCC") <- isSpike(DataSpikes)
  expect_error(BASiCS_MCMC(Data = sce, N = 50, 
                           Thin = 5, Burn = 25, Regression = FALSE, 
                           WithSpikes = TRUE),
               regexp = ".*does not contain the \'SpikeInput\' slot.*")
  
  # Wrong SpikeInput
  metadata(sce)$SpikeInput <- 1
  expect_error(BASiCS_MCMC(Data = sce, N = 50, 
                           Thin = 5, Burn = 25, Regression = FALSE, 
                           WithSpikes = TRUE),
               regexp = ".*Spike-in assignment is not compatible with data.*")
  
  metadata(sce)$SpikeInput <- c(1,2,3)
  expect_error(BASiCS_MCMC(Data = sce, N = 50, 
                           Thin = 5, Burn = 25, Regression = FALSE, 
                           WithSpikes = TRUE),
               regexp = ".*Spike-in assignment is not compatible with data.*")
  
  # Right SpikeInfo assignment  
  metadata(sce)$SpikeInput <- metadata(DataSpikes)$SpikeInput
  expect_error(BASiCS_MCMC(Data = sce, N = 50, 
                           Thin = 5, Burn = 25, Regression = FALSE, 
                           WithSpikes = TRUE))
})



