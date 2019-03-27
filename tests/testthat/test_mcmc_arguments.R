context("MCMC arguments\n")

test_that("MCMC fails for one or multiple arguments", {
  DataSpikes <- makeExampleBASiCS_Data(WithSpikes = TRUE, 
                                       WithBatch = TRUE)
  DataSpikesNoBatch <- makeExampleBASiCS_Data(WithSpikes = TRUE, 
                                              WithBatch = FALSE)
  DataNoSpikes <- makeExampleBASiCS_Data(WithSpikes = FALSE)
  
  # Check standard MCMC setting - missing arguments
  expect_error(run_MCMC(Data = DataSpikes),
               regexp = "argument \"N\" is missing, with no default")
  expect_error(run_MCMC(Data = DataSpikes, N = 50),
               regexp = "argument \"Thin\" is missing, with no default")
  expect_error(run_MCMC(Data = DataSpikes, N = 50, 
                           Thin = 5),
               regexp = "argument \"Burn\" is missing, with no default")
  expect_error(run_MCMC(Data = DataSpikes, N = 50, 
                           Thin = 5, Burn = 25),
               regexp = "argument \"Regression\" is missing, with no default")
  
  
  DataSpikesBatchNotFactor <- DataSpikes
  SummarizedExperiment::colData(DataSpikesBatchNotFactor)$BatchInfo <- paste("Batch",
    SummarizedExperiment::colData(DataSpikesBatchNotFactor)$BatchInfo
  )
  expect_error(run_MCMC(Data = DataSpikesBatchNotFactor, N = 20, Thin = 2,
                           Burn = 4, Regression = TRUE),
               NA)
  # Check standard MCMC setting - additional parameters
  expect_error(run_MCMC(Data = DataSpikes, N = 50, 
                           Thin = 5, Burn = 25, Regression = FALSE), NA)
  expect_error(run_MCMC(Data = DataSpikes, N = 50, 
                           Thin = 5, Burn = 25, Regression = FALSE, 
                           WithSpikes = TRUE), NA)
  expect_error(run_MCMC(Data = DataNoSpikes, N = 50, 
                           Thin = 5, Burn = 25, Regression = FALSE, 
                           WithSpikes = TRUE),
               regexp = ".*does not contain TRUE values.*")
  
  Data2 <- DataSpikes
  S4Vectors::metadata(Data2)$SpikeInput <- NULL 
  expect_error(run_MCMC(Data = Data2, N = 50, 
                           Thin = 5, Burn = 25, Regression = FALSE, 
                           WithSpikes = TRUE),
               regexp = ".*does not contain the \'SpikeInput\' slot.*")
  
  Data2 <- DataSpikes
  Data2 <- clearSpikes(Data2)  
  expect_error(run_MCMC(Data = Data2, N = 50, 
                           Thin = 5, Burn = 25, Regression = FALSE, 
                           WithSpikes = TRUE),
               regexp = ".*does not contain a logical vector to indicate.*")
               
  # Check standard MCMC setting - without spikes
  expect_error(run_MCMC(Data = DataSpikes, N = 50, 
                           Thin = 5, Burn = 25, Regression = FALSE, 
                           WithSpikes = FALSE), NA)
  expect_error(run_MCMC(Data = DataSpikesNoBatch, N = 50, 
                           Thin = 5, Burn = 25, Regression = FALSE, 
                           WithSpikes = FALSE),
               regexp = ".*requires the data to contain at least 2 batches*")
  
  Data2 <- DataNoSpikes
  S4Vectors::metadata(Data2)$SpikeInput <- NULL 
  expect_error(run_MCMC(Data = Data2, N = 50, 
                           Thin = 5, Burn = 25, Regression = FALSE, 
                           WithSpikes = FALSE), NA)
  
  Data2 <- DataNoSpikes
  Data2 <- clearSpikes(Data2)  
  expect_error(run_MCMC(Data = Data2, N = 50, 
                           Thin = 5, Burn = 25, Regression = FALSE, 
                           WithSpikes = FALSE), NA)
  
  Data2 <- DataNoSpikes
  SummarizedExperiment::colData(Data2)$BatchInfo <- rep(1, length(SummarizedExperiment::colData(Data2)$BatchInfo))  
  expect_error(run_MCMC(Data = Data2, N = 50, 
                           Thin = 5, Burn = 25, Regression = FALSE, 
                           WithSpikes = FALSE),
               regexp = ".*requires the data to contain at least 2 batches*")
  
  Data2 <- DataNoSpikes
  SummarizedExperiment::colData(Data2)$BatchInfo <- NULL
  expect_error(run_MCMC(Data = Data2, N = 50, 
                           Thin = 5, Burn = 25, Regression = FALSE, 
                           WithSpikes = FALSE),
               regexp = ".*does not contain a BatchInfo vector*")
  
  # Regression BASiCS
  expect_error(run_MCMC(Data = DataSpikes, N = 50, 
                           Thin = 5, Burn = 25, Regression = TRUE), NA)
  
  expect_error(run_MCMC(Data = DataSpikes, N = 50, 
                           Thin = 5, Burn = 25, Regression = TRUE,
                           k = 2), 
               regexp = "The number of basis functions needs to be >= 4.")
  
  expect_error(run_MCMC(Data = DataSpikes, N = 50, 
                           Thin = 5, Burn = 25, Regression = TRUE, 
                           WithSpikes = TRUE), NA)
  expect_error(run_MCMC(Data = DataNoSpikes, N = 50, 
                           Thin = 5, Burn = 25, Regression = TRUE, 
                           WithSpikes = TRUE),
               regexp = ".*does not contain TRUE values.*")
  
  Data2 <- DataSpikes
  S4Vectors::metadata(Data2)$SpikeInput <- NULL 
  expect_error(run_MCMC(Data = Data2, N = 50, 
                           Thin = 5, Burn = 25, Regression = TRUE, 
                           WithSpikes = TRUE),
               regexp = ".*does not contain the \'SpikeInput\' slot.*")
  
  Data2 <- DataSpikes
  Data2 <- clearSpikes(Data2)  
  expect_error(run_MCMC(Data = Data2, N = 50, 
                           Thin = 5, Burn = 25, Regression = TRUE, 
                           WithSpikes = TRUE),
               regexp = ".*does not contain a logical vector to indicate.*")
  
  # Check standard MCMC setting - without spikes
  expect_error(run_MCMC(Data = DataSpikes, N = 50, 
                           Thin = 5, Burn = 25, Regression = TRUE, 
                           WithSpikes = FALSE), NA)
  expect_error(run_MCMC(Data = DataSpikesNoBatch, N = 50, 
                           Thin = 5, Burn = 25, Regression = TRUE, 
                           WithSpikes = FALSE),
               regexp = ".*requires the data to contain at least 2 batches*")
  
  Data2 <- DataNoSpikes
  S4Vectors::metadata(Data2)$SpikeInput <- NULL 
  expect_error(run_MCMC(Data = Data2, N = 50, 
                           Thin = 5, Burn = 25, Regression = TRUE, 
                           WithSpikes = FALSE), NA)
  
  Data2 <- DataNoSpikes
  Data2 <- clearSpikes(Data2)  
  expect_error(run_MCMC(Data = Data2, N = 50, 
                           Thin = 5, Burn = 25, Regression = TRUE, 
                           WithSpikes = FALSE), NA)
  
  Data2 <- DataNoSpikes
  SummarizedExperiment::colData(Data2)$BatchInfo <- rep(1, length(SummarizedExperiment::colData(Data2)$BatchInfo))  
  expect_error(run_MCMC(Data = Data2, N = 50, 
                           Thin = 5, Burn = 25, Regression = TRUE, 
                           WithSpikes = FALSE),
               regexp = ".*requires the data to contain at least 2 batches*")
  
  Data2 <- DataNoSpikes
  SummarizedExperiment::colData(Data2)$BatchInfo <- NULL
  expect_error(run_MCMC(Data = Data2, N = 50, 
                           Thin = 5, Burn = 25, Regression = TRUE, 
                           WithSpikes = FALSE),
               regexp = ".*does not contain a BatchInfo vector*")
  
  # Self-generated SCE object
  
  ### Without spikes
  
  # Test if it contains a 'counts' slot
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(test = counts(DataSpikes)),
    colData = SummarizedExperiment::colData(DataSpikes)
  )
  expect_error(run_MCMC(Data = sce, N = 50, 
                           Thin = 5, Burn = 25, Regression = FALSE, 
                           WithSpikes = FALSE),
               regexp = ".*does not contain a \'counts\' slot*")
  
  # Test if it is a SingleCellExperimentObject
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts(DataSpikes)),
    colData = SummarizedExperiment::colData(DataSpikes)
  )
  expect_error(run_MCMC(Data = sce, N = 50, 
                           Thin = 5, Burn = 25, Regression = FALSE, 
                           WithSpikes = FALSE),
               regexp = ".*is not a SummarizedExperiment::lass object.*")
  
  # Test if it contains a batch vector
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts(DataSpikes))
  )
  expect_error(run_MCMC(Data = sce, N = 50, 
                           Thin = 5, Burn = 25, Regression = FALSE, 
                           WithSpikes = FALSE),
               regexp = ".*does not contain a BatchInfo vector.*")

  # incorrect batch vector
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts(DataSpikes)),
    colData = data.frame(BatchInfo = SummarizedExperiment::colData(DataSpikesNoBatch)$BatchInfo)
  )
  expect_error(run_MCMC(Data = sce, N = 50, 
                           Thin = 5, Burn = 25, Regression = FALSE, 
                           WithSpikes = FALSE), 
               regexp = ".*requires the data to contain at least 2 batches.*")
    
  # Incorporate a batch vector
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts(DataSpikes)),
    colData = data.frame(BatchInfo = SummarizedExperiment::colData(DataSpikes)$BatchInfo)
  )
  expect_error(run_MCMC(Data = sce, N = 50, 
                           Thin = 5, Burn = 25, Regression = FALSE, 
                           WithSpikes = FALSE), NA)
  
  ### With spikes
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts(DataSpikes))
  )
  expect_error(run_MCMC(Data = sce, N = 50, 
                           Thin = 5, Burn = 25, Regression = FALSE, 
                           WithSpikes = TRUE),
               regexp = ".*does not contain a logical vector to indicate.*")
  
  isSpike(sce, "ERCC") <- isSpike(DataSpikes)
  expect_error(run_MCMC(Data = sce, N = 50, 
                           Thin = 5, Burn = 25, Regression = FALSE, 
                           WithSpikes = TRUE),
               regexp = ".*does not contain the \'SpikeInput\' slot.*")
  
  # Wrong SpikeInput
  S4Vectors::metadata(sce)$SpikeInput <- 1
  expect_error(run_MCMC(Data = sce, N = 50, 
                           Thin = 5, Burn = 25, Regression = FALSE, 
                           WithSpikes = TRUE),
               regexp = ".*Spike-in assignment is not compatible with data.*")
  
  S4Vectors::metadata(sce)$SpikeInput <- c(1,2,3)
  expect_error(run_MCMC(Data = sce, N = 50, 
                           Thin = 5, Burn = 25, Regression = FALSE, 
                           WithSpikes = TRUE),
               regexp = ".*Spike-in assignment is not compatible with data.*")
  
  # Right SpikeInfo assignment  
  S4Vectors::metadata(sce)$SpikeInput <- S4Vectors::metadata(DataSpikes)$SpikeInput
  expect_error(run_MCMC(Data = sce, N = 50, 
                           Thin = 5, Burn = 25, Regression = FALSE, 
                           WithSpikes = TRUE), NA)

  
})



