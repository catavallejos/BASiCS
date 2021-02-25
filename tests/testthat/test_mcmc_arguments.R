context("MCMC arguments")

set.seed(1)
DataSpikes <- makeExampleBASiCS_Data(WithSpikes = TRUE, WithBatch = TRUE)
DataSpikesNoBatch <- makeExampleBASiCS_Data(WithSpikes = TRUE, WithBatch = FALSE)
DataNoSpikes <- makeExampleBASiCS_Data(WithSpikes = FALSE)

test_that("Errors in basic MCMC arguments", {
  
  expect_error(run_MCMC(Data = DataSpikes),
               regexp = "argument \"N\" is missing, with no default")
  expect_error(run_MCMC(Data = DataSpikes, N = 10),
               regexp = "argument \"Thin\" is missing, with no default")
  expect_error(run_MCMC(Data = DataSpikes, N = 10, Thin = 2),
               regexp = "argument \"Burn\" is missing, with no default")
  expect_error(run_MCMC(Data = DataSpikes, N = 10, Thin = 2, Burn = 4),
               regexp = "argument \"Regression\" is missing, with no default")
  pp <- BASiCS_PriorParam(DataSpikes, MinGenesPerRBF = NA)
  expect_error(run_MCMC(Data = DataSpikes, 
                        N = 10, Thin = 2, Burn = 4, 
                        Regression = TRUE, PriorParam = pp),
               NA)
  expect_error(run_MCMC(Data = DataSpikes, 
                        N = 10, Thin = 2, Burn = 4, 
                        Regression = TRUE, WithSpikes = TRUE, PriorParam = pp),
               NA)
  # Check batch info when it contains text
  DataSpikesBatchNotFactor <- DataSpikes
  colData(DataSpikesBatchNotFactor)$BatchInfo <- 
    paste("Batch", colData(DataSpikesBatchNotFactor)$BatchInfo)
  expect_error(run_MCMC(Data = DataSpikesBatchNotFactor, 
                        N = 20, Thin = 2, Burn = 4, Regression = TRUE), NA)
})

test_that("MCMC arguments fail (spikes; no-regression)", {
  
  # Check standard MCMC setting - additional parameters
  expect_error(run_MCMC(Data = DataSpikes, 
                        N = 10, Thin = 2, Burn = 4, 
                        Regression = FALSE), NA)
  pp <- BASiCS_PriorParam(DataSpikes, PriorDelta="gamma")
  expect_error(run_MCMC(Data = DataSpikes, 
                        N = 10, Thin = 2, Burn = 4,
                        PriorParam = pp,
                        Regression = FALSE), NA)

  expect_error(run_MCMC(Data = DataSpikes, 
                        N = 10, Thin = 2, Burn = 4, 
                        Regression = FALSE, WithSpikes = TRUE), NA)
  expect_error(run_MCMC(Data = DataNoSpikes, 
                        N = 10, Thin = 2, Burn = 4, 
                        Regression = FALSE, WithSpikes = TRUE),
               regexp = ".*does not contain information about spike-in genes*")
  
  # Checks after removing spike-ins metadata but not data
  Data2 <- DataSpikes
  rowData(altExp(Data2)) <- NULL
  expect_error(run_MCMC(Data = Data2, 
                        N = 20, Thin = 2, Burn = 4, 
                        Regression = FALSE, WithSpikes = TRUE),
               regexp = "rowData of altExp is missing.*")
  # Checks after removing spike-ins data but not metadata
  Data2 <- DataSpikes
  altExp(Data2) <- NULL
  expect_error(run_MCMC(Data = Data2, 
                        N = 20, Thin = 2, Burn = 4, 
                        Regression = FALSE, WithSpikes = TRUE),
               regexp = ".*does not contain information about spike-in genes*")
               
  # Check if MCMC is to ignore spike-ins
  expect_error(run_MCMC(Data = DataSpikes, 
                        N = 20, Thin = 2, Burn = 4, 
                        Regression = FALSE, WithSpikes = FALSE), NA)
  expect_warning(run_MCMC(Data = DataSpikesNoBatch, 
                        N = 20, Thin = 2, Burn = 4, 
                        Regression = FALSE, WithSpikes = FALSE),
               regexp = ".*recommends that the data contain at least 2 batches*")
  Data2 <- DataNoSpikes
  altExp(Data2) <- NULL  
  expect_error(run_MCMC(Data = Data2, 
                        N = 10, Thin = 2, Burn = 4, 
                        Regression = FALSE, WithSpikes = FALSE), NA)
  # Checks if data does not have spike-ins
  expect_error(run_MCMC(Data = DataNoSpikes, 
                         N = 10, Thin = 2, Burn = 4, 
                        Regression = TRUE, WithSpikes = TRUE),
               regexp = ".*does not contain information about spike-in genes.*")
})

test_that("MCMC arguments fail (no-spikes; no-regression)", {  
  
  # Checks if only one batch
  Data2 <- DataNoSpikes
  colData(Data2)$BatchInfo <- rep(1, length(colData(Data2)$BatchInfo))  
  expect_warning(run_MCMC(Data = Data2, 
                        N = 10, Thin = 2, Burn = 4, 
                        Regression = FALSE, WithSpikes = FALSE),
               regexp = ".*recommends that the data contain at least 2 batches*")
  # Checks if only batch info is not available
  Data2 <- DataNoSpikes
  SummarizedExperiment::colData(Data2)$BatchInfo <- NULL
  expect_warning(run_MCMC(Data = Data2, 
                        N = 10, Thin = 2, Burn = 4, 
                        Regression = FALSE, WithSpikes = FALSE),
               regexp = ".*should contain a BatchInfo vector*")

})

test_that("MCMC arguments fail (regression)", {   

  expect_error(run_MCMC(Data = DataSpikes, 
                        N = 10, Thin = 2, Burn = 4, 
                        PriorParam = BASiCS_PriorParam(DataSpikes, k = 2),
                        Regression = TRUE), 
               regexp = "The number of basis functions needs to be >= 4.")
  
  Data2 <- DataSpikes
  rowData(altExp(Data2)) <- NULL
  expect_error(run_MCMC(Data = Data2, 
                        N = 10, Thin = 2, Burn = 4, 
                        Regression = TRUE, WithSpikes = TRUE),
               regexp = ".*rowData of altExp is missing!.*")
  
  Data2 <- DataSpikes
  altExp(Data2) <- NULL   
  expect_error(run_MCMC(Data = Data2, 
                        N = 10, Thin = 2, Burn = 4, 
                        Regression = TRUE, WithSpikes = TRUE),
               regexp = ".*does not contain information about spike-in genes*")
  
  # Check standard MCMC setting - without spikes
  expect_error(run_MCMC(Data = DataSpikes, 
                        N = 10, Thin = 2, Burn = 4, 
                        Regression = TRUE, WithSpikes = FALSE), NA)
  expect_warning(run_MCMC(Data = DataSpikesNoBatch, 
                        N = 10, Thin = 2, Burn = 4, 
                        Regression = TRUE, WithSpikes = FALSE),
               regexp = ".*recommends that the data contain at least 2 batches*")
  
  Data2 <- DataNoSpikes
  expect_error(run_MCMC(Data = Data2, 
                        N = 10, Thin = 2, Burn = 4, 
                        Regression = TRUE, WithSpikes = FALSE), NA)
  
  Data2 <- DataNoSpikes
  altExp(Data2) <- NULL 
  expect_error(run_MCMC(Data = Data2, 
                        N = 10, Thin = 2, Burn = 4, 
                        Regression = TRUE, WithSpikes = FALSE), NA)
  
  Data2 <- DataNoSpikes
  colData(Data2)$BatchInfo <- rep(1, length(colData(Data2)$BatchInfo))  
  expect_warning(run_MCMC(Data = Data2, 
                        N = 10, Thin = 2, Burn = 4, 
                        Regression = TRUE, WithSpikes = FALSE),
               regexp = ".*recommends that the data contain at least 2 batches*")
  
  Data2 <- DataNoSpikes
  colData(Data2)$BatchInfo <- NULL
  expect_warning(run_MCMC(Data = Data2, 
                        N = 10, Thin = 2, Burn = 4, 
                        Regression = TRUE, WithSpikes = FALSE),
               regexp = ".*should contain a BatchInfo vector*")

  expect_error(
    run_MCMC(DataNoSpikes, N = 10, Thin = 2, Burn = 4, 
             PrintProgress = FALSE, Regression = TRUE,
             WithSpikes = FALSE,
             PriorParam = BASiCS_PriorParam(DataNoSpikes, k = 2, MinGenesPerRBF = 100)),
    "Consider setting MinGenesPerRBF to NA or a lower positive integer."
  )
})

test_that("Checks for user generated SCE object", {    
  
  # Test if it contains a 'counts' slot
  sce <- SingleCellExperiment(assays = list(test = counts(DataSpikes)),
                              colData = colData(DataSpikes))
  expect_error(run_MCMC(Data = sce, 
                         N = 20, Thin = 2, Burn = 4, 
                        Regression = FALSE, WithSpikes = FALSE),
               regexp = ".*does not contain a \'counts\' slot*")
  
  # Test if it is a SingleCellExperimentObject
  sce <- SummarizedExperiment(assays = list(counts = counts(DataSpikes)),
                              colData = colData(DataSpikes))
  expect_error(run_MCMC(Data = sce, 
                         N = 20, Thin = 2, Burn = 4, 
                        Regression = FALSE, WithSpikes = FALSE),
               regexp = ".*is not a SingleCellExperiment class object.*")
  
  # Test if it contains a batch vector
  sce <- SingleCellExperiment(assays = list(counts = counts(DataSpikes)))
  expect_warning(run_MCMC(Data = sce, 
                         N = 20, Thin = 2, Burn = 4, 
                        Regression = FALSE, WithSpikes = FALSE),
               regexp = "should contain a BatchInfo vector")

  # incorrect batch vector
  sce <- SingleCellExperiment(
    assays = list(counts = counts(DataSpikes)),
    colData = data.frame(BatchInfo = colData(DataSpikesNoBatch)$BatchInfo)
  )
  expect_warning(run_MCMC(Data = sce, 
                         N = 20, Thin = 2, Burn = 4, 
                        Regression = FALSE, WithSpikes = FALSE), 
               regexp = ".*recommends that the data contain at least 2 batches.*")
    
  # Incorporate a batch vector
  sce <- SingleCellExperiment(
    assays = list(counts = counts(DataSpikes)),
    colData = data.frame(BatchInfo = colData(DataSpikes)$BatchInfo)
  )
  expect_error(run_MCMC(Data = sce, 
                         N = 20, Thin = 2, Burn = 4, 
                        Regression = FALSE, WithSpikes = FALSE), NA)
  
  ### With spikes
  sce <- SingleCellExperiment(assays = list(counts = counts(DataSpikes)))
  expect_error(run_MCMC(Data = sce, 
                         N = 20, Thin = 2, Burn = 4, 
                        Regression = FALSE, WithSpikes = TRUE),
               regexp = ".*does not contain information about spike-in genes*")
  
  altExp(sce, "spike-ins") <- altExp(DataSpikes)
  rowData(altExp(sce)) <- NULL
  
  # Wrong SpikeInput
  rowData(altExp(sce)) <- NULL
  expect_error(run_MCMC(Data = sce, 
                         N = 20, Thin = 2, Burn = 4, 
                        Regression = FALSE, WithSpikes = TRUE),
               regexp = ".*rowData of altExp is missing.*")
  
  rowData(altExp(sce))[[1]] <- 1:2
  expect_error(run_MCMC(Data = sce, 
                         N = 20, Thin = 2, Burn = 4, 
                        Regression = FALSE, WithSpikes = TRUE),
               regexp = ".*rowData of altExp must have two columns only.*")
  
  # Right SpikeInfo assignment  
  rowData(altExp(sce)) <- rowData(altExp(DataSpikes))
  expect_error(
    run_MCMC(
      Data = sce,
      N = 10,
      Thin = 2,
      Burn = 4,
      Regression = FALSE,
      WithSpikes = TRUE
    ),
    NA
  )
})


test_that("MCMC works with different input classes", {
  counts(DataSpikes) <- Matrix::Matrix(counts(DataSpikes))
  expect_error(
    run_MCMC(
      Data = DataSpikes,
      N = 10,
      Thin = 2,
      Burn = 4, 
      Regression = TRUE
    ),
    NA
  )
  counts(DataSpikes) <- Matrix::Matrix(counts(DataSpikes), sparse = TRUE)
  expect_error(
    run_MCMC(
      Data = DataSpikes, 
      N = 10,
      Thin = 2,
      Burn = 4, 
      Regression = TRUE
    ),
    NA
  )
})

test_that("PriorMu", {
  set.seed(1)
  expect_error(
    BASiCS_PriorParam(DataSpikes, PriorMu = "das"),
    "'arg' should be one of"
  )
  expect_error(
    BASiCS_PriorParam(DataSpikes, PriorMu = "default"),
    NA
  )
  expect_error(
    BASiCS_PriorParam(DataSpikes, PriorMu = "EmpiricalBayes"),
    NA
  )
})
