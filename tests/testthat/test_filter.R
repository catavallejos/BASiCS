context("BASiCS_Filter\n")

test_that("BASiCS_Filter", 
{

  set.seed(1)
  Counts <- matrix(rpois(50*10, 2), ncol = 10)
  rownames(Counts) <- c(paste0('Gene', 1:40), paste0('Spike', 1:10))
  # Two genes with zero total counts
  Counts[3, ] <- 0; Counts[42, ] <- 0
  # Once cell with zero total counts
  Counts[, 7] <- 0
  Tech <- c(rep(FALSE,40),rep(TRUE,10))
  set.seed(2)
  SpikeInput <- rgamma(10,1,1)
  SpikeInfo <- data.frame('SpikeID' = paste0('Spike', 1:10), 
                          'SpikeInput' = SpikeInput)
  Filter <- BASiCS_Filter(Counts, Tech, SpikeInput,
                          MinTotalCountsPerCell = 2, 
                          MinTotalCountsPerGene = 2,
                          MinCellsWithExpression = 2, 
                          MinAvCountsPerCellsWithExpression = 2)
  expect_true(all.equal(names(Filter), 
                        c("Counts", "Tech", "SpikeInput", "BatchInfo", 
                          "IncludeGenes", "IncludeCells")))
  
  IncludeCells <- rep(TRUE, times = 10); IncludeCells[7] <- FALSE
  expect_true(all.equal(Filter$IncludeCells, IncludeCells))
  
  IncludeGenes <- rep(TRUE, times = 50); 
  IncludeGenes[c(3, 12, 24, 25, 32, 36, 41, 42, 45, 47)] <- FALSE
  expect_true(all.equal(Filter$IncludeGenes, IncludeGenes))
  
  expect_true(all.equal(rownames(Filter$Counts)[Filter$Tech],
                        c("Spike3", "Spike4", "Spike6", 
                          "Spike8", "Spike9", "Spike10")))
})

