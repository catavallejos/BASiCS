HiddenChecksBASiCS_Data <- function(Data,
                                    WithSpikes)
{
  # Checks for creating the SingleCellExperiment class
  errors <- character()
  
  CountsBio <- counts(Data)
  GeneNames <- rownames(counts(Data)) 
  BatchInfo <- colData(Data)$BatchInfo

  if (sum(matrixStats::colSums2(CountsBio) == 0) > 0) 
    errors <- c(errors, "Some cells have zero reads mapping back to the 
                intrinsic genes. Please remove them before running the MCMC.\n")
  if (sum(matrixStats::rowSums2(CountsBio) == 0) > 0) 
    warning("Some genes have zero counts across all cells. \n",
            "If comparing 2 groups, use `PriorDelta = 'log-normal' in BASiCS_MCMC.\n",
            "If not, please remove those genes.\n")
  if (!is.null(BatchInfo) & length(BatchInfo) != ncol(Data)) 
    errors <- c(errors, "BatchInfo slot is not compatible with the number of 
                cells contained in the data.\n")
  
  # Checks simplified as this should be already a SingleCellExperiment
  # Also, no longer need to check spike-ins are at the bottom 
  
  if( WithSpikes | (length(altExpNames(Data)) > 0) ) {
    CountsTech <- assay(altExp(Data))
    SpikeInput <- metadata(Data)$SpikeInput
    
    if (!(is.null(SpikeInput)) & !(is.numeric(SpikeInput) & all(SpikeInput > 0) & 
                                   sum(!is.finite(SpikeInput)) == 0)) 
    errors <- c(errors, "Invalid value for 'SpikeInput'.\n")
    
    if (sum(matrixStats::colSums2(CountsTech) == 0) > 0) 
      errors <- c(errors, "Some cells have zero reads mapping back to the 
                  spike-in genes. Please remove these before running the MCMC.\n")
    if(length(SpikeInput) != nrow(assay(altExp(Data)))) 
      errors <- c(errors, "Information provided in 'SpikeInput' is not compatible 
                  with the data provided as 'altExp'.\n")
    if( ncol(Data) != ncol(assay(altExp(Data))) ) 
      errors <- c(errors, "the dimensions of the spike-in data does not match
                  the intrinsic genes dimension")
  } else {
    if (length(unique(BatchInfo)) == 1) 
      errors <- c(errors, "If spike-in genes are not available, BASiCS 
                  requires the data to contain at least 2 batches of cells 
                  (for the same population)\n")
  }

  return(errors)
}