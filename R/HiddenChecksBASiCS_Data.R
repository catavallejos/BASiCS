HiddenChecksBASiCS_Data <- function(Counts, 
                                    Tech, 
                                    SpikeInput, 
                                    GeneNames,
                                    BatchInfo)
{
  # Checks for creating the SingleCellExperiment class
  errors <- character()
  
  if (!(is.numeric(Counts) & all(Counts >= 0) & sum(!is.finite(Counts)) == 0)) 
    errors <- c(errors, "Invalid value for 'Counts'")
  
  if (sum(Counts%%1) > 0) 
    errors <- c(errors, "Invalid 'Counts' (must be positive integers)")
  
  if (!(is.logical(Tech))) 
    errors <- c(errors, "Invalid value for 'Tech'")
  
  if (!(is.numeric(SpikeInput) & all(SpikeInput > 0) & 
        sum(!is.finite(SpikeInput)) == 0)) 
    errors <- c(errors, "Invalid value for 'SpikeInput'.")
  
  q <- nrow(Counts)
  q.bio <- q - length(SpikeInput)
  n <- ncol(Counts)
  
  # Checks valid for datasets with spikes only
  if (length(SpikeInput) > 1) 
  {
    if (!(length(Tech) == q & sum(!Tech) == q.bio)) 
      errors <- c(errors, "Argument's dimensions are not compatible.")
    
    if (sum(matrixStats::colSums2(Counts[Tech, ]) == 0) > 0) 
      errors <- c(errors, "Some cells have zero reads mapping back to the 
                  spike-in genes. Please remove them before creating the 
                  SingleCellExperiment object.")
    
    if (sum(matrixStats::colSums2(Counts[!Tech, ]) == 0) > 0) 
      errors <- c(errors, "Some cells have zero reads mapping back to the 
                  intrinsic genes. Please remove them before creating the 
                  SingleCellExperiment object.")
    
    if (!(sum(Tech[seq_len(q.bio)]) == 0 & 
          sum(Tech[seq(q.bio + 1, q)]) == q - q.bio)) 
      errors <- c(errors, "Expression counts are not in the right format 
                  (spike-in genes must be at the bottom of the matrix).")
  } 
  else 
  {
    # Checks valid for datasets with no spikes only
    
    if (sum(matrixStats::colSums2(Counts) == 0) > 0) 
      errors <- c(errors, "Some cells have zero reads mapping back to the 
                  intrinsic genes. Please remove them before creating the 
                  SingleCellExperiment object.")
    
    if (length(unique(BatchInfo)) == 1) 
      errors <- c(errors, "If spike-in genes are not available, BASiCS 
                  requires the data to contain at least 2 batches of cells 
                  (for the same population)")
  }
  
  # Checks valid for any data
  if (length(Tech) != q) 
    errors <- c(errors, "Argument's dimensions are not compatible.")
  
  if (length(GeneNames) != q) 
    errors <- c(errors, "Incorrect length of GeneNames.")
  
  if (sum(matrixStats::rowSums2(Counts) == 0) > 0) 
    warning("Some genes have zero counts across all cells. \n",
            "If comparing 2 groups, use `PriorDelta = 'log-normal' in BASiCS_MCMC.\n",
            "If not, please remove those genes.")
  
  if (length(BatchInfo) != n) 
    errors <- c(errors, "BatchInfo slot is not compatible with the number of 
                cells contained in Counts slot.")
  
  return(errors)
}