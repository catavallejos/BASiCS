.ChecksBASiCS_Data <- function(Data, WithSpikes) {

  # Checks for creating the SingleCellExperiment class
  errors <- character()
  
  CountsBio <- counts(Data)
  GeneNames <- rownames(counts(Data)) 
  BatchInfo <- colData(Data)$BatchInfo

  if (sum(Matrix::colSums(CountsBio) == 0)) {
    errors <- c(errors,
      "Some cells have zero reads mapping back to the 
      intrinsic genes. Please remove them before running the MCMC.\n"
    )
  }
  if (sum(Matrix::rowSums(CountsBio) == 0)) {
    warning(
      "Some genes have zero counts across all cells. \n",
      "If comparing 2 groups, use `PriorDelta = 'log-normal' in BASiCS_MCMC.\n",
      "If not, please remove those genes.\n"
    )
  }
  if (!is.null(BatchInfo) & length(BatchInfo) != ncol(Data)) {
    errors <- c(errors,
      "BatchInfo slot is not compatible with the number of 
      cells contained in the data.\n"
    )
  }
  
  # Checks simplified as this should be already a SingleCellExperiment
  # Also, no longer need to check spike-ins are at the bottom 
  
  if (!WithSpikes) {
    if (length(unique(BatchInfo)) == 1) {
      errors <- c(errors,
        "If spike-in genes are not available, BASiCS 
        requires the data to contain at least 2 batches of cells 
        (for the same population)\n"
      )
    }
  } else {
    if (length(altExpNames(Data)) != 1) {
      if (length(altExpNames(Data)) > 1) {
        errors <- c(errors,
          "More than one 'altExp' provided; only one allowed. \n
          'altExp' must contain observations for spike-in genes. \n"
        )
      } else {
        errors <- c(errors, "spike-in counts must be provided as 'altExp'. \n")
      }
    } else {
      if (!("SpikeInput" %in% names(metadata(Data)))) {
        errors <- c(errors, "'SpikeInput' was not provided as metadata.\n")
      }
      # Extract spike-ins
      CountsTech <- assay(altExp(Data))
      SpikeInput <- metadata(Data)$SpikeInput
      # Validity checks for SpikeInput
      if (!is.data.frame(SpikeInput)) {
        errors <- c(errors, "'SpikeInput' must be a 'data.frame'.\n")
      }
      if (ncol(SpikeInput) != 2) {
        errors <- c(errors, "'SpikeInput' must have two columns only.\n")
      }
      if (nrow(SpikeInput) != nrow(CountsTech)) {
        errors <- c(errors,
          "'SpikeInput' dimensions not compatible with 'altExp'. \n"
        )
      }
      # Validity checks for input concentrations
      if (!(is.numeric(SpikeInput[, 2]) & all(SpikeInput[, 2] > 0) & 
            sum(!is.finite(SpikeInput[, 2])) == 0)) {
        errors <- c(errors,
          "Invalid value in the 2nd column of 'SpikeInput'.\n"
        )
      }
      # Check order in SpikeInput
      if (any(SpikeInput[, 1] != rownames(altExp(Data)))) {
        errors <- c(errors, "'SpikeInput' row order does not match 'altExp'.\n")
      }
      # Check all cells have non-zero total count
      if (sum(Matrix::colSums(CountsTech) == 0) > 0) {
        errors <- c(errors,
          "Some cells have zero reads mapping back to the 
          spike-in genes. Please remove these before running the MCMC.\n"
        )
      }
      # Check all genes have non-zero total count
      if (sum(Matrix::rowSums(CountsTech) == 0) > 0) {
        errors <- c(errors,
          "Some spike-in genes have zero total reads
          across all cells. Please remove these before running the MCMC.\n"
        )
      }
    }
  }

  return(errors)
}
