HiddenBASiCS_MCMC_GlobalParams <- function(Data) {
  
  # If data contains spike-ins
  if (!is.null(SingleCellExperiment::isSpike(Data))) {
    q <- nrow(Data)
    q.bio <- sum(!SingleCellExperiment::isSpike(Data))
    spikes <- SingleCellExperiment::isSpike(Data)
    BioCounts <- counts(Data)[!spikes, ]
  }
  else {
    BioCounts <- counts(Data)
    q.bio <- q <- nrow(Data)
    spikes <- rep(FALSE, q)
  }
  
  n <- dim(counts(Data))[2]
  
  # If Data contains batches
  if(!is.null(colData(Data)$BatchInfo)){
    # Store the correct number of levels in batch vector
    BatchInfo <- as.factor(colData(Data)$BatchInfo)
    BatchInfo <- factor(BatchInfo, levels = unique(BatchInfo))
    nBatch <- length(unique(BatchInfo))
  }
  else{
    BatchInfo <- rep(1, times = n)
    nBatch <- 1
  }
  
  # Parameters associated to the presence of batches
  if(nBatch > 1) {
    BatchDesign <- model.matrix(~BatchInfo - 1)  
  } else { 
    # If there are no batches or the user asked to ignore them
    BatchDesign <- matrix(1, nrow = n, ncol = 1) 
  }
  
  list(BioCounts = as.matrix(BioCounts), q = q, q.bio = q.bio, n = n,
       nBatch = nBatch, BatchInfo = BatchInfo, BatchDesign = BatchDesign)
}