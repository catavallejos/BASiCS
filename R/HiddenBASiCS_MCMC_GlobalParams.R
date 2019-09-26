HiddenBASiCS_MCMC_GlobalParams <- function(Data, BatchInfoColumn) {
  
  BioCounts <- counts(Data)
  n <- ncol(Data)
  q.bio <- nrow(Data)
  
  # If data contains spike-ins
  if (length(altExpNames(Data)) > 1){
    q <- q.bio + nrow(altExp(Data))
  } else { q <- q.bio }
  
  # If Data contains batches
  if(!is.null(colData(Data)[[BatchInfoColumn]])){
    # Store the correct number of levels in batch vector
    BatchInfo <- as.factor(colData(Data)[[BatchInfoColumn]])
    BatchInfo <- factor(BatchInfo, levels = unique(BatchInfo))
    nBatch <- length(unique(BatchInfo))
  } else {
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
  
  list(
    BioCounts = as.matrix(BioCounts),
    q = q,
    q.bio = q.bio,
    n = n,
    nBatch = nBatch,
    BatchInfo = BatchInfo,
    BatchDesign = BatchDesign
  )
}
