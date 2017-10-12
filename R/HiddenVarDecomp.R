# Used in BASiCS_VarianceDecomp
HiddenVarDecomp <- function(Chain) 
{
  if (!is(Chain, "BASiCS_Chain")) 
    stop("'Chain' is not a BASiCS_Chain class object.")
    
  N <- nrow(Chain@parameters$delta)
  q.bio <- ncol(Chain@parameters$delta)
  UniqueBatch <- colnames(Chain@parameters$theta)
  nBatch <- length(UniqueBatch)
  CellName <- colnames(Chain@parameters$phi)
    
  if (nBatch > 1) { Theta <- matrixStats::rowMedians(Chain@parameters$theta) } 
  else { Theta <- as.vector(Chain@parameters$theta) }
    
  # To store global values (uses median values across all cells)
  PhiS <- matrixStats::rowMedians(Chain@parameters$phi * Chain@parameters$s)
  Aux <- (1/(PhiS * Chain@parameters$mu)) + Chain@parameters$delta * (Theta + 1)
  TechVarGlobal <- Theta/(Aux + Theta)
  BioVarGlobal <- (Chain@parameters$delta * (Theta + 1))/(Aux + Theta)
    
  # To store batch specific values (in arrays)
  TechVarBatch <- array(0, dim = c(N, q.bio, nBatch))  # Technical
  BioVarBatch <- array(0, dim = c(N, q.bio, nBatch))  # Biological
    
  if (nBatch > 1) 
  {
    for (Batch in seq_len(nBatch)) 
    {
      PhiBatch <- Chain@parameters$phi[, grep(UniqueBatch[Batch], CellName)]
      SBatch <- Chain@parameters$s[, grep(UniqueBatch[Batch], CellName)]
      PhiSBatch <- matrixStats::rowMedians(PhiBatch * SBatch)
      
      Aux <- (1/(PhiSBatch * Chain@parameters$mu)) + 
                Chain@parameters$delta * (Chain@parameters$theta[,Batch] + 1)
      TechVarBatch[,,Batch] <- Chain@parameters$theta[,Batch] / (Aux + Chain@parameters$theta[,Batch])
      BioVarBatch[,,Batch] <- (Chain@parameters$delta * (Chain@parameters$theta[,Batch] + 1)) / 
                                    (Aux + Chain@parameters$theta[,Batch])
    }
  }
    
  if (nBatch > 1) 
  {
    list(TechVarGlobal = TechVarGlobal, 
         BioVarGlobal = BioVarGlobal, 
         TechVarBatch = TechVarBatch, 
         BioVarBatch = BioVarBatch)
  } 
  else { list(TechVarGlobal = TechVarGlobal, BioVarGlobal = BioVarGlobal) }
}
