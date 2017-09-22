# Used in BASiCS_VarianceDecomp
HiddenVarDecomp <- function(Chain) 
{
  if (!is(Chain, "BASiCS_Chain")) 
    stop("'Chain' is not a BASiCS_Chain class object.")
    
  N <- nrow(Chain@delta)
  q.bio <- ncol(Chain@delta)
  UniqueBatch <- colnames(Chain@theta)
  nBatch <- length(UniqueBatch)
  CellName <- colnames(Chain@phi)
    
  if (nBatch > 1) { Theta <- matrixStats::rowMedians(Chain@theta) } 
  else { Theta <- as.vector(Chain@theta) }
    
  # To store global values (uses median values across all cells)
  PhiS <- matrixStats::rowMedians(Chain@phi * Chain@s)
  Aux <- (1/(PhiS * Chain@mu)) + Chain@delta * (Theta + 1)
  TechVarGlobal <- Theta/(Aux + Theta)
  BioVarGlobal <- (Chain@delta * (Theta + 1))/(Aux + Theta)
    
  # To store batch specific values (in arrays)
  TechVarBatch <- array(0, dim = c(N, q.bio, nBatch))  # Technical
  BioVarBatch <- array(0, dim = c(N, q.bio, nBatch))  # Biological
    
  if (nBatch > 1) 
  {
    for (Batch in seq_len(nBatch)) 
    {
      PhiBatch <- Chain@phi[, grep(UniqueBatch[Batch], CellName)]
      SBatch <- Chain@s[, grep(UniqueBatch[Batch], CellName)]
      PhiSBatch <- matrixStats::rowMedians(PhiBatch * SBatch)
      
      Aux <- (1/(PhiSBatch * Chain@mu)) + 
                Chain@delta * (Chain@theta[,Batch] + 1)
      TechVarBatch[,,Batch] <- Chain@theta[,Batch] / (Aux + Chain@theta[,Batch])
      BioVarBatch[,,Batch] <- (Chain@delta * (Chain@theta[,Batch] + 1)) / 
                                    (Aux + Chain@theta[,Batch])
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
