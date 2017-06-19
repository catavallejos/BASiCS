# Combine all Hidden functions in one file

HiddenVarDecomp <- function(Data, object)
{
  
  if(!is(object,"BASiCS_Chain")) stop("'object' is not a BASiCS_Chain class object.")
  
  N = nrow(object@delta); q.bio = ncol(object@delta)
  UniqueBatch = unique(Data@BatchInfo)
  nBatch = length(UniqueBatch)
  
  if(nBatch > 1) {Theta = apply(object@theta, 1, median)}
  else{ Theta = as.vector(object@theta)}
  
  # To store global values (uses median values across all cells)
  PhiS = apply(object@phi *object@s, 1, median)
  Aux = (1/(PhiS * object@mu[,1:q.bio])) + object@delta * (Theta+1)
  TechVarGlobal = Theta / ( Aux + Theta )
  BioVarGlobal = (object@delta * (Theta + 1)) / (Aux + Theta)
  
  # To store batch specific values (in arrays)
  TechVarBatch = array(0, dim = c(N, q.bio, nBatch)) # Technical
  BioVarBatch = array(0, dim = c(N, q.bio, nBatch)) # Biological
  
  if(nBatch > 1)
  {
    for(Batch in 1:nBatch)
    {
      PhiSBatch = apply(object@phi[, Data@BatchInfo == UniqueBatch[Batch]] *
                          object@s[, Data@BatchInfo == UniqueBatch[Batch]], 1, median)
      Aux = (1/(PhiSBatch * object@mu[,1:q.bio])) + object@delta *(object@theta[,Batch]+1)
      TechVarBatch[,,Batch] = object@theta[,Batch] / ( Aux + object@theta[,Batch] )
      BioVarBatch[,,Batch] = (object@delta * (object@theta[,Batch] + 1)) / (Aux + object@theta[,Batch])
    }
  }
  
  if(nBatch > 1)
  {
    list("TechVarGlobal"=TechVarGlobal,
         "BioVarGlobal"=BioVarGlobal,
         "TechVarBatch"=TechVarBatch,
         "BioVarBatch"=BioVarBatch)
  }
  else
  {
    list("TechVarGlobal"=TechVarGlobal,
         "BioVarGlobal"=BioVarGlobal)
  }
}

HiddenTailProbUp<-function(chain,threshold){return(sum(chain>threshold)/length(chain))}
HiddenTailProbLow<-function(chain,threshold){return(sum(chain<threshold)/length(chain))}

HiddenProbHVG<-function(
  VarThreshold, # Variance contribution threshold for HVG/LVG detection. Value must be between 0 and 1.
  VarDecomp) # Output of the variance decomposition obtained using BASiCS_Variance
{
  return(apply(VarDecomp$BioVarGlobal, 2, HiddenTailProbUp, threshold = VarThreshold))
}

HiddenProbLVG<-function(
  VarThreshold, # Variance contribution threshold for HVG/LVG detection. Value must be between 0 and 1.
  VarDecomp) # Output of the variance decomposition obtained using BASiCS_Variance
{
  return(apply(VarDecomp$BioVarGlobal, 2, HiddenTailProbLow, threshold = VarThreshold))
}

HiddenEFDR<-function(
  EviThreshold, # Evidence threshold (it must be contained in (0,1))
  VarThreshold, # Variance contribution threshold choosen by the user (it must be contained in (0,1))
  Prob) # Output of the variance decomposition obtained using BASiCS_Variance
{
  return(sum((1-Prob)*I(Prob>EviThreshold))/sum(I(Prob>EviThreshold)))
}

HiddenEFNR<-function(
  EviThreshold, # Evidence threshold (it must be contained in (0,1))
  VarThreshold, # Variance contribution threshold choosen by the user (it must be contained in (0,1))
  Prob) # Output of the variance decomposition obtained using BASiCS_Variance
{
  return(sum(Prob*I(EviThreshold>=Prob))/sum(I(EviThreshold>=Prob)))
}