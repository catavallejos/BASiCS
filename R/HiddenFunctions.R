#### Helper function for the MCMC simulations

HiddenBASiCS_MCMC_Start<-function(
  Data,
  ...)
{
  if(!is(Data,"SummarizedExperiment")) stop("'Data' is not a SummarizedExperiment class object.")
  
  # Number of instrinsic genes
  q <- length(rowData(Data)$Tech)
  q.bio<-sum(!rowData(Data)$Tech)
  # Number of cells
  n <- dim(assay(Data))[2]
  
  # Initialize normalization as the 'scran' estimates
  sizes.aux = c(20, 40, 60, 80, 100)
  if(n < 200) {sizes.aux = c(20, 40, 60, 80)}
  if(n < 160) {sizes.aux = c(20, 40, 60)}
  if(n < 120) {sizes.aux = c(20, 40)}
  if(n < 80) {sizes.aux = c(20)}
  if(n < 40) {sizes.aux = c(10)}
  size_scran <- scran::computeSumFactors(assay(Data)[!colData(Data)$Tech], sizes = sizes.aux)
  
  if(length(metadata(Data)$SpikeInput) > 1)
  {
    # Initialize s as the empirical capture efficiency rates
    s0 = colSums(assay(Data)[colData(Data)$Tech]) / sum(metadata(Data)$SpikeInput); nu0=s0
    phi0 = size_scran / s0
    phi0 = n * phi0 / sum(phi0)   
    
    # Initialize mu using average 'normalised counts' across cells
    # and true input values for spike-in genes
    nCountsBio <- t( t(assay(Data)[!rowData(Data)$Tech,]) / (phi0*s0) )
    meansBio <- rowMeans( nCountsBio )
    mu0<-c(meansBio + 1,metadata(Data)$SpikeInput) # +1 to avoid zeros as starting values
  }
  else
  {
    phi0 = size_scran
    phi0 = n * phi0 / sum(phi0); 
    for(B in unique(metadata(Data)$BatchInfo))
    {
      aux = metadata(Data)$BatchInfo == B
      phi0[aux] = sum(aux) * phi0[aux] / sum(phi0[aux])
    }
    
    nu0=phi0; s0 = NULL   
    
    # Initialize mu using average 'normalised counts' across cells
    nCountsBio <- t( t(assay(Data)) / phi0 )
    meansBio <- rowMeans( nCountsBio )
    mu0 <- meansBio + 1 # +1 to avoid zeros as starting values    
  }
  
  # Random stating value for delta
  delta0 = rgamma(q.bio,1,1) + 1
  
  # Random stating value for theta (within typically observed range)
  theta0=runif(1, min = 0.2, max = 1)
  
  # If given, load default values for adaptive proposal variances
  args <- list(...)
  ls.mu0 = ifelse("ls.mu0" %in% names(args),args$ls.mu0,-4)
  ls.delta0 = ifelse("ls.delta0" %in% names(args),args$ls.delta0,-2)
  ls.phi0 = ifelse("ls.phi0" %in% names(args),args$ls.phi0,11)
  ls.nu0 = ifelse("ls.nu0" %in% names(args),args$ls.nu0,-10)
  ls.theta0 = ifelse("ls.theta0" %in% names(args),args$ls.theta0,-4)
  
  # Starting values for the proposal variances
  #  ls.mu0 =  pmax(2 * log (0.02 * abs(log(mu0))),ls.mu0)
  #  ls.delta0 =  pmax(2 * log (0.02 * abs(log(delta0))),ls.delta0)
  if(length(metadata(Data)$SpikeInput) > 1) {ls.mu0 = rep(ls.mu0, q)}
  else{ls.mu0 = rep(ls.mu0, q.bio)}
  ls.delta0 = rep(ls.delta0, q.bio)
  ls.phi0 = ifelse(n<200, pmax(2*log(n),ls.phi0), 11) # 6
  ls.nu0 =  pmax(2 * log (0.02 * abs(log(nu0))),ls.nu0)
  ls.theta0 =  pmax(2 * log (0.02 * abs(log(theta0))),ls.theta0)
  
  return(list("mu0"=mu0, "delta0"=delta0, "phi0"=phi0, "s0"=s0, "nu0"=nu0, "theta0"=theta0,
              "ls.mu0"=ls.mu0, "ls.delta0"=ls.delta0, "ls.phi0"=ls.phi0, "ls.nu0"=ls.nu0, "ls.theta0"=ls.theta0))
}


#### Helper functions for differential testing

HiddenTailProbUpDV<-function(chain,threshold){return(sum(chain>threshold)/length(chain))}
HiddenTailProbLowDV<-function(chain,threshold){return(sum(chain<threshold)/length(chain))}

HiddenProbDE<-function(
  chain, # MCMC chain for log-fold change parameter (in absolute value)
  tol) # Minimum tolerance difference for log-fold change
{  
  return(apply(chain, 2, HiddenTailProbUpDV, threshold = tol))
}

HiddenEFDRDV<-function(
  EviThreshold, # Evidence threshold (it must be contained in (0,1))
  Prob) # Probability of changes in expression (for a given tolerance)
{  
  return(sum((1-Prob)*I(Prob > EviThreshold))/sum(I(Prob > EviThreshold)))
}

HiddenEFNRDV<-function(
  EviThreshold, # Evidence threshold (it must be contained in (0,1))
  Prob) # Probability of changes in expression (for a given tolerance) 
{
  return(sum(Prob*I(Prob <= EviThreshold))/sum(I(Prob <= EviThreshold)))
}


#### Helper functions for variance decomposition 

HiddenVarDecomp <- function(Data, object)
{
  
  if(!is(object,"BASiCS_Chain")) stop("'object' is not a BASiCS_Chain class object.")
  
  N = nrow(object@delta); q.bio = ncol(object@delta)
  UniqueBatch = unique(metadata(Data)$BatchInfo)
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
      PhiSBatch = apply(object@phi[, metadata(Data)$BatchInfo == UniqueBatch[Batch]] *
                          object@s[, metadata(Data)$BatchInfo == UniqueBatch[Batch]], 1, median)
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

#### Helper functions for detecting highly and lowly variable genes

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