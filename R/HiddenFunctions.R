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
  size_scran <- scran::computeSumFactors(as.matrix(assay(Data)[!rowData(Data)$Tech,,drop=FALSE]), sizes = sizes.aux)
  
  if(length(metadata(Data)$SpikeInput) > 1)
  {
    # Initialize s as the empirical capture efficiency rates
    s0 = colSums(assay(Data)[rowData(Data)$Tech,]) / sum(metadata(Data)$SpikeInput); nu0=s0
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
    # +1 to avoid zeros as starting values 
    meansBio <- ifelse(meansBio == 0, meansBio + 1, meansBio)
    mu0 <- meansBio 
  }
  
  # Random stating value for delta
  # Defined by the CV for high- and mid-expressed genes
  # This comes from equation (2) in Vallejos et al (2016)
  varsBio <- apply( nCountsBio, 1, var )
  cv2Bio <- varsBio / (meansBio)^2
  delta0 = rgamma(q.bio,1,1) + 1
  Aux = which(meansBio > quantile(meansBio, 0.10))
  delta0[Aux] <- cv2Bio[Aux]
  # 1e-3 added to be coherent with tolerance used within MCMC sampler
  delta0 = delta0 + 1e-3
  
  # Random stating value for theta (within typically observed range)
  theta0 = runif(1, min = 0.2, max = 1)
  
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

HiddenHeaderDetectHVG_LVG <- function(Data,
                                      object,
                                      VarThreshold,
                                      EviThreshold = NULL,
                                      EFDR = 0.05, 
                                      OrderVariable = "Prob",
                                      Plot = FALSE)
{
  if(!is(Data,"SummarizedExperiment")) stop("'Data' is not a SummarizedExperiment class object. Please use the 'newBASiCS_Data' function to create a SummarizedExperiment object.")
  if(!is(object,"BASiCS_Chain")) stop("'object' is not a BASiCS_Chain class object.")
  if(VarThreshold<0 | VarThreshold>1 | !is.finite(VarThreshold)) stop("Variance contribution thresholds for HVG/LVG detection must be contained in (0,1)")
  if(!is.logical(Plot) | length(Plot)!=1) stop("Please insert TRUE or FALSE for Plot parameter")
  if(!is.null(EviThreshold))
  {
    if(EviThreshold<0 | EviThreshold>1 | !is.finite(EviThreshold))
      stop("Evidence thresholds for HVG and LVG detection must be contained in (0,1) \n For automatic threshold search use EviThreshold = NULL.")
  }
  if(!(OrderVariable %in% c("GeneNames", "Mu", "Delta", "Sigma", "Prob"))) stop("Invalid 'OrderVariable' value")
  if(is.null(EviThreshold)) {message(paste("Posterior probability threshold to be defined by EFDR = ", 100*EFDR, "% (+-2.5% tolerance) ..."))}
}


HiddenThresholdSearchDetectHVG_LVG <- function(EviThreshold,
                                               VarThreshold,
                                               Prob, 
                                               EFDR)
{
  # If EviThreshold is not set a priori (search)
  if(length(EviThreshold) == 0)
  {
    EviThresholds <- seq(0.5,0.9995,by=0.0005)
    
    EFDRgrid <- sapply(EviThresholds, HiddenEFDR, VarThreshold = VarThreshold, Prob = Prob)
    EFNRgrid <- sapply(EviThresholds, HiddenEFNR, VarThreshold = VarThreshold, Prob = Prob)
    
    above <- abs(EFDRgrid - EFDR)
    
    if(sum(!is.na(above)) > 0)
    {
      # Search EFDR closest to the desired value
      EFDRopt <- EFDRgrid[above == min(above, na.rm = TRUE) & !is.na(above)] 
      # If multiple threholds lead to same EFDR, choose the one with the lowest EFDR
      EFNRopt <- EFNRgrid[EFDRgrid == mean(EFDRopt) & !is.na(EFDRgrid)]
      if(sum(!is.na(EFNRopt)) > 0)
      {
        optimal <- which(EFDRgrid == mean(EFDRopt) & EFNRgrid == mean(EFNRopt))
      } 
      else
      {
        optimal <- which(EFDRgrid == mean(EFDRopt))
      }
      # Quick fix for EFDR/EFNR ties; possibly not an issue in real datasets
      optimal <- median(round(median(optimal)))
      OptThreshold <- c(EviThresholds[optimal], EFDRgrid[optimal], EFNRgrid[optimal])
      
      if(OptThreshold[1] < 0.5) 
      {
        message("For the given variance contribution threshold, the evidence threshold 
                that achieves the desired EFDR is below 0.5. By default, the evidence
                threshold will be set at 0.5 (corresponding EFDR/EFNR reported)")
        OptThreshold <- c(EviThresholds[1], EFDRgrid[1], EFNRgrid[1])
      }
      
      # Message when different to desired EFDR is large
      if( abs(OptThreshold[2] - EFDR) > 0.025 )
      {
        message("For the given variance contribution threshold, it is not possible 
                to find an evidence threshold (>0.5) that achieves the desired EFDR level 
                (tolerance +- 0.025). The output of this function reflects the 
                closest possible value. \n")         
      }  
      }
    else
    {
      message("For the given variance contribution threshold, it is not possible 
              to estimate EFDR. By default, the evidence
              threshold will be set at 0.5. \n")    
      OptThreshold <- c(0.5, NA, NA)  
    }
  }
  # If EviThreshold is set a priori
  else
  {
    EFDR = HiddenEFDR(EviThreshold, VarThreshold, Prob)
    EFNR = HiddenEFNR(EviThreshold, VarThreshold, Prob)
    OptThreshold <- c(EviThreshold, EFDR, EFNR)
  }  
  list("EviThresholds" = EviThresholds,
       "EFDRgrid" = EFDRgrid,
       "EFNRgrid" = EFNRgrid,
       "OptThreshold" = OptThreshold)
}

HiddenPlot1DetectHVG_LVG <- function(EviThresholds, 
                                     EFDRgrid,
                                     EFNRgrid,
                                     OptThreshold,
                                     EviThreshold,
                                     EFDR)
{
  plot(EviThresholds, EFDRgrid, type = "l", lty = 1, bty = "n", lwd = 2,  
       ylab = "Error rate", xlab = "Evidence threshold", ylim = c(0,1))
  lines(EviThresholds, EFNRgrid, lty = 2, lwd = 2)
  if(length(EviThreshold) == 0) {abline(h = EFDR, col = "blue", lwd = 2, lty = 1)}
  abline(v = OptThreshold[1], col = "red", lwd = 2, lty = 1)
  if(length(EviThreshold) == 0)
  {
    legend('topleft', c("EFDR", "EFNR", "Target EFDR"), lty = c(1:2,1), 
           col = c("black", "black", "blue"), bty = "n", lwd = 2)  
  }
}

HiddenPlot2DetectHVG_LVG <- function(args, 
                                     Task, 
                                     Mu, 
                                     Prob,
                                     OptThreshold,
                                     Hits)
{
  if("ylim" %in% names(args)) {ylim = args$ylim} else{ylim = c(0, 1)}
  if("xlim" %in% names(args)) {xlim = args$xlim} else{xlim = c(min(Mu),max(Mu))}
  cex = ifelse("cex" %in% names(args),args$cex, 1.5)
  pch = ifelse("pch" %in% names(args),args$pch, 16)
  col = ifelse("col" %in% names(args),args$col, 8)
  bty = ifelse("bty" %in% names(args),args$bty, "n")
  cex.lab = ifelse("cex.lab" %in% names(args),args$cex.lab, 1)
  cex.axis = ifelse("cex.axis" %in% names(args),args$cex.axis, 1)
  cex.main = ifelse("cex.main" %in% names(args),args$cex.main, 1)
  xlab = ifelse("xlab" %in% names(args),args$xlab, "Mean expression")
  if(Task == "HVG") {ylab = ifelse("ylab" %in% names(args),args$ylab, "HVG probability")}
  else {ylab = ifelse("ylab" %in% names(args),args$ylab, "LVG probability")}
  main = ifelse("main" %in% names(args),args$main, "")
  
  plot(Mu, Prob, log="x", pch = pch, ylim = ylim, xlim = xlim, col = col, cex = cex,
       bty = bty, cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main,
       xlab = xlab, ylab = ylab, main = main)
  abline(h = OptThreshold[1], lty = 2, col = "black")
  points(Mu[Hits], Prob[Hits], pch = pch, col = "red", cex = cex)
}

# Used within BASiCS_TestDE function
HiddenThresholdSearchTestDE <- function(ChainLFC, 
                                        Epsilon,
                                        EviThreshold,
                                        GenesSelect,
                                        EFDR,
                                        Task)
{
  # Calculating posterior probabilities
  if(Epsilon > 0) {Prob = HiddenProbDE(chain = abs(ChainLFC), tol = Epsilon)}
  else 
  {
    Prob_aux = HiddenProbDE(chain = ChainLFC, tol = 0)
    Prob = 2*pmax(Prob_aux, 1-Prob_aux) - 1
  }
  
  # Posterior probability threshold search
  if(is.null(EviThreshold))
  {
    EviThresholds <- seq(0.5,0.9995,by=0.00025)
    
    if(is.null(GenesSelect))
    {
      EFDRgrid <- sapply(EviThresholds, HiddenEFDRDV, Prob = Prob)
      EFNRgrid <- sapply(EviThresholds, HiddenEFNRDV, Prob = Prob)      
    }
    else
    {
      EFDRgrid <- sapply(EviThresholds, HiddenEFDRDV, Prob = Prob[GenesSelect])
      EFNRgrid <- sapply(EviThresholds, HiddenEFNRDV, Prob = Prob[GenesSelect])        
    }
    
    optimal = round(median(which(abs(EFDRgrid - EFDR) == min(abs(EFDRgrid - EFDR)))))
    OptThreshold <- c(EviThresholds[optimal], EFDRgrid[optimal], EFNRgrid[optimal])
    
    EviThreshold = OptThreshold[1]
    if(is.na(EviThreshold)) 
    { 
      message(paste(Task, ": EFDR calibration failed. Probability threshold automatically set equal to 0.90 \n"))
      EviThreshold = 0.90
    }
  }
  else
  {
    if(is.null(GenesSelect))
    {
      EFDRgrid = HiddenEFDRDV(EviThreshold, Prob)
      EFNRgrid = HiddenEFNRDV(EviThreshold, Prob)
    }
    else
    {
      EFDRgrid = HiddenEFDRDV(EviThreshold, Prob[GenesSelect])
      EFNRgrid = HiddenEFNRDV(EviThreshold, Prob[GenesSelect])      
    }
    OptThreshold <- c(EviThreshold, EFDRgrid, EFNRgrid)
  }  
  list("Prob" = Prob,
       "OptThreshold" = OptThreshold,
       "EFDRgrid" = EFDRgrid,
       "EFNRgrid" = EFNRgrid)
}
