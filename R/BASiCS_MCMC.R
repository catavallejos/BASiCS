#' @title BASiCS MCMC sampler
#'
#' @description MCMC sampler to perform Bayesian inference for single-cell mRNA sequencing datasets using the model described in Vallejos et al (2015).
#'
#' @param Data A \code{\link[SummarizedExperiment]{SummarizedExperiment}} object.
#' @param N Total number of iterations for the MCMC sampler. Use \code{N>=max(4,Thin)}, \code{N} being a multiple of \code{Thin}.
#' @param Thin Thining period for the MCMC sampler. Use \code{Thin>=2}.
#' @param Burn Burn-in period for the MCMC sampler. Use \code{Burn>=1}, \code{Burn<N}, \code{Burn} being a multiple of \code{Thin}.
#' @param ... Optional parameters.
#' \describe{
#' \item{\code{PriorParam}}{List of 7 elements, containing the hyper-parameter values required for the adopted prior (see Vallejos et al, 2015). All elements must be positive real numbers.
#' \describe{
#'   \item{\code{s2.mu}}{Scale hyper-parameter for the log-Normal(\code{0},\code{s2.mu}) prior that is shared by all gene-specific expression rate parameters \eqn{\mu[i]}.
#'   Default: \code{s2.mu = 0.5}.}
#'   \item{\code{a.delta}}{Only used when `PriorDelta == 'gamma'`. Shape hyper-parameter for the Gamma(\code{a.delta},\code{b.delta}) prior that is shared by all gene-specific biological cell-to-cell heterogeneity hyper-parameters \eqn{\delta[i]}.
#'   Default: \code{a.delta = 1}.}
#'   \item{\code{b.delta}}{Only used when `PriorDelta == 'gamma'`. Rate hyper-parameter for the Gamma(\code{a.delta},\code{b.delta}) prior that is shared by all gene-specific biological cell-to-cell heterogeneity hyper-parameters \eqn{\delta[i]}.
#'   Default: \code{b.delta = 1}.}
#'   \item{\code{p.phi}}{Dirichlet hyper-parameter for the joint of all (scaled by \code{n}) cell-specific mRNA content normalising constants \eqn{\phi[j] / n}.
#'   Default: \code{p.phi = rep(1, n)}.}
#'   \item{\code{a.s}}{Shape hyper-parameter for the Gamma(\code{a.s},\code{b.s}) prior that is shared by all cell-specific capture efficiency normalising constants \eqn{s[j]}.
#'   Default: \code{a.s = 1}.}
#'   \item{\code{b.s}}{Rate hyper-parameter for the Gamma(\code{a.s},\code{b.s}) prior that is shared by all cell-specific capture efficiency normalising constants \eqn{s[j]}.
#'   Default: \code{b.s = 1}.}
#'   \item{\code{a.theta}}{Shape hyper-parameter for the Gamma(\code{a.theta},\code{b.theta}) prior for technical noise hyper-parameter \eqn{\theta}.
#'   Default: \code{a.theta = 1}.}
#'   \item{\code{b.theta}}{Rate hyper-parameter for the Gamma(\code{a.theta},\code{b.theta}) prior for technical noise hyper-parameter \eqn{\theta}.
#'   Default: \code{b.theta = 1}.}
#'   \item{\code{s2.delta}}{Only used when `PriorDelta == 'log-normal'`. Scale hyper-parameter for the log-Normal(\code{0},\code{s2.delta}) prior that is shared by all gene-specific expression rate parameters \eqn{\delta[i]}.
#'   Default: \code{s2.delta = 0.5}. }
#' }}
#' \item{\code{AR}}{Optimal acceptance rate for adaptive Metropolis Hastings updates. It must be a positive number between 0 and 1. Default (and recommended): \code{ar = 0.44}}.
#'
#' \item{\code{StopAdapt}}{Iteration at which adaptive proposals are not longer adapted. Use \code{StopAdapt>=1}. Default: \code{StopAdapt = Burn}.}
#'
#' \item{\code{StoreChains}}{If \code{StoreChains = T}, the slots of the generated \code{BASiCS_Chain} object are stored in separate .txt files. Each row of the output file containing an interation (\code{RunName} argument used to index file names). Default: \code{StoreChains = F}.}
#' \item{\code{StoreAdapt}}{If \code{StoreAdapt = T}, trajectory of adaptive proposal variances (in log-scale) for each parameter are stored in separate .txt files. Each row of the output file containing an interation (\code{RunName} argument used to index file names). Default: \code{StoreAdapt = F}.}
#' \item{\code{StoreDir}}{Directory where output files are stored. Only required if \code{StoreChains = TRUE} and/or \code{StoreAdapt = TRUE}). Default: \code{StoreDir = getwd()}.}
#' \item{\code{RunName}}{String used to index `.txt` files storing chains and/or adaptive proposal variances.}
#' \item{\code{PrintProgress}}{If \code{PrintProgress = FALSE}, console-based progress report is suppressed.}
#' \item{\code{ls.phi0}}{Starting value for the adaptive concentration parameter of the Metropolis proposals for \code{phi}.}
#' \item{\code{PriorDelta}}{Specifies the prior used for \code{delta}. Possible values are 'gamma' (Gamma(\code{a.theta},\code{b.theta}) prior) and 'log-normal' (log-Normal(\code{0},\code{s2.delta}) prior) .}
#' \item{\code{Start}}{In general, we do not advise to specify this argument. Default options have been tuned to facilitate convergence. It can be used to set user defined starting points for the MCMC algorithm. If used, it must be a list containing the following elements: \code{mu0},
#' \code{delta0}, \code{phi0}, \code{s0}, \code{nu0}, \code{theta0}, \code{ls.mu0}, \code{ls.delta0}, \code{ls.phi0}, \code{ls.nu0}, \code{ls.theta0}}
#' }
#'
#' @return An object of class \code{\link[BASiCS]{BASiCS_Chain-class}}.
#'
#' @examples
#'
#' # Built-in simulated dataset
#' Data = makeExampleBASiCS_Data()
#' # For real data use the newBASiCS_Data(Counts, Tech, SpikeInfo, BatchInfo) function
#'
#' # Only a short run of the MCMC algorithm for illustration purposes
#' # Longer runs migth be required to reach convergence
#' MCMC_Output <- BASiCS_MCMC(Data, N = 10000, Thin = 10, Burn = 5000, PrintProgress = FALSE)
#' head(displayChainBASiCS(MCMC_Output, Param = "mu"))
#' head(displayChainBASiCS(MCMC_Output, Param = "delta"))
#' head(displayChainBASiCS(MCMC_Output, Param = "phi"))
#' head(displayChainBASiCS(MCMC_Output, Param = "s"))
#' head(displayChainBASiCS(MCMC_Output, Param = "nu"))
#' head(displayChainBASiCS(MCMC_Output, Param = "theta"))
#'
#' # Traceplots
#' plot(MCMC_Output, Param = "mu", Gene = 1)
#' plot(MCMC_Output, Param = "delta", Gene = 1)
#' plot(MCMC_Output, Param = "phi", Cell = 1)
#' plot(MCMC_Output, Param = "s", Cell = 1)
#' plot(MCMC_Output, Param = "nu", Cell = 1)
#' plot(MCMC_Output, Param = "theta", Batch = 1)
#'
#' # Calculating posterior medians and 95% HPD intervals
#' MCMC_Summary <- Summary(MCMC_Output)
#' head(displaySummaryBASiCS(MCMC_Summary, Param = "mu"))
#' head(displaySummaryBASiCS(MCMC_Summary, Param = "delta"))
#' head(displaySummaryBASiCS(MCMC_Summary, Param = "phi"))
#' head(displaySummaryBASiCS(MCMC_Summary, Param = "s"))
#' head(displaySummaryBASiCS(MCMC_Summary, Param = "nu"))
#' head(displaySummaryBASiCS(MCMC_Summary, Param = "theta"))
#'
#' # Graphical display of posterior medians and 95% HPD intervals
#' plot(MCMC_Summary, Param = "mu", main = "All genes")
#' plot(MCMC_Summary, Param = "mu", Genes = 1:10, main = "First 10 genes")
#' plot(MCMC_Summary, Param = "delta", main = "All genes")
#' plot(MCMC_Summary, Param = "delta", Genes = c(2,5,10,50), main = "5 customized genes")
#' plot(MCMC_Summary, Param = "phi", main = "All cells")
#' plot(MCMC_Summary, Param = "phi", Cells = 1:5, main = "First 5 cells")
#' plot(MCMC_Summary, Param = "s", main = "All cells")
#' plot(MCMC_Summary, Param = "s", Cells = 1:5, main = "First 5 cells")
#' plot(MCMC_Summary, Param = "nu", main = "All cells")
#' plot(MCMC_Summary, Param = "nu", Cells = 1:5, main = "First 5 cells")
#' plot(MCMC_Summary, Param = "theta")
#'
#' # To constrast posterior medians of cell-specific parameters
#' plot(MCMC_Summary, Param = "phi", Param2 = "s")
#' plot(MCMC_Summary, Param = "phi", Param2 = "nu")
#' plot(MCMC_Summary, Param = "s", Param2 = "nu")
#'
#' # To constrast posterior medians of gene-specific parameters
#' plot(MCMC_Summary, Param = "mu", Param2 = "delta", log = "x")
#'
#' # Highly and lowly variable genes detection
#' #DetectHVG <- BASiCS_DetectHVG(Data, MCMC_Output, VarThreshold = 0.70, Plot = TRUE)
#' #DetectLVG <- BASiCS_DetectLVG(Data, MCMC_Output, VarThreshold = 0.40, Plot = TRUE)
#'
#' #plot(MCMC_Summary, Param = "mu", Param2 = "delta", log = "x", col = 8)
#' #points(DetectHVG$Table[DetectHVG$Table$HVG==1,2], DetectHVG$Table[DetectHVG$Table$HVG==1,3],
#' #       pch = 16, col = "red", cex = 1)
#' #points(DetectLVG$Table[DetectLVG$Table$LVG==1,2], DetectLVG$Table[DetectLVG$Table$LVG==1,3],
#' #       pch = 16, col = "blue", cex = 1)
#'
#' # If variance thresholds are not fixed
#' #BASiCS_VarThresholdSearchHVG(Data, MCMC_Output, VarThresholdsGrid = seq(0.70,0.75,by=0.01))
#' #BASiCS_VarThresholdSearchLVG(Data, MCMC_Output, VarThresholdsGrid = seq(0.40,0.45,by=0.01))
#'
#' @author Catalina A. Vallejos \email{cnvallej@@uc.cl}
#'
#' @references Vallejos, Marioni and Richardson (2015). Bayesian Analysis of Single-Cell Sequencing data. PloS Computational Biology.
#' @references Vallejos, Richardson and Marioni (2016). Beyond comparisons of means: understanding changes in gene expression at the single cell level. Genome Biology.

# Change Data class to SE class

BASiCS_MCMC <- function(
  Data,
  N,
  Thin,
  Burn,
  ...)
{
  
  
  if(!is(Data,"SummarizedExperiment")) stop("'Data' is not a SummarizedExperiment class object.")
  
  # SOME QUANTITIES USED THROUGHOUT THE MCMC ALGORITHM
  q=length(rowData(Data)$Tech); q.bio=sum(!rowData(Data)$Tech); n=dim(assay(Data))[2]
  
  args <- list(...)
  
  if(!("PriorDelta" %in% names(args))) {message(" --------------------------------------------------------------------------- \n IMPORTANT: by default, the argument PriorDelta was set equal to 'gamma'. \n When performing a differential over-dispersion test between two populations, \n we recommend to change this value to PriorDelta = 'log-normal' \n --------------------------------------------------------------------------- \n" )}
  
  
  if("PriorParam" %in% names(args)) {PriorParam = args$PriorParam}
  else { PriorParam = list(s2.mu = 0.5, s2.delta = 0.5, a.delta = 1, b.delta = 1, p.phi = rep(1, times = n), a.phi = 1, b.phi = 1, a.s = 1, b.s = 1, a.theta = 1, b.theta = 1)}
  AR = ifelse("AR" %in% names(args),args$AR, 0.44)
  StopAdapt = ifelse("StopAdapt" %in% names(args),args$StopAdapt, Burn)
  StoreChains = ifelse("StoreChains" %in% names(args),args$StoreChains, F)
  StoreAdapt = ifelse("StoreAdapt" %in% names(args),args$StoreAdapt, F)
  StoreDir = ifelse("StoreDir" %in% names(args),args$StoreDir, getwd())
  RunName = ifelse("RunName" %in% names(args),args$RunName, "")
  PrintProgress = ifelse("PrintProgress" %in% names(args),args$PrintProgress, TRUE)
  PriorDelta = ifelse("PriorDelta" %in% names(args), args$PriorDelta, "gamma")
  
  if(!(length(N) == 1 | length(Thin) == 1 | length(Burn) == 1)) stop("Invalid parameter values.")
  if(!(N%%Thin==0 & N>=max(4,Thin))) stop("Please use an integer value for N. It must also be a multiple of thin (N>=4)).")
  if(!(Thin%%1==0 & Thin>=2)) stop("Please use an integer value for Thin (Thin>=2).")
  if(!(Burn%%Thin==0 & Burn<N & Burn>=1)) stop("Please use an integer value for Burn. It must also be lower than N and a multiple of thin (Burn>=1).")
  
  if(!(PriorParam$s2.mu>0  & length(PriorParam$s2.mu) == 1 &
       PriorParam$s2.delta>0  & length(PriorParam$s2.delta) == 1 &
       PriorParam$a.delta>0  & length(PriorParam$a.delta) == 1 &
       PriorParam$b.delta>0  & length(PriorParam$b.delta) == 1 &
       all(PriorParam$p.phi>0) & length(PriorParam$p.phi) == n &
       PriorParam$a.s>0      & length(PriorParam$a.s) == 1 &
       PriorParam$b.s>0      & length(PriorParam$b.s) == 1 &
       PriorParam$a.theta>0  & length(PriorParam$a.theta) == 1 &
       PriorParam$b.theta>0) & length(PriorParam$b.theta) == 1) stop("Invalid prior hyper-parameter values.")
  
  if(!(AR>0 & AR<1 & length(AR) == 1)) stop("Invalid AR value. Recommended value: AR = 0.44.")
  if(!(StopAdapt>0)) stop("Invalid StopAdapt value.")
  if(!(is.logical(StoreChains) & length(StoreChains) == 1)) stop("Invalid StoreChains value.")
  if(!(is.logical(StoreAdapt) & length(StoreAdapt) == 1)) stop("Invalid StoreAdapt value.")
  if(!(file.info(StoreDir)["isdir"])) stop("Invalid StoreDir value.")
  if(!(PriorDelta %in% c("gamma","log-normal"))) stop("Invalid PriorDelta value.")
  
  PriorDeltaNum = ifelse(PriorDelta == "gamma", 1, 2)
  
  # SOME SUMS USED THROUGHOUT THE MCMC ALGORITHM
  sum.bycell.all<-apply(assay(Data),1,sum)
  sum.bycell.bio<-apply(assay(Data)[1:q.bio,],1,sum)
  sum.bygene.all<-apply(assay(Data),2,sum)
  sum.bygene.bio<-apply(assay(Data)[1:q.bio,],2,sum)
  
  ls.phi0 = ifelse("ls.phi0" %in% names(args), args$ls.phi0, 11)
  
  # GENERATING STARTING VALUES
  if("Start" %in% names(args)) {Start = args$Start}
  else{Start=HiddenBASiCS_MCMC_Start(Data)}
  # Starting values for MCMC chains
  mu0=as.vector(Start$mu0); delta0=as.vector(Start$delta0)
  phi0=as.vector(Start$phi0); s0=as.vector(Start$s0)
  nu0=as.vector(Start$nu0); theta0=as.numeric(Start$theta0)
  # Starting values for adaptive proposal variances
  ls.mu0=as.vector(Start$ls.mu0); ls.delta0=as.vector(Start$ls.delta0)
  ls.phi0=as.numeric(Start$ls.phi0)
  ls.nu0=as.vector(Start$ls.nu0); ls.theta0=as.numeric(Start$ls.theta0)
  
  StoreAdaptNumber = as.numeric(StoreAdapt)
  nBatch = length(unique(metadata(Data)$BatchInfo))
  
  # If spikes are available
  if(length(metadata(Data)$SpikeInput) > 1)
  {
    if(nBatch > 1)
    {
      BatchDesign = model.matrix(~as.factor(metadata(Data)$BatchInfo)-1)
      
      # MCMC SAMPLER (FUNCTION IMPLEMENTED IN C++)
      Time = system.time(Chain <- HiddenBASiCS_MCMCcppBatch(
        N,
        Thin,
        Burn,
        as.matrix(assay(Data)),
        BatchDesign,
        mu0, delta0, phi0, s0, nu0, theta0,
        PriorParam$s2.mu,
        PriorParam$a.delta, PriorParam$b.delta,
        PriorParam$p.phi,
        PriorParam$a.s, PriorParam$b.s,
        PriorParam$a.theta, PriorParam$b.theta,
        AR,
        ls.mu0, ls.delta0, ls.phi0, ls.nu0, ls.theta0,
        sum.bycell.all, sum.bycell.bio, sum.bygene.all, sum.bygene.bio,
        StoreAdaptNumber,StopAdapt,as.numeric(PrintProgress),
        PriorParam$s2.delta, PriorDeltaNum))
    }
    else
    {
      # MCMC SAMPLER (FUNCTION IMPLEMENTED IN C++)
      Time = system.time(Chain <- HiddenBASiCS_MCMCcpp(
        N,
        Thin,
        Burn,
        as.matrix(assay(Data)),
        mu0, delta0, phi0, s0, nu0, theta0,
        PriorParam$s2.mu,
        PriorParam$a.delta, PriorParam$b.delta,
        PriorParam$p.phi,
        PriorParam$a.s, PriorParam$b.s,
        PriorParam$a.theta, PriorParam$b.theta,
        AR,
        ls.mu0, ls.delta0, ls.phi0, ls.nu0, ls.theta0,
        sum.bycell.all, sum.bycell.bio, sum.bygene.all, sum.bygene.bio,
        StoreAdaptNumber,StopAdapt,as.numeric(PrintProgress),
        PriorParam$s2.delta, PriorDeltaNum))
    }  
  }
  # If spikes are not available
  else
  {
    message("--------------------------------------------------------------------", "\n",
            'IMPORTANT: this part of the code is under development. DO NOT USE. \n',
            "--------------------------------------------------------------------", "\n")
    
    #if(PriorDelta == "gamma") stop("PriorDelta = 'gamma' is not supported for the no-spikes case")
    
    # 1: Full constrain; 2: Non-zero genes only
    ConstrainType = ifelse("ConstrainType" %in% names(args),args$ConstrainType, 2)
    ConstrainLimit = ifelse("ConstrainLimit" %in% names(args),args$ConstrainLimit, 1)
    ConstrainAlpha = ifelse("ConstrainAlpha" %in% names(args),args$ConstrainAlpha, 0.05)
    ConstrainProb = ifelse("ConstrainProb" %in% names(args),args$ConstrainProb, 0.95)
    
    BatchDesign = model.matrix(~as.factor(metadata(Data)$BatchInfo)-1)
    BatchSizes = table(metadata(Data)$BatchInfo)
    BatchIds = as.numeric(names(BatchSizes))
    BatchOffSet = rep(1, times = nBatch)
    for(k in 2:nBatch)
    {
      BatchOffSet[k] = median(colSums(assay(Data)[,metadata(Data)$BatchInfo == BatchIds[k]])) / 
        median(colSums(assay(Data)[,metadata(Data)$BatchInfo == BatchIds[1]]))
    }
    # Auxiliary vector contaning a gene index
    Index = (1:q.bio) - 1
    # In the following '+1' is used as c++ vector indexes vectors setting '0' as its first element
    # Constrain for gene-specific expression rates
    if(ConstrainType == 1) # Full constrain
    {
      # Note we use 'ConstrainLimit + 1' as 1 pseudo-count was added when computing 'mu0' (to avoid numerical issues)
      ConstrainGene = (1:q.bio) - 1
      NotConstrainGene = 0
      Constrain = mean(log(mu0[ConstrainGene+1]))
    }
    if(ConstrainType == 2) # Trimmed constrain based on mean
    {
      # Note we use 'ConstrainLimit + 1' as 1 pseudo-count was added when computing 'mu0' (to avoid numerical issues)
      ConstrainGene = which(mu0 >= ConstrainLimit + 1) - 1
      NotConstrainGene = which(mu0 < ConstrainLimit + 1) - 1
      Constrain = mean(log(mu0[ConstrainGene+1]))
    }
    if(ConstrainType == 3) # Trimmed constrain based on detection
    {
      Detection = rowMeans(assay(Data) > 0)
      ConstrainGene = which(Detection >= ConstrainLimit) - 1
      NotConstrainGene = which(Detection < ConstrainLimit) - 1
      Constrain = mean(log(mu0[ConstrainGene+1]))
    }
    
    StochasticRef = ifelse("StochasticRef" %in% names(args),args$StochasticRef, FALSE)
    
    if(StochasticRef == TRUE)
    {
      aux.ref = cbind(ConstrainGene, abs(log(mu0[ConstrainGene+1]) - Constrain))
      aux.ref = aux.ref[order(aux.ref[,2]),]
      RefGenes = aux.ref[1:200, 1]
      RefGene = RefGenes[1]
    }
    else
    {
      aux.ref = which(abs(log(mu0[ConstrainGene+1]) - Constrain) == min(abs(log(mu0[ConstrainGene+1]) - Constrain)))[1]
      RefGene = ConstrainGene[aux.ref]   
      RefGenes = RefGene    
    }
    
    # MCMC SAMPLER (FUNCTION IMPLEMENTED IN C++)
    Time = system.time(Chain <- HiddenBASiCS_MCMCcppNoSpikes(
      N,
      Thin,
      Burn,
      as.matrix(assay(Data)),
      BatchDesign,
      mu0, delta0, phi0, nu0, theta0,
      PriorParam$s2.mu,
      PriorParam$a.delta, PriorParam$b.delta,
      PriorParam$a.phi, PriorParam$b.phi,
      PriorParam$a.theta, PriorParam$b.theta,
      AR,
      ls.mu0, ls.delta0, ls.nu0, ls.theta0,
      sum.bycell.all, sum.bygene.all, 
      StoreAdaptNumber,StopAdapt,as.numeric(PrintProgress),
      PriorParam$s2.delta, PriorDeltaNum, 
      metadata(Data)$BatchInfo, BatchIds, as.vector(BatchSizes), BatchOffSet,
      Constrain, Index, RefGene, RefGenes, ConstrainGene, NotConstrainGene, 
      ConstrainType))
  }
  
  
  Chain$mu = Chain$mu[,1:q.bio]
  colnames(Chain$mu) = rownames(assay(Data))[!rowData(Data)$Tech]
  colnames(Chain$delta) = rownames(assay(Data))[!rowData(Data)$Tech]
  colnames(Chain$phi) = paste0("Cell",1:n)
  if(length(metadata(Data)$SpikeInput) > 1) {colnames(Chain$s) = paste0("Cell",1:n)}
  colnames(Chain$nu) = paste0("Cell",1:n)
  colnames(Chain$theta) = paste0("Batch",1:nBatch)
  
  cat("--------------------------------------------------------------------", "\n")
  cat("MCMC running time", "\n")
  cat("--------------------------------------------------------------------", "\n")
  cat(Time)
  cat("\n")
  
  OldDir = getwd()
  
  if(StoreChains)
  {
    setwd(StoreDir)
    
    message("--------------------------------------------------------------------", "\n",
            "Storing MCMC chains of model parameters as .txt files in", "\n",
            paste0("'",StoreDir,"' directory ... "), "\n",
            "--------------------------------------------------------------------", "\n")
    
    write.table(Chain$mu[,1:q.bio],paste0("chain_mu_",RunName,".txt"),col.names=TRUE,row.names=FALSE)
    write.table(Chain$delta,paste0("chain_delta_",RunName,".txt"),col.names=TRUE,row.names=FALSE)
    write.table(Chain$phi,paste0("chain_phi_",RunName,".txt"),col.names=TRUE,row.names=FALSE)
    if(length(metadata(Data)$SpikeInput) > 1){write.table(Chain$s,paste0("chain_s_",RunName,".txt"),col.names=TRUE,row.names=FALSE)}
    write.table(Chain$nu,paste0("chain_nu_",RunName,".txt"),col.names=TRUE,row.names=FALSE)
    write.table(Chain$theta,paste0("chain_theta_",RunName,".txt"),col.names=TRUE,row.names=FALSE)
    
    setwd(OldDir)
  }
  
  if(StoreAdapt)
  {
    setwd(StoreDir)
    
    message("--------------------------------------------------------------------", "\n",
            "Storing trajectories of adaptive proposal variances (log-scale) as .txt files in", "\n",
            paste0("'",StoreDir,"' directory ... "), "\n",
            "--------------------------------------------------------------------", "\n")
    
    colnames(Chain$ls.mu) = rownames(assay(Data))[!rowData(Data)$Tech]
    colnames(Chain$ls.delta) = rownames(assay(Data))[!rowData(Data)$Tech]
    colnames(Chain$ls.phi) = "AllCells"
    colnames(Chain$ls.nu) = paste0("Cell",1:n)
    colnames(Chain$ls.theta) = paste0("Batch",1:nBatch)
    
    write.table(Chain$ls.mu,paste0("chain_ls.mu_",RunName,".txt"),col.names=T,row.names=F)
    write.table(Chain$ls.delta,paste0("chain_ls.delta_",RunName,".txt"),col.names=T,row.names=F)
    write.table(Chain$ls.phi,paste0("chain_ls.phi_",RunName,".txt"),col.names=T,row.names=F)
    write.table(Chain$ls.nu,paste0("chain_ls.nu_",RunName,".txt"),col.names=T,row.names=F)
    write.table(Chain$ls.theta,paste0("chain_ls.theta_",RunName,".txt"),col.names=T,row.names=F)
    
    setwd(OldDir)
  }
  
  message("--------------------------------------------------------------------", "\n",
          "Output", "\n",
          "--------------------------------------------------------------------", "\n")
  
  if(length(metadata(Data)$SpikeInput) == 1) {Chain$s <- matrix(1, ncol = ncol(Chain$phi), nrow = nrow(Chain$phi))}
  
  ChainClass <- newBASiCS_Chain(mu = Chain$mu, delta = Chain$delta, phi = Chain$phi,
                                s = Chain$s, nu = Chain$nu, theta = Chain$theta)
  
  if(length(metadata(Data)$SpikeInput) == 1)
  {
    message("\n",
            "--------------------------------------------------------------------", "\n", 
            paste("BASiCS version", packageVersion("BASiCS"), ": horizontal integration (no-spikes case)"), "\n",
            "--------------------------------------------------------------------", "\n",
            paste("ConstrainType:", ConstrainType), "\n")  
    if(length(RefGenes) == 1) 
    {
      message(paste("Reference gene:", RefGene + 1), "\n",
              paste("Information stored as a .txt file in"), "\n",
              paste0("'",StoreDir,"' directory ... "), "\n",
              "--------------------------------------------------------------------", "\n")
      
      setwd(StoreDir)
      
      TableRef = cbind.data.frame("GeneNames" = rownames(assay(Data))[RefGene+1], 
                                  "GeneIndex" = RefGene+1, 
                                  stringsAsFactors = FALSE)
      write.table(TableRef,paste0("TableRef_",RunName,".txt"), col.names = T, row.names = F)
      
      setwd(OldDir)      
    }
    else
    {
      setwd(StoreDir)
      
      TableRef = cbind.data.frame("GeneNames" = rownames(assay(Data))[RefGenes+1], 
                                  "GeneIndex" = RefGenes+1, 
                                  "ReferenceFreq" = Chain$RefFreq[RefGenes+1],
                                  stringsAsFactors = FALSE)
      write.table(TableRef,paste0("TableRef_",RunName,".txt"), col.names = T, row.names = F)
      
      setwd(OldDir)
      
      message(paste("Randomly, 1 out of", length(RefGenes), "genes was left as reference at each iteration"), "\n",
              paste("List of reference genes and their associated frequencies stored as a .txt file in"), "\n",
              paste0("'",StoreDir,"' directory ... "), "\n",
              "--------------------------------------------------------------------", "\n")
    }
  }
  else
  {
    message("--------------------------------------------------------------------", "\n",
            paste("BASiCS version", packageVersion("BASiCS"), ": vertical integration (spikes case)"), "\n",
            "--------------------------------------------------------------------", "\n")
  }
  
  return(ChainClass)
}
