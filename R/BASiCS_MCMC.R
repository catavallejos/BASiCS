#' @title BASiCS MCMC sampler
#'
#' @description MCMC sampler to perform Bayesian inference for single-cell 
#' mRNA sequencing datasets using the model described in Vallejos et al (2015).
#'
#' @param Data A \code{\linkS4class{SingleCellExperiment}} object. 
#' This MUST be formatted to include the spike-ins information (see vignette). 
#' @param N Total number of iterations for the MCMC sampler. 
#' Use \code{N>=max(4,Thin)}, \code{N} being a multiple of \code{Thin}.
#' @param Thin Thining period for the MCMC sampler. Use \code{Thin>=2}.
#' @param Burn Burn-in period for the MCMC sampler. Use \code{Burn>=1}, 
#' \code{Burn<N}, \code{Burn} being a multiple of \code{Thin}.
#' @param ... Optional parameters.
#' \describe{
#' \item{\code{PriorDelta}}{Specifies the prior used for \code{delta}. 
#' Possible values are 'gamma' (Gamma(\code{a.theta},\code{b.theta}) prior) and 
#' 'log-normal' (log-Normal(\code{0},\code{s2.delta}) prior) .}. 
#' Default value: \code{PriorDelta = 'log-normal'}. 
#' \item{\code{PriorParam}}{List of 7 elements, containing the hyper-parameter 
#' values required for the adopted prior (see Vallejos et al, 2015, 2016). 
#' All elements must be positive real numbers.
#' \describe{
#'   \item{\code{s2.mu}}{Scale hyper-parameter for the 
#'         log-Normal(\code{0},\code{s2.mu}) prior that is shared by all 
#'         gene-specific expression rate parameters \eqn{\mu_i}.
#'         Default: \code{s2.mu = 0.5}.}
#'   \item{\code{s2.delta}}{Only used when `PriorDelta == 'log-normal'`. 
#'         Scale hyper-parameter for the log-Normal(\code{0},\code{s2.delta}) 
#'         prior that is shared by all gene-specific over-dispersion parameters 
#'         \eqn{\delta_i}. Default: \code{s2.delta = 0.5}. }
#'   \item{\code{a.delta}}{Only used when `PriorDelta == 'gamma'`. 
#'         Shape hyper-parameter for the Gamma(\code{a.delta},\code{b.delta}) 
#'         prior that is shared by all gene-specific biological over-dispersion 
#'         parameters \eqn{\delta_i}. Default: \code{a.delta = 1}.}
#'   \item{\code{b.delta}}{Only used when `PriorDelta == 'gamma'`. 
#'         Rate hyper-parameter for the Gamma(\code{a.delta},\code{b.delta}) 
#'         prior that is shared by all gene-specific biological over-dispersion 
#'         hyper-parameters \eqn{\delta_i}. Default: \code{b.delta = 1}.}
#'   \item{\code{p.phi}}{Dirichlet hyper-parameter for the joint of all 
#'         (scaled by \code{n}) cell-specific mRNA content normalising constants 
#'         \eqn{\phi_j / n}. Default: \code{p.phi} \code{= rep(1, n)}.}
#'   \item{\code{a.s}}{Shape hyper-parameter for the 
#'         Gamma(\code{a.s},\code{b.s}) prior that is shared by all 
#'         cell-specific capture efficiency normalising constants \eqn{s_j}.
#'         Default: \code{a.s = 1}.}
#'   \item{\code{b.s}}{Rate hyper-parameter for the Gamma(\code{a.s},\code{b.s}) 
#'         prior that is shared by all cell-specific capture efficiency 
#'         normalising constants \eqn{s_j}. Default: \code{b.s = 1}.}
#'   \item{\code{a.theta}}{Shape hyper-parameter for the 
#'         Gamma(\code{a.theta},\code{b.theta}) prior for technical noise 
#'         parameter \eqn{\theta}. Default: \code{a.theta = 1}.}
#'   \item{\code{b.theta}}{Rate hyper-parameter for the 
#'         Gamma(\code{a.theta},\code{b.theta}) prior for technical noise 
#'         parameter \eqn{\theta}. Default: \code{b.theta = 1}.}
#'
#' }}
#' \item{\code{AR}}{Optimal acceptance rate for adaptive Metropolis Hastings 
#'       updates. It must be a positive number between 0 and 1. Default 
#'       (and recommended): \code{AR = 0.44}}.
#'
#' \item{\code{StopAdapt}}{Iteration at which adaptive proposals are not longer 
#'       adapted. Use \code{StopAdapt>=1}. Default: \code{StopAdapt = Burn}.}
#'
#' \item{\code{StoreChains}}{If \code{StoreChains = TRUE}, the generated 
#'       \code{BASiCS_Chain} object is stored as a `.Rds` file (\code{RunName} 
#'       argument used to index the file name). 
#'       Default: \code{StoreChains = FALSE}.}
#' \item{\code{StoreAdapt}}{If \code{StoreAdapt = TRUE}, trajectory of 
#'       adaptive proposal variances (in log-scale) for all parameters is 
#'       stored as a list in a `.Rds` file (\code{RunName} argument used to 
#'       index file name). Default: \code{StoreAdapt = FALSE}.}
#' \item{\code{StoreDir}}{Directory where output files are stored. 
#'       Only required if \code{StoreChains = TRUE} and/or 
#'       \code{StoreAdapt = TRUE}). Default: \code{StoreDir = getwd()}.}
#' \item{\code{RunName}}{String used to index `.Rds` files storing chains 
#'       and/or adaptive proposal variances.}
#' \item{\code{PrintProgress}}{If \code{PrintProgress = FALSE}, console-based 
#'       progress report is suppressed.}
#' \item{\code{ls.phi0}}{Starting value for the adaptive concentration parameter 
#'       of the Metropolis proposals for 
#'       \eqn{\phi = (\phi_1, \ldots, \phi_n)'}.}
#' \item{\code{Start}}{Starting values for the MCMC sampler. We do not advise to 
#'       specify this argument. Default options have been tuned to facilitate 
#'       convergence. If changed, it must be a list containing the following 
#'       elements: \code{mu0}, \code{delta0}, \code{phi0}, \code{s0}, 
#'       \code{nu0}, \code{theta0}, \code{ls.mu0}, \code{ls.delta0}, 
#'       \code{ls.phi0}, \code{ls.nu0} and \code{ls.theta0}}
#' }
#'
#' @return An object of class \code{\link[BASiCS]{BASiCS_Chain}}.
#'
#' @examples
#'
#' # Built-in simulated dataset
#' Data = makeExampleBASiCS_Data()
#' # To analyse real data, please refer to the instructions in: 
#' # https://github.com/catavallejos/BASiCS/wiki/2.-Input-preparation
#'
#' # Only a short run of the MCMC algorithm for illustration purposes
#' # Longer runs migth be required to reach convergence
#' Chain <- BASiCS_MCMC(Data, N = 50, Thin = 2, Burn = 10, 
#'                      PrintProgress = FALSE)
#' 
#' # For illustration purposes we load a built-in 'BASiCS_Chain' object 
#' # (obtained using the 'BASiCS_MCMC' function)
#' data(ChainSC)
#' 
#' # `displayChainBASiCS` can be used to extract information from this output. 
#' # For example:
#' head(displayChainBASiCS(ChainSC, Param = 'mu'))
#'
#' # Traceplot (examples only)
#' plot(ChainSC, Param = 'mu', Gene = 1)
#' plot(ChainSC, Param = 'phi', Cell = 1)
#' plot(ChainSC, Param = 'theta', Batch = 1)
#'
#' # Calculating posterior medians and 95% HPD intervals
#' ChainSummary <- Summary(ChainSC)
#' 
#' # `displaySummaryBASiCS` can be used to extract information from this output. 
#' # For example:
#' head(displaySummaryBASiCS(ChainSummary, Param = 'mu'))
#'
#' # Graphical display of posterior medians and 95% HPD intervals 
#' # For example:
#' plot(ChainSummary, Param = 'mu', main = 'All genes')
#' plot(ChainSummary, Param = 'mu', Genes = 1:10, main = 'First 10 genes')
#' plot(ChainSummary, Param = 'phi', main = 'All cells')
#' plot(ChainSummary, Param = 'phi', Cells = 1:5, main = 'First 5 cells')
#' plot(ChainSummary, Param = 'theta')
#'
#' # To constrast posterior medians of cell-specific parameters 
#' # For example:
#' par(mfrow = c(1,2))
#' plot(ChainSummary, Param = 'phi', Param2 = 's', SmoothPlot = FALSE)
#' # Recommended for large numbers of cells
#' plot(ChainSummary, Param = 'phi', Param2 = 's', SmoothPlot = TRUE) 
#'
#' # To constrast posterior medians of gene-specific parameters
#' par(mfrow = c(1,2))
#' plot(ChainSummary, Param = 'mu', Param2 = 'delta', log = 'x', 
#'      SmoothPlot = FALSE)
#' # Recommended
#' plot(ChainSummary, Param = 'mu', Param2 = 'delta', log = 'x', 
#'      SmoothPlot = TRUE) 
#'
#' # Highly and lowly variable genes detection (within a single group of cells)
#' DetectHVG <- BASiCS_DetectHVG(ChainSC, VarThreshold = 0.60, 
#'                               EFDR = 0.10, Plot = TRUE)
#' DetectLVG <- BASiCS_DetectLVG(ChainSC, VarThreshold = 0.40, 
#'                               EFDR = 0.10, Plot = TRUE)
#'
#' plot(ChainSummary, Param = 'mu', Param2 = 'delta', log = 'x', col = 8)
#' with(DetectHVG$Table, points(Mu[HVG == TRUE], Delta[HVG == TRUE],
#'        pch = 16, col = 'red', cex = 1))
#' with(DetectLVG$Table, points(Mu[LVG == TRUE], Delta[LVG == TRUE],
#'        pch = 16, col = 'blue', cex = 1))
#'
#' # If variance thresholds are not fixed
#' BASiCS_VarThresholdSearchHVG(ChainSC, 
#'                              VarThresholdsGrid = seq(0.55,0.65,by=0.01), 
#'                              EFDR = 0.10)
#' BASiCS_VarThresholdSearchLVG(ChainSC, 
#'                              VarThresholdsGrid = seq(0.35,0.45,by=0.01), 
#'                              EFDR = 0.10)
#'                              
#' # To obtain denoised rates / counts, see:
#' help(BASiCS_DenoisedRates)
#' help(BASiCS_DenoisedCounts)
#' 
#' # For examples of differential analyses between 2 populations of cells see:
#' help(BASiCS_TestDE)
#' 
#'
#' @author Catalina A. Vallejos \email{cnvallej@@uc.cl} and Nils Eling
#'
#' @references 
#' Vallejos, Marioni and Richardson (2015). PLoS Computational Biology. 
#' 
#' Vallejos, Richardson and Marioni (2016). Genome Biology.
BASiCS_MCMC <- function(Data, N, Thin, Burn, ...) {
    
  if (!is(Data, "SingleCellExperiment")) 
    stop("'Data' is not a SingleCellExperiment class object.")
  # The following checks are only relevant when the input data was
  # not created using the `newBASiCS_Data` function
  if(!("SpikeInput" %in% names(metadata(Data))))
    stop("'Data' does not contained all the required information \n",
         "See: https://github.com/catavallejos/BASiCS/wiki/2.-Input-preparation")
  if(!("BatchInfo" %in% names(metadata(Data))))
    stop("'Data' does not contained all the required information \n",
         "See: https://github.com/catavallejos/BASiCS/wiki/2.-Input-preparation")
  errors <- HiddenChecksBASiCS_Data(assay(Data), isSpike(Data), 
                                    metadata(Data)$SpikeInput, 
                                    rownames(assay(Data)), 
                                    metadata(Data)$BatchInfo)
  if (length(errors) > 0) stop(errors) 
    
  # SOME QUANTITIES USED THROUGHOUT THE MCMC ALGORITHM
  q <- length(isSpike(Data))
  q.bio <- sum(!isSpike(Data))
  n <- dim(assay(Data))[2]
    
  args <- list(...)
    
  if (!("PriorDelta" %in% names(args)) & !("Regression" %in% names(args))) 
  {
    message("-------------------------------------------------------------\n", 
    "NOTE: default choice PriorDelta = 'log-normal'  (recommended value). \n",
    "Vallejos et al (2015) used a 'gamma' prior instead.\n", 
            "-------------------------------------------------------------\n")
  }
    
  if ("PriorParam" %in% names(args)) { PriorParam <- args$PriorParam } 
  else 
  {
    PriorParam <- list(s2.mu = 0.5, s2.delta = 0.5, a.delta = 1, 
                       b.delta = 1, p.phi = rep(1, times = n), 
                       a.phi = 1, b.phi = 1, a.s = 1, b.s = 1, 
                       a.theta = 1, b.theta = 1)
  }
  AR <- ifelse("AR" %in% names(args), args$AR, 0.44)
  StopAdapt <- ifelse("StopAdapt" %in% names(args), args$StopAdapt, Burn)
  StoreChains <- ifelse("StoreChains" %in% names(args), args$StoreChains, FALSE)
  StoreAdapt <- ifelse("StoreAdapt" %in% names(args), args$StoreAdapt, FALSE)
  StoreDir <- ifelse("StoreDir" %in% names(args), args$StoreDir, getwd())
  RunName <- ifelse("RunName" %in% names(args), args$RunName, "")
  PrintProgress <- ifelse("PrintProgress" %in% names(args), 
                          args$PrintProgress, TRUE)
  PriorDelta <- ifelse("PriorDelta" %in% names(args), 
                       args$PriorDelta, "log-normal")
  CollapsedMCMC <- ifelse("CollapsedMCMC" %in% names(args), args$CollapsedMCMC, FALSE)
  
  # This is specific for the regression case
  if("Regression" %in% names(args)){
    if(args$Regression == TRUE){
      k <- ifelse("k" %in% names(args), args$k, 12)
      if(k <= 3) stop("The number of basis functions needs to be >= 4.")
    }
  }
  variance <- ifelse("Var" %in% names(args), args$Var, 1.2)
  eta <- ifelse("eta" %in% names(args), args$eta, 5)
    
  if (!(length(N) == 1 | length(Thin) == 1 | length(Burn) == 1)) 
    stop("Invalid parameter values.")
  if (!(N%%Thin == 0 & N >= max(4, Thin))) 
    stop("Please use an integer value for N (N>=4); multiple of thin.")
  if (!(Thin%%1 == 0 & Thin >= 2)) 
    stop("Please use an integer value for Thin (Thin>=2).")
  if (!(Burn%%Thin == 0 & Burn < N & Burn >= 1)) 
    stop("Please use an integer value for Burn (1<=Burn<N); multiple of thin.")
    
  if (!(PriorParam$s2.mu > 0 & length(PriorParam$s2.mu) == 1 & 
        PriorParam$s2.delta > 0 & length(PriorParam$s2.delta) == 1 & 
        PriorParam$a.delta > 0 & length(PriorParam$a.delta) == 1 & 
        PriorParam$b.delta > 0 & length(PriorParam$b.delta) == 1 & 
        all(PriorParam$p.phi > 0) & length(PriorParam$p.phi) == n & 
        PriorParam$a.s > 0 & length(PriorParam$a.s) == 1 & 
        PriorParam$b.s > 0 & length(PriorParam$b.s) == 1 & 
        PriorParam$a.theta > 0 & length(PriorParam$a.theta) == 1 & 
        PriorParam$b.theta > 0) & length(PriorParam$b.theta) == 1) 
    stop("Invalid prior hyper-parameter values.")
    
  if (!(AR > 0 & AR < 1 & length(AR) == 1)) 
    stop("Invalid AR value. Recommended value: AR = 0.44.")
  if (!(StopAdapt > 0)) 
    stop("Invalid StopAdapt value.")
  if (!(is.logical(StoreChains) & length(StoreChains) == 1)) 
    stop("Invalid StoreChains value.")
  if (!(is.logical(StoreAdapt) & length(StoreAdapt) == 1)) 
    stop("Invalid StoreAdapt value.")
  if (!(file.info(StoreDir)["isdir"])) 
    stop("Invalid StoreDir value.")
  if (!(PriorDelta %in% c("gamma", "log-normal"))) 
    stop("Invalid PriorDelta value.")
    
  PriorDeltaNum <- ifelse(PriorDelta == "gamma", 1, 2)
    
  # SOME SUMS USED THROUGHOUT THE MCMC ALGORITHM
  sum.bycell.all <- matrixStats::rowSums2(assay(Data))
  sum.bycell.bio <- matrixStats::rowSums2(assay(Data)[1:q.bio, ])
  sum.bygene.all <- matrixStats::colSums2(assay(Data))
  sum.bygene.bio <- matrixStats::colSums2(assay(Data)[1:q.bio, ])
    
  ls.phi0 <- ifelse("ls.phi0" %in% names(args), args$ls.phi0, 11)
    
  # GENERATING STARTING VALUES
  if ("Start" %in% names(args)) { Start = args$Start } 
  else
  {
    if("Regression" %in% names(args))
    {
      if(args$Regression == TRUE)
      {
        Start <- HiddenBASiCS_MCMC_Start(Data, k = k, eta = eta)
      }
    }
    else { Start <- HiddenBASiCS_MCMC_Start(Data) }
  }
  
  # Starting values for MCMC chains
  mu0 <- as.vector(Start$mu0)
  delta0 <- as.vector(Start$delta0)
  phi0 <- as.vector(Start$phi0)
  s0 <- as.vector(Start$s0)
  nu0 <- as.vector(Start$nu0)
  theta0 <- as.numeric(Start$theta0)
  
  # Starting values for adaptive proposal variances
  ls.mu0 <- as.vector(Start$ls.mu0)
  ls.delta0 <- as.vector(Start$ls.delta0)
  ls.phi0 <- as.numeric(Start$ls.phi0)
  ls.nu0 <- as.vector(Start$ls.nu0)
  ls.theta0 <- as.numeric(Start$ls.theta0)
  
  # Starting values for Regression
  if("Regression" %in% names(args))
  {
    if(args$Regression == TRUE)
    {
      m0 <- as.vector(Start$m0); V0 <- as.matrix(Start$V0)
      sigma2.a0 <- Start$sigma2.a0; sigma2.b0 <- Start$sigma2.b0
      beta0 <- Start$beta0; sigma20 <- Start$sigma20
      lambda0 <- Start$lambda0
    }
  }
    
  StoreAdaptNumber <- as.numeric(StoreAdapt)
  nBatch <- length(unique(metadata(Data)$BatchInfo))
    
  if(nBatch > 1)
  {
    BatchDesign <- model.matrix(~as.factor(metadata(Data)$BatchInfo) - 1)  
    BatchInfo <- as.numeric(metadata(Data)$BatchInfo)
  }
  else
  { 
    # If there are no batches or the user asked to ignore them
    BatchDesign <- matrix(1, nrow = n, ncol = 1) 
    BatchInfo <- rep(1, times = n)
    nBatch <- 1
  }
    
  # If spikes are available (stable version)
  if (length(metadata(Data)$SpikeInput) > 1) 
  {
    # If regression case is chosen
    if("Regression" %in% names(args)){
      if(args$Regression == TRUE){
      Time = system.time(Chain <- HiddenBASiCS_MCMCcppReg(N, Thin, Burn, 
                                                           as.matrix(assay(Data)), 
                                                           BatchDesign,
                                                           mu0[(q.bio+1):q],
                                                           mu0[1:q.bio], delta0, 
                                                           phi0, s0, 
                                                           nu0, rep(theta0, nBatch), 
                                                           PriorParam$s2.mu, 
                                                           PriorParam$p.phi, 
                                                           PriorParam$a.s, 
                                                           PriorParam$b.s, 
                                                           PriorParam$a.theta, 
                                                           PriorParam$b.theta, 
                                                           AR, 
                                                           ls.mu0[1:q.bio], 
                                                           ls.delta0, 
                                                           ls.phi0, ls.nu0, 
                                                           rep(ls.theta0, nBatch), 
                                                           sum.bycell.all, 
                                                           sum.bycell.bio, 
                                                           sum.bygene.all, 
                                                           sum.bygene.bio, 
                                                           StoreAdaptNumber, 
                                                           StopAdapt, 
                                                           as.numeric(PrintProgress),
                                                           k, m0, V0,
                                                           sigma2.a0, sigma2.b0, 
                                                           beta0, sigma20, 
                                                           eta, lambda0, 
                                                           variance))
      
      # Remove epsilons for genes that are not expressed in at least 2 cells
      genes <- apply(assay(Data)[!isSpike(Data),], 1, function(n){length(which(n>0))}) > 1
      
      Chain$epsilon[,!genes] <- NA
      }
    }
    else
    {
      # MCMC SAMPLER (FUNCTION IMPLEMENTED IN C++)
      if(!CollapsedMCMC)
      {
        Time <- system.time(Chain <- HiddenBASiCS_MCMCcpp(N, Thin, Burn, 
                                                          as.matrix(assay(Data)), 
                                                          BatchDesign,
                                                          mu0[(q.bio+1):q],
                                                          mu0[1:q.bio], delta0, 
                                                          phi0, s0, 
                                                          nu0, rep(theta0, nBatch), 
                                                          PriorParam$s2.mu, 
                                                          PriorParam$a.delta, 
                                                          PriorParam$b.delta, 
                                                          PriorParam$s2.delta,
                                                          PriorDeltaNum,
                                                          PriorParam$p.phi, 
                                                          PriorParam$a.s, 
                                                          PriorParam$b.s, 
                                                          PriorParam$a.theta, 
                                                          PriorParam$b.theta, 
                                                          AR, 
                                                          ls.mu0[1:q.bio], 
                                                          ls.delta0, 
                                                          ls.phi0, ls.nu0, 
                                                          rep(ls.theta0, nBatch), 
                                                          sum.bycell.all, 
                                                          sum.bycell.bio, 
                                                          sum.bygene.all, 
                                                          sum.bygene.bio, 
                                                          StoreAdaptNumber, 
                                                          StopAdapt, 
                                                          as.numeric(PrintProgress)))       
      }
      else
      {
        Time <- system.time(Chain <- HiddenBASiCS_Cata_MCMCcpp(N, Thin, Burn, 
                                                          as.matrix(assay(Data)), 
                                                          BatchDesign,
                                                          mu0[(q.bio+1):q],
                                                          mu0[1:q.bio], delta0, 
                                                          phi0, 
                                                          nu0, rep(theta0, nBatch), 
                                                          PriorParam$s2.mu, 
                                                          PriorParam$a.delta, 
                                                          PriorParam$b.delta, 
                                                          PriorParam$s2.delta,
                                                          PriorDeltaNum,
                                                          PriorParam$p.phi, 
                                                          PriorParam$a.s, 
                                                          PriorParam$b.s, 
                                                          PriorParam$a.theta, 
                                                          PriorParam$b.theta, 
                                                          AR, 
                                                          ls.mu0[1:q.bio], 
                                                          ls.delta0, 
                                                          ls.phi0, ls.nu0, 
                                                          rep(ls.theta0, nBatch), 
                                                          sum.bycell.all, 
                                                          sum.bycell.bio, 
                                                          sum.bygene.all, 
                                                          sum.bygene.bio, 
                                                          StoreAdaptNumber, 
                                                          StopAdapt, 
                                                          as.numeric(PrintProgress)))
      }

    }
  } 
  else 
  {
    # If spikes are not available
    message("-------------------------------------------------------------\n",  
            "IMPORTANT: this code is under development. DO NOT USE \n", 
            "This part of the code is just a place-holder \n", 
            "-------------------------------------------------------------\n")
        
    if (PriorDelta == "gamma") 
      stop("PriorDelta = 'gamma' is not supported for the no-spikes case")
        
    # 1: Full constrain; 2: Genes with average count >= 1
    ConstrainType <- ifelse("ConstrainType" %in% names(args), 
                            args$ConstrainType, 2)
    if (ConstrainType == 1) 
    {
      ConstrainGene <- (1:q.bio) - 1
      NotConstrainGene <- 0
    }
    if (ConstrainType == 2) 
    {
      ConstrainGene <- which(matrixStats::rowMeans2(assay(Data)) >= 1) - 1
      NotConstrainGene <- which(matrixStats::rowMeans2(assay(Data)) < 1) - 1
    }
    Constrain <- mean(log(mu0[ConstrainGene + 1]))
    
    # Whether or not a stochatic reference is used
    # If stochastic, range of possible reference values only includes 
    # the nearest 200 genes located around the constrain
    StochasticRef <- ifelse("StochasticRef" %in% names(args), 
                            args$StochasticRef, TRUE)
    if (StochasticRef == TRUE) 
    {
      aux.ref <- cbind(ConstrainGene, 
                       abs(log(mu0[ConstrainGene + 1]) - Constrain))
      aux.ref <- aux.ref[order(aux.ref[, 2]), ]
      # Fix for the code to run on the synthetic small dataset
      # generated by makeExample_BASiCS function (less than 200 genes)
      if(length(ConstrainGene) > 200) { RefGenes <- aux.ref[1:200, 1] }
      else { RefGenes <- aux.ref[, 1] }
      RefGene <- RefGenes[1]
    } 
    else 
    {
      aux.ref <- which(abs(log(mu0[ConstrainGene + 1]) - Constrain) == 
                         min(abs(log(mu0[ConstrainGene + 1]) - Constrain)))[1]
      RefGene <- ConstrainGene[aux.ref]
      RefGenes <- RefGene
    }

    # Auxiliary vector contaning a gene index
    Index <- (1:q.bio) - 1
        
    # MCMC SAMPLER (FUNCTION IMPLEMENTED IN C++)
    Time <- system.time(Chain <- HiddenBASiCS_MCMCcppNoSpikes(
                                                     N, Thin, Burn, 
                                                     as.matrix(assay(Data)), 
                                                     BatchDesign, 
                                                     mu0, delta0, 
                                                     phi0, nu0, 
                                                     rep(theta0, nBatch), 
                                                     PriorParam$s2.mu, 
                                                     PriorParam$a.delta, 
                                                     PriorParam$b.delta, 
                                                     PriorParam$s2.delta, 
                                                     PriorDeltaNum, 
                                                     PriorParam$a.phi, 
                                                     PriorParam$b.phi, 
                                                     PriorParam$a.theta, 
                                                     PriorParam$b.theta, 
                                                     AR, 
                                                     ls.mu0, ls.delta0, 
                                                     ls.nu0, 
                                                     rep(ls.theta0, nBatch), 
                                                     sum.bycell.all, 
                                                     sum.bygene.all, 
                                                     StoreAdaptNumber, 
                                                     StopAdapt, 
                                                     as.numeric(PrintProgress), 
                                                     Constrain, 
                                                     Index, 
                                                     RefGene, RefGenes, 
                                                     ConstrainGene, 
                                                     NotConstrainGene, 
                                                     ConstrainType,
                                                     as.numeric(StochasticRef)))
  }
    
  Chain$mu <- Chain$mu[, 1:q.bio]
  colnames(Chain$mu) <- rownames(assay(Data))[!isSpike(Data)]
  colnames(Chain$delta) <- rownames(assay(Data))[!isSpike(Data)]
  if(!is.null(Chain$epsilon) & !is.null(Chain$lambda)) 
    { colnames(Chain$epsilon) = colnames(Chain$lambda) = colnames(Chain$mu) }
  ###################### Change this to keep the cell labels the same
  CellLabels <- paste0(colnames(assay(Data)), "_Batch", 
                       metadata(Data)$BatchInfo)
  colnames(Chain$phi) <- CellLabels
  if (length(metadata(Data)$SpikeInput) > 1) { colnames(Chain$s) <- CellLabels }
  colnames(Chain$nu) <- CellLabels
  colnames(Chain$theta) <- paste0("Batch", unique(metadata(Data)$BatchInfo))
    
  message("-------------------------------------------------------------\n",
          "MCMC running time \n",
          "-------------------------------------------------------------\n",
          "user: ", round(Time['user.self'], 3), "\n", 
          "system: ", round(Time['sys.self'], 3), "\n", 
          "elapsed: ", round(Time['elapsed'], 3), "\n")
  
  message("-------------------------------------------------------------\n",
          "Output \n",
          "-------------------------------------------------------------\n")

    if (length(metadata(Data)$SpikeInput) == 1) 
    {
        Chain$s <- matrix(1, ncol = ncol(Chain$phi), nrow = nrow(Chain$phi))
    }
    ChainClass <- newBASiCS_Chain(parameters = Chain)
    
    OldDir <- getwd()
    
    if (StoreChains) 
    {
      setwd(StoreDir)
        
      message("-------------------------------------------------------------\n", 
              "BASiCS_Chain object stored as ", 
              paste0("chain_", RunName, ".Rds"), " file in", "\n", 
              paste0("'", StoreDir, "' directory ... "), "\n", 
              "-------------------------------------------------------------\n")
        
      saveRDS(ChainClass, file = paste0("chain_", RunName, ".Rds"))
        
      setwd(OldDir)
    }
    
    if (StoreAdapt)
    {
      setwd(StoreDir)
        
      message("-------------------------------------------------------------\n", 
              "Storing trajectories of adaptive proposal variances (log-scale) as", 
              "chain_ls_", RunName, ".Rds file in \n", 
              "'", StoreDir, "' directory ... \n", 
              "-------------------------------------------------------------\n")
        
      ChainLS <- list(ls.mu = Chain$ls.mu, ls.delta = Chain$ls.delta, 
                      ls.phi = Chain$ls.phi, ls.nu = Chain$ls.nu, 
                      ls.theta = Chain$ls.theta)
      saveRDS(ChainLS, file = paste0("chain_ls_", RunName, ".Rds"))
        
      setwd(OldDir)
    }
    
    # This 'if' refers to the no-spikes case 
    # Still under development - ignore for now!
    if (length(metadata(Data)$SpikeInput) == 1) 
    {
      message("-------------------------------------------------------------\n", 
              paste("BASiCS version", packageVersion("BASiCS"), 
                    ": horizontal integration (no-spikes case)"), "\n", 
              "-------------------------------------------------------------\n", 
              "ConstrainType: ", ConstrainType, "\n")
      if (length(RefGenes) == 1) 
      {
        message("Reference gene:", RefGene + 1, "\n", 
                "Information stored as a .txt file in \n", 
                "'", StoreDir, "' directory ... \n", 
                "-------------------------------------------------------------\n")
            
        setwd(StoreDir)
            
        TableRef <- cbind.data.frame(GeneNames = rownames(assay(Data))[RefGene + 1], 
                                     GeneIndex = RefGene + 1, 
                                     stringsAsFactors = FALSE)
        write.table(TableRef, paste0("TableRef_", RunName, ".txt"), 
                    col.names = TRUE, row.names = FALSE)
            
        setwd(OldDir)
      } 
      else 
      {
        setwd(StoreDir)
            
        TableRef = cbind.data.frame(GeneNames = rownames(assay(Data))[RefGenes + 1], 
                                    GeneIndex = RefGenes + 1, 
                                    ReferenceFreq = Chain$RefFreq[RefGenes + 1], 
                                    stringsAsFactors = FALSE)
        write.table(TableRef, paste0("TableRef_", RunName, ".txt"), 
                    col.names = TRUE, row.names = FALSE)
            
        setwd(OldDir)
            
        message("Randomly, 1 out of ", length(RefGenes), "\n",
                "genes was left as reference at each iteration \n", 
                "List of reference genes and frequencies stored as a .txt file in\n", 
                "'", StoreDir, "' directory ... \n", 
                "-------------------------------------------------------------\n")
      }
    } 
    else 
    {
      message("-------------------------------------------------------------\n", 
              "BASiCS version ", packageVersion("BASiCS"), " : \n", 
              "vertical integration (spikes case) \n", 
              "-------------------------------------------------------------\n") 
    }
    
    return(ChainClass)
}
