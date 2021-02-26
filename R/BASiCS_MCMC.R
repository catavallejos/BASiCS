#' @title BASiCS MCMC sampler
#'
#' @description MCMC sampler to perform Bayesian inference for single-cell
#' mRNA sequencing datasets using the model described in Vallejos et al (2015).
#'
#' @param Data A \code{\linkS4class{SingleCellExperiment}} object.
#' If \code{WithSpikes = TRUE}, this MUST be formatted to include
#' the spike-ins and/or batch information (see vignette).
#' @param N Total number of iterations for the MCMC sampler.
#' Use \code{N>=max(4,Thin)}, \code{N} being a multiple of \code{Thin}.
#' @param Thin Thining period for the MCMC sampler. Use \code{Thin>=2}.
#' @param Burn Burn-in period for the MCMC sampler. Use \code{Burn>=1},
#' \code{Burn<N}, \code{Burn} being a multiple of \code{Thin}.
#' @param Regression  If \code{Regression = TRUE}, BASiCS exploits a joint prior
#' formulation for mean and over-dispersion parameters to estimate a measure of
#' residual over-dispersion is not confounded by mean expression. Recommended
#' setting is \code{Regression = TRUE}.
#' @param WithSpikes  If \code{WithSpikes = TRUE}, BASiCS will use reads from
#' added spike-ins to estimate technical variability. 
#' If \code{WithSpikess = FALSE},
#' BASiCS depends on replicated experiments (batches) to estimate
#' technical variability. In this case, please supply the BatchInfo vector
#' in \code{colData(Data)}. Default: \code{WithSpikes = TRUE}.
#' @param PriorParam List of prior parameters for BASiCS_MCMC.
#' Should be created using \code{\link{BASiCS_PriorParam}}.
#' @param SubsetBy Character value specifying whether a divide and
#' conquer inference strategy should be used. When this is set to \code{"gene"},
#' inference is performed on batches of genes separately, and when it is set to
#' \code{"cell"}, inference is performed on batches of cells separately.
#' Posterior distributions are combined using posterior interval estimation
#' (see Li et al., 2016).
#' @param NSubsets If \code{SubsetBy="gene"} or 
#' \code{SubsetBy="cell"}, \code{NSubsets} specifies the number of batches 
#' to create and perform divide and conquer inference with.
#' @param CombineMethod The method used to combine 
#' subposteriors if \code{SubsetBy} is set to \code{"gene"} or 
#' \code{"cell"}. Options are \code{"pie"} corresponding to
#' posterior interval estimation (see Li et al., 2016) or
#' \code{"consensus"} (see Scott et al., 2016).
#' Both of these methods use a form of weighted average to
#' combine subposterior draws into the final posterior.
#' @param Weighting The weighting method used in the weighted
#' average chosen using \code{CombineMethod}. Available
#' options are \code{"naive"} (unweighted), \code{"n_weight"}
#' (weights are chosen based on the size of each partition)
#' and \code{"inverse_variance"} (subposteriors are weighted based
#' on the inverse of the variance of the subposterior for each 
#' parameter).
#' @param Threads Integer specifying the number of threads to be used to 
#' parallelise parameter updates. Default value is the globally set
#' \code{"Ncpus"} option, or 1 if this option is not set.
#' @param BPPARAM A \code{\link{BiocParallelParam}} instance,
#' used for divide and conquer inference.
#' @param ... Optional parameters.
#' \describe{
#'   \item{
#'     \code{AR}
#'   }{
#'     Optimal acceptance rate for adaptive Metropolis Hastings
#'     updates. It must be a positive number between 0 and 1. Default
#'     (and recommended): \code{AR = 0.44}.
#'   }
#'   \item{
#'     \code{StopAdapt}
#'   }{
#'     Iteration at which adaptive proposals are not longer
#'     adapted. Use \code{StopAdapt>=1}. Default: \code{StopAdapt = Burn}.
#'   }
#'   \item{
#'     \code{StoreChains}
#'   }{
#'     If \code{StoreChains = TRUE}, the generated
#'     \code{BASiCS_Chain} object is stored as a `.Rds` file (\code{RunName}
#'     argument used to index the file name).
#'     Default: \code{StoreChains = FALSE}.}
#'   \item{
#'     \code{StoreAdapt}
#'   }{
#'     If \code{StoreAdapt = TRUE}, trajectory of
#'     adaptive proposal variances (in log-scale) for all parameters is
#'     stored as a list in a `.Rds` file (\code{RunName} argument used to
#'     index file name). Default: \code{StoreAdapt = FALSE}.
#'   }
#'   \item{
#'     \code{StoreDir}
#'   }{
#'     Directory where output files are stored.
#'     Only required if \code{StoreChains = TRUE} and/or
#'     \code{StoreAdapt = TRUE}. Default: \code{StoreDir = getwd()}.
#'   }
#'   \item{
#'     \code{RunName}
#'   }{
#'     String used to index `.Rds` files storing chains
#'     and/or adaptive proposal variances.
#'   }
#'   \item{
#'     \code{PrintProgress}
#'   }{
#'     If \code{PrintProgress = FALSE}, console-based
#'     progress report is suppressed.}
#'   \item{
#'     \code{Start}
#'   }{
#'     Starting values for the MCMC sampler. We do not advise
#'     to use this argument. Default options have been tuned to facilitate
#'     convergence. If changed, it must be a list containing the following
#'     elements: 
#'     \code{mu0}, \code{delta0}, \code{phi0}, \code{s0}, \code{nu0}, 
#'     \code{theta0}, \code{ls.mu0}, \code{ls.delta0}, \code{ls.phi0},
#'     \code{ls.nu0} and \code{ls.theta0}
#'   }
#'   \item{
#'     \code{GeneExponent/CellExponent}
#'   }{
#'     Exponents applied to the prior for MCMC updates. Intended for use only 
#'     when performing divide & conquer MCMC strategies.
#'   }
#' }
#' @return An object of class \code{\link[BASiCS]{BASiCS_Chain}}.
#'
#' @examples
#'
#' # Built-in simulated dataset
#' set.seed(1) 
#' Data <- makeExampleBASiCS_Data()
#' # To analyse real data, please refer to the instructions in:
#' # https://github.com/catavallejos/BASiCS/wiki/2.-Input-preparation
#'
#' # Only a short run of the MCMC algorithm for illustration purposes
#' # Longer runs migth be required to reach convergence
#' Chain <- BASiCS_MCMC(Data, N = 50, Thin = 2, Burn = 10, Regression = FALSE,
#'                      PrintProgress = FALSE, WithSpikes = TRUE)
#'
#' # To run the regression version of BASiCS, use:
#' Chain <- BASiCS_MCMC(Data, N = 50, Thin = 2, Burn = 10, Regression = TRUE,
#'                      PrintProgress = FALSE, WithSpikes = TRUE)
#'
#' # To run the non-spike version BASiCS requires the data to contain at least
#' # 2 batches:
#' set.seed(2)
#' Data <- makeExampleBASiCS_Data(WithBatch = TRUE)
#' Chain <- BASiCS_MCMC(Data, N = 50, Thin = 2, Burn = 10, Regression = TRUE,
#'                      PrintProgress = FALSE, WithSpikes = FALSE)
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
#' # `displaySummaryBASiCS` can be used to extract information from this output
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
#' # To obtain denoised rates / counts, see:
#' # help(BASiCS_DenoisedRates)
#' # and
#' # help(BASiCS_DenoisedCounts)
#'
#' # For examples of differential analyses between 2 populations of cells see:
#' # help(BASiCS_TestDE)
#'
#' @author Catalina A. Vallejos \email{cnvallej@@uc.cl}
#' @author Nils Eling \email{eling@@ebi.ac.uk}
#'
#' @references
#' Vallejos, Marioni and Richardson (2015). PLoS Computational Biology.
#'
#' Vallejos, Richardson and Marioni (2016). Genome Biology.
#'
#' Eling et al (2018). Cell Systems
#'
#' Simple, Scalable and Accurate Posterior Interval Estimation
#' Cheng Li and Sanvesh Srivastava and David B. Dunson
#' arXiv (2016)
#' 
#' Bayes and Big Data:  The Consensus Monte Carlo Algorithm
#' Steven L. Scott, Alexander W. Blocker, Fernando V. Bonassi,
#' Hugh A. Chipman, Edward I. George and Robert E. McCulloch
#' International Journal of Management Science and Engineering 
#' Management (2016)
#' @export
BASiCS_MCMC <- function(
    Data,
    N,
    Thin,
    Burn,
    Regression,
    WithSpikes = TRUE,
    PriorParam = BASiCS_PriorParam(Data, PriorMu = "EmpiricalBayes"),
    SubsetBy = c("none", "gene", "cell"),
    NSubsets = 1,
    CombineMethod = c("pie", "consensus"),
    Weighting = c("naive", "n_weight", "inverse_variance"),
    Threads = getOption("Ncpus", default = 1L),
    BPPARAM = BiocParallel::bpparam(),
    ...) {


  # Checks to ensure input arguments are valid
  .BASiCS_MCMC_InputCheck(Data, N, Thin, Burn, Regression, WithSpikes, Threads,
    NSubsets)
  SubsetBy <- match.arg(SubsetBy)

  if (SubsetBy != "none" && NSubsets != 1) {
    if (SubsetBy == "cell") {
      warning("Divide and conquer inference using cell-wise partitions is not recommended.")
    }
    CombineMethod <- match.arg(CombineMethod)
    Weighting <- match.arg(Weighting)
    Chains <- BASiCS_DivideAndConquer(
      Data,
      N = N,
      NSubsets = NSubsets,
      Thin = Thin,
      Burn = Burn,
      Regression = Regression,
      WithSpikes = WithSpikes,
      SubsetBy = SubsetBy,
      PriorParam = PriorParam,
      BPPARAM = BPPARAM,
      ...
    )
    Chain <- .combine_subposteriors(
      Chains,
      CombineMethod = CombineMethod,
      SubsetBy = SubsetBy,
      Weighting = Weighting,
      BPPARAM = BPPARAM
    )
    return(Chain)
  }

  # Some global values used throughout the MCMC algorithm and checks
  # Numbers of cells/genes/batches + count table + design matrix for batch
  GPar <- .BASiCS_MCMC_GlobalParams(Data)

  # Total counts per cell and/or gene
  sum.bycell.bio <- Matrix::rowSums(counts(Data))
  sum.bygene.bio <- Matrix::colSums(counts(Data))
  if (WithSpikes) {
    sum.bycell.tech <- Matrix::rowSums(assay(altExp(Data)))
    sum.bygene.tech <- Matrix::colSums(assay(altExp(Data))) 
    sum.bycell.all <- c(sum.bycell.bio, sum.bycell.tech)
    sum.bygene.all <- sum.bygene.bio + sum.bygene.tech
  } else {
    sum.bycell.all <- sum.bycell.bio
    sum.bygene.all <- sum.bygene.bio
  }

  # Assignment of optional parameters (default values if input not provided)
  ArgsDef <- .BASiCS_MCMC_ExtraArgs(
    Data = Data,
    Burn = Burn,
    GPar = GPar,
    Regression = Regression,
    WithSpikes = WithSpikes,
    PriorParam = PriorParam,
    ...
  )
  
  # Starting values for the MCMC (parameters and adaptive variances)
  # Loaded separately to simplify calling to the elements of its list
  # Same for prior parameters
  Start <- ArgsDef$Start
  PriorParam <- ArgsDef$PriorParam

  MCMCArgs <- list(
    N = N,
    Thin = Thin,
    Burn = Burn,
    Counts = GPar$BioCounts,
    BatchDesign = GPar$BatchDesign,
    mu0 = Start$mu0,
    delta0 = Start$delta0,
    s0 = Start$s0,
    nu0 = Start$nu0,
    theta0 = rep(Start$theta0, GPar$nBatch),
    mu_mu = PriorParam$mu.mu,
    s2mu = PriorParam$s2.mu,
    as = PriorParam$a.s,
    bs = PriorParam$b.s,
    atheta = PriorParam$a.theta,
    btheta = PriorParam$b.theta,
    ar = ArgsDef$AR,
    LSmu0 = Start$ls.mu0,
    LSdelta0 = Start$ls.delta0,
    LSnu0 = Start$ls.nu0,
    LStheta0 = rep(Start$ls.theta0, GPar$nBatch),
    sumByGeneAll = sum.bygene.all,
    StoreAdapt = as.numeric(ArgsDef$StoreAdapt),
    EndAdapt = ArgsDef$StopAdapt,
    PrintProgress = as.numeric(ArgsDef$PrintProgress),
    mintol_mu = ArgsDef$mintol_mu,
    mintol_delta = ArgsDef$mintol_delta,
    mintol_nu = ArgsDef$mintol_nu,
    mintol_theta = ArgsDef$mintol_theta,
    geneExponent = PriorParam$GeneExponent,
    cellExponent = PriorParam$CellExponent,
    threads = Threads
  )


  NonRegressionArgs <- list(
    prior_delta = (PriorParam$PriorDelta == "log-normal") + 1,
    s2delta = PriorParam$s2.delta,
    adelta = PriorParam$a.delta,
    bdelta = PriorParam$b.delta
  )

  RegressionArgs <- list(
    k = PriorParam$k,
    m0 = PriorParam$m,
    V0 = PriorParam$V,
    sigma2_a0 = PriorParam$a.sigma2,
    sigma2_b0 = PriorParam$b.sigma2,
    beta0 = Start$beta0,
    sigma20 = Start$sigma20,
    eta0 = PriorParam$eta,
    lambda0 = Start$lambda0,
    variance = PriorParam$variance,
    FixLocations = PriorParam$FixLocations,
    RBFMinMax = PriorParam$RBFMinMax,
    RBFLocations = PriorParam$RBFLocations
  )

  NoSpikeArgs <- list(
    sumByCellAll = sum.bycell.all,
    Constrain = ArgsDef$Constrain,
    Index = ArgsDef$Index,
    RefGene = ArgsDef$RefGene,
    RefGenes = ArgsDef$RefGenes,
    ConstrainGene = ArgsDef$ConstrainGene,
    NotConstrainGene = ArgsDef$NotConstrainGene,
    StochasticRef = as.numeric(PriorParam$StochasticRef)
  )


  # If spikes are available
  if (WithSpikes) {
    SpikeArgs <- list(
      sumByCellBio = sum.bycell.bio,
      sumByGeneBio = sum.bygene.bio,
      phi0 = Start$phi0,
      aphi = PriorParam$p.phi,
      LSphi0 = Start$ls.phi0,
      muSpikes = rowData(altExp(Data))[, 2]
    )
    MCMCArgs <- c(MCMCArgs, SpikeArgs)
    # If regression case is chosen
    if (Regression) {
      message("Running with spikes BASiCS sampler (regression case) ... \n")
      MCMCArgs <- c(MCMCArgs, RegressionArgs)
      MCMCFun <- ".BASiCS_MCMCcppReg"
    } else {
      message("Running with spikes BASiCS sampler (no regression) ... \n")
      MCMCArgs <- c(MCMCArgs, NonRegressionArgs)
      MCMCFun <- ".BASiCS_MCMCcpp"
    }
  } else {
    # If spikes are not available
    MCMCArgs <- c(MCMCArgs, NoSpikeArgs)
    if (Regression) {
      message("Running no spikes BASiCS sampler (regression case) ... \n")
      MCMCArgs <- c(MCMCArgs, RegressionArgs)
      MCMCFun <- ".BASiCS_MCMCcppRegNoSpikes"
    } else {
      message("Running no spikes BASiCS sampler (no regression) ... \n")
      MCMCArgs <- c(MCMCArgs, NonRegressionArgs)
      MCMCFun <- ".BASiCS_MCMCcppNoSpikes"
    }
  }
  Time <- system.time(Chain <- do.call(what = MCMCFun, args = MCMCArgs))
  if (Regression) {
    # Remove epsilons for genes that are not expressed in at least 2 cells
    # Discuss this with John (potentially include an optional arg about this)
    AtLeast2Cells <- Matrix::rowSums(GPar$BioCounts > 0) > 1
    Chain$epsilon[, !AtLeast2Cells] <- NA
    Chain$delta[, !AtLeast2Cells] <- NA
  }

  # Format column names of MCMC chains
  colnames(Chain$mu) <- colnames(Chain$delta) <- rownames(GPar$BioCounts)
  if (Regression) {
    colnames(Chain$epsilon) <- colnames(Chain$mu)
    Chain$lambda <- NULL # Remove to reduce storage
  }
  CellLabels <- paste0(colnames(counts(Data)), "_Batch", GPar$BatchInfo)
  colnames(Chain$s) <- CellLabels
  if (WithSpikes) {
    colnames(Chain$phi) <- CellLabels
  }
  colnames(Chain$nu) <- CellLabels
  colnames(Chain$theta) <- paste0("Batch", unique(GPar$BatchInfo))

  message(
    "-------------------------------------------------------------\n",
    "MCMC running time \n",
    "-------------------------------------------------------------\n",
    "user: ", round(Time['user.self'], 3), "\n",
    "system: ", round(Time['sys.self'], 3), "\n",
    "elapsed: ", round(Time['elapsed'], 3), "\n"
  )

  message(
    "-------------------------------------------------------------\n",
    "Output \n",
    "-------------------------------------------------------------\n"
  )

  # Convert output into a `BASiCS_Chain` object
  ChainClass <- newBASiCS_Chain(parameters = Chain[!grepl("ls.", names(Chain))])
  # Store chain and/or adaptive variances
  .BASiCS_MCMC_OutputStore(
    ChainClass = ChainClass,
    Chain = Chain,
    StoreChains = ArgsDef$StoreChains,
    StoreAdapt = ArgsDef$StoreAdapt,
    StoreDir = ArgsDef$StoreDir,
    RunName = ArgsDef$RunName
  )
  # Store reference gene information (no spikes case only)
  if (!WithSpikes) {
    if (ArgsDef$StoreChains) {
      .BASiCS_MCMC_RefFreqStore(
        Data = Data,
        Chain = Chain,
        RefGene = ArgsDef$RefGene,
        RefGenes = ArgsDef$RefGenes,
        StoreDir = ArgsDef$StoreDir,
        RunName = ArgsDef$RunName
      )
    }
  } else {
    message(
      "-------------------------------------------------------------\n",
      "BASiCS version ", packageVersion("BASiCS"), " : \n",
      "vertical integration (spikes case) \n",
      "-------------------------------------------------------------\n"
    )
  }

  # Return `BASiCS_MCMC` object
  return(ChainClass)
}
