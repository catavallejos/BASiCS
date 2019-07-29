#' @title BASiCS MCMC sampler
#'
#' @description MCMC sampler to perform Bayesian inference for single-cell
#' mRNA sequencing datasets using the model described in Vallejos et al (2015).
#'
#' @param Data A \code{\linkS4class{SingleCellExperiment}} object.
#' If \code{WithSpikes = TRUE}, this MUST be formatted to include
#' the spike-ins information (see vignette).
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
#' added spike-ins to estimate technical variability. If \code{WithSpikess = FALSE},
#' BASiCS depends on replicated experiments (batches) to estimate
#' technical variability. In this case, please supply the BatchInfo vector
#' in \code{colData(Data)}. Default: \code{WithSpikes = TRUE}.
#' @param Verbose Verbosity flag. Set FALSE to prevent printing of messages.
#' @param ... Optional parameters.
#' \describe{
#'   \item{
#'     \code{PriorDelta}
#'   }{
#'     Specifies the prior used for \code{delta}.
#'     Possible values are 'gamma' (Gamma(\code{a.theta},\code{b.theta}) prior) and
#'     'log-normal' (log-Normal(\code{0},\code{s2.delta}) prior).
#'     Default value: \code{PriorDelta = 'log-normal'}.
#'   }
#'   \item{
#'     \code{PriorParam}
#'   }{
#'     List of 7 elements, containing the hyper-parameter
#'     values required for the adopted prior (see Vallejos et al, 2015, 2016).
#'     All elements must be positive real numbers.
#'     \describe{
#'       \item{
#'         \code{s2.mu}
#'       }{
#'         Scale hyper-parameter for the
#'         log-Normal(\code{0},\code{s2.mu}) prior that is shared by all
#'         gene-specific expression rate parameters \eqn{\mu_i}.
#'         Default: \code{s2.mu = 0.5}.
#'       }
#'       \item{
#'         \code{s2.delta}
#'       }{
#'         Only used when `PriorDelta == 'log-normal'`.
#'         Scale hyper-parameter for the log-Normal(\code{0},\code{s2.delta})
#'         prior that is shared by all gene-specific over-dispersion parameters
#'         \eqn{\delta_i}. Default: \code{s2.delta = 0.5}. 
#'       }
#'       \item{
#'         \code{a.delta}
#'       }{
#'         Only used when `PriorDelta == 'gamma'`.
#'         Shape hyper-parameter for the Gamma(\code{a.delta},\code{b.delta})
#'         prior that is shared by all gene-specific biological over-dispersion
#'         parameters \eqn{\delta_i}. Default: \code{a.delta = 1}.
#'       }
#'       \item{
#'         \code{b.delta}
#'       }{
#'         Only used when `PriorDelta == 'gamma'`.
#'         Rate hyper-parameter for the Gamma(\code{a.delta},\code{b.delta})
#'         prior that is shared by all gene-specific biological over-dispersion
#'         hyper-parameters \eqn{\delta_i}. Default: \code{b.delta = 1}.
#'       }
#'       \item{
#'         \code{p.phi}
#'       }{
#'         Dirichlet hyper-parameter for the joint of all
#'         (scaled by \code{n}) cell-specific mRNA content normalising
#'         constants \eqn{\phi_j / n}.
#'         Default: \code{p.phi} \code{= rep(1, n)}.
#'       }
#'       \item{
#'         \code{a.s}
#'       }{
#'         Shape hyper-parameter for the
#'         Gamma(\code{a.s},\code{b.s}) prior that is shared by all
#'         cell-specific capture efficiency normalising constants \eqn{s_j}.
#'         Default: \code{a.s = 1}.
#'       }
#'       \item{
#'         \code{b.s}
#'       }{
#'         Rate hyper-parameter for the Gamma(\code{a.s},
#'         \code{b.s}) prior that is shared by all cell-specific capture
#'         efficiency normalising constants \eqn{s_j}.
#'         Default: \code{b.s = 1}.
#'       }
#'       \item{
#'         \code{a.theta}
#'       }{
#'         Shape hyper-parameter for the
#'         Gamma(\code{a.theta},\code{b.theta}) prior for technical noise
#'         parameter \eqn{\theta}. Default: \code{a.theta = 1}.
#'       }
#'       \item{
#'         \code{b.theta}
#'       }{
#'         Rate hyper-parameter for the
#'         Gamma(\code{a.theta},\code{b.theta}) prior for technical noise
#'         parameter \eqn{\theta}. Default: \code{b.theta = 1}.
#'       }
#'       \item{
#'         \code{eta}
#'       }{
#'         Only used when \code{Regression = TRUE}. \code{eta}
#'         specifies the degress of freedom for the residual term.
#'         Default: \code{eta = 5}.
#'       }
#'     }
#'   }
#'   \item{
#'     \code{k}
#'   }{
#'     Only used when \code{Regression = TRUE}. \code{k} specifies
#'     the number of regression Gaussian Radial Basis Functions (GRBF) used
#'     within the correlated prior adopted for gene-specific over-dispersion
#'     and mean expression paramters. Default: \code{k = 12}. 
#'   }
#'   \item{
#'     \code{FixML}
#'   }{
#'     Only used when \code{Regression = TRUE}. \code{FixML} specifies
#'     whether the locations of the \code{k - 2} Gaussian Radial Basis Functions 
#'     (GRBF) used within the correlated prior for gene-specific mean and
#'     over-dispersion are fixed a priori.
#'   }
#'   \item{
#'     \code{ml}
#'   }{
#'     Only used when \code{Regression = TRUE} and \code{FixML = TRUE}. 
#'     \code{ml} specifies
#'     the location of the \code{k - 2} Gaussian Radial Basis Functions (GRBF) 
#'     used within the correlated prior adopted for gene-specific 
#'     over-dispersion and mean expression parameters. If not supplied, 
#'     it is estimated from the starting mu values, \code{mu0}.
#'   }
#'   \item{
#'     \code{Var}
#'   }{
#'     Only used when \code{Regression = TRUE}. \code{Var}
#'     specifies the GRBF scaling parameter. Default: \code{Var = 1.2}.
#'   }
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
#'     Default: \code{StoreChains = FALSE}.
#'   }
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
#'     \code{StoreAdapt = TRUE}). Default: \code{StoreDir = getwd()}.
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
#'     progress report is suppressed.
#'   }
#'   \item{
#'     \code{Start}
#'   }{
#'     Starting values for the MCMC sampler. We do not advise
#'     to use this argument. Default options have been tuned to facilitate
#'     convergence. If changed, it must be a list containing the following
#'     elements: \code{mu0}, \code{delta0}, \code{phi0}, \code{s0},
#'     \code{nu0}, \code{theta0}, \code{ls.mu0}, \code{ls.delta0},
#'     \code{ls.phi0}, \code{ls.nu0} and \code{ls.theta0}
#'   }
#' }
#'
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
#' @export
BASiCS_MCMC <- function(Data,
                        N,
                        Thin,
                        Burn,
                        Regression,
                        WithSpikes = TRUE,
                        geneExponent = 1,
                        cellExponent = 1,
                        Verbose = TRUE,
                        ...)
{
  # Checks to ensure input arguments are valid
  HiddenBASiCS_MCMC_InputCheck(Data, N, Thin, Burn, Regression, WithSpikes)
  # Some global values used throughout the MCMC algorithm and checks
  # Numbers of cells/genes/batches + count table + design matrix for batch
  GPar <- HiddenBASiCS_MCMC_GlobalParams(Data)

  # Total counts per cell and/or gene
  sum.bycell.all <- matrixStats::rowSums2(counts(Data))
  sum.bygene.all <- matrixStats::colSums2(counts(Data))
  sum.bycell.bio <- matrixStats::rowSums2(GPar$BioCounts)
  sum.bygene.bio <- matrixStats::colSums2(GPar$BioCounts)

  # Assignment of optional parameters (default values if input not provided)
  ArgsDef <- HiddenBASiCS_MCMC_ExtraArgs(Data,
                                         Burn,
                                         GPar,
                                         Regression,
                                         WithSpikes,
                                         ...)
  # Starting values for the MCMC (parameters and adaptive variances)
  # Loaded separately to simplify calling to the elements of its list
  # Same for prior parameters
  Start <- ArgsDef$Start
  PriorParam <- ArgsDef$PriorParam
  if (WithSpikes) {
    SpikeInput <- metadata(Data)$SpikeInput
    if (Regression) {
      if (Verbose) {
        message("Running with spikes BASiCS sampler (regression case) ... \n")
      }
      Time <- system.time(
        Chain <- HiddenBASiCS_MCMCcppReg(
          N = N,
          Thin = Thin,
          Burn = Burn,
          Counts = GPar$BioCounts,
          BatchDesign = GPar$BatchDesign,
          muSpikes = SpikeInput,
          mu0 = Start$mu0,
          delta0 = Start$delta0,
          phi0 = Start$phi0,
          s0 = Start$s0,
          nu0 = Start$nu0,
          theta0 = rep(Start$theta0, GPar$nBatch),
          s2mu = PriorParam$s2.mu,
          aphi = PriorParam$p.phi,
          as = PriorParam$a.s,
          bs = PriorParam$b.s,
          atheta = PriorParam$a.theta,
          btheta = PriorParam$b.theta,
          ar = ArgsDef$AR,
          LSmu0 = Start$ls.mu0,
          LSdelta0 = Start$ls.delta0,
          LSphi0 = Start$ls.phi0,
          LSnu0 = Start$ls.nu0,
          LStheta0 = rep(Start$ls.theta0, GPar$nBatch),
          sumByCellAll = sum.bycell.all,
          sumByCellBio = sum.bycell.bio,
          sumByGeneBio = sum.bygene.bio,
          sumByGeneAll = sum.bygene.all,
          StoreAdapt = as.numeric(ArgsDef$StoreAdapt),
          EndAdapt = ArgsDef$StopAdapt,
          PrintProgress = as.numeric(ArgsDef$PrintProgress),
          k = ArgsDef$k,
          m0 = PriorParam$m,
          V0 = PriorParam$V,
          sigma2_a0 = PriorParam$a.sigma2,
          sigma2_b0 = PriorParam$b.sigma2,
          beta0 = Start$beta0,
          sigma20 = Start$sigma20,
          eta0 = PriorParam$eta,
          lambda0 = Start$lambda0,
          variance = ArgsDef$variance,
          ml = ArgsDef$ml,
          FixML = ArgsDef$FixML,
          geneExponent = geneExponent,
          cellExponent = cellExponent,
          mintol_mu = ArgsDef$mintol_mu, 
          mintol_delta = ArgsDef$mintol_delta,
          mintol_nu = ArgsDef$mintol_nu,
          mintol_theta = ArgsDef$mintol_theta
        )
      )
      # Remove epsilons for genes that are not expressed in at least 2 cells
      # Discuss this with John (potentially include an optional arg about this)
      AtLeast2Cells <- matrixStats::rowSums2(ifelse(GPar$BioCounts > 0, 1, 0)) > 1
      Chain$epsilon[,!AtLeast2Cells] <- NA
    } else {
      if (Verbose) {
        message("Running with spikes BASiCS sampler (no regression) ... \n")
      }
      Time <- system.time(
        Chain <- HiddenBASiCS_MCMCcpp(
          N = N,
          Thin = Thin,
          Burn = Burn,
          Counts = GPar$BioCounts,
          BatchDesign = GPar$BatchDesign,
          muSpikes = SpikeInput,
          mu0 = Start$mu0,
          delta0 = Start$delta0,
          phi0 = Start$phi0,
          s0 = Start$s0,
          nu0 = Start$nu0,
          theta0 = rep(Start$theta0, GPar$nBatch),
          s2mu = PriorParam$s2.mu,
          adelta = PriorParam$a.delta,
          bdelta = PriorParam$b.delta,
          s2delta = PriorParam$s2.delta,
          prior_delta = ArgsDef$PriorDeltaNum,
          aphi = PriorParam$p.phi,
          as = PriorParam$a.s,
          bs = PriorParam$b.s,
          atheta = PriorParam$a.theta,
          btheta = PriorParam$b.theta,
          ar = ArgsDef$AR,
          LSmu0 = Start$ls.mu0,
          LSdelta0 = Start$ls.delta0,
          LSphi0 = Start$ls.phi0,
          LSnu0 = Start$ls.nu0,
          LStheta0 = rep(Start$ls.theta0, GPar$nBatch),
          sumByCellAll = sum.bycell.all,
          sumByCellBio = sum.bycell.bio ,
          sumByGeneAll = sum.bygene.all,
          sumByGeneBio = sum.bygene.bio,
          StoreAdapt = as.numeric(ArgsDef$StoreAdapt),
          EndAdapt = ArgsDef$StopAdapt,
          PrintProgress = as.numeric(ArgsDef$PrintProgress),
          geneExponent = geneExponent,
          cellExponent = cellExponent,
          mintol_mu = ArgsDef$mintol_mu, 
          mintol_delta = ArgsDef$mintol_delta,
          mintol_nu = ArgsDef$mintol_nu,
          mintol_theta = ArgsDef$mintol_theta
        )
      )
    }
  } else {
    # Definition of parameters that are specific to the no-spikes case
    if (Regression) {
      if (Verbose) {
        message("Running no spikes BASiCS sampler (regression case) ... \n")
      }
      Time <- system.time(
        Chain <- HiddenBASiCS_MCMCcppRegNoSpikes(
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
          sumByCellAll = sum.bycell.bio,
          sumByGeneAll = sum.bygene.bio,
          StoreAdapt = as.numeric(ArgsDef$StoreAdapt),
          EndAdapt = ArgsDef$StopAdapt,
          PrintProgress = as.numeric(ArgsDef$PrintProgress),
          k = ArgsDef$k,
          m0 = PriorParam$m,
          V0 = PriorParam$V,
          sigma2_a0 = PriorParam$a.sigma2,
          sigma2_b0 = PriorParam$b.sigma2,
          beta0 = Start$beta0,
          sigma20 = Start$sigma20,
          eta0 = PriorParam$eta,
          lambda0 = Start$lambda0,
          variance = ArgsDef$variance,
          Constrain = ArgsDef$Constrain,
          Index = ArgsDef$Index,
          RefGene = ArgsDef$RefGene,
          RefGenes = ArgsDef$RefGenes,
          ConstrainGene = ArgsDef$ConstrainGene,
          NotConstrainGene = ArgsDef$NotConstrainGene,
          ConstrainType = ArgsDef$ConstrainType,
          StochasticRef = as.numeric(ArgsDef$StochasticRef),
          ml = ArgsDef$ml,
          FixML = ArgsDef$FixML,
          geneExponent = geneExponent,
          cellExponent = cellExponent,
          mintol_mu = ArgsDef$mintol_mu, 
          mintol_delta = ArgsDef$mintol_delta,
          mintol_nu = ArgsDef$mintol_nu,
          mintol_theta = ArgsDef$mintol_theta
        )
      )
      # Remove epsilons for genes that are not expressed in at least 2 cells
      # Discuss this with John (potentially include an optional arg about this)
      AtLeast2Cells <- matrixStats::rowSums2(ifelse(GPar$BioCounts > 0, 1, 0)) > 1
      Chain$epsilon[,!AtLeast2Cells] <- NA
    }
    else {
      if (Verbose) {
        message("Running no spikes BASiCS sampler (no regression) ... \n")
      }
      Time <- system.time(
        Chain <- HiddenBASiCS_MCMCcppNoSpikes(
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
          s2mu = PriorParam$s2.mu,
          adelta = PriorParam$a.delta,
          bdelta = PriorParam$b.delta,
          s2delta = PriorParam$s2.delta,
          prior_delta = ArgsDef$PriorDeltaNum,
          as = PriorParam$a.s,
          bs = PriorParam$b.s,
          atheta = PriorParam$a.theta,
          btheta = PriorParam$b.theta,
          ar = ArgsDef$AR,
          LSmu0 = Start$ls.mu0,
          LSdelta0 = Start$ls.delta0,
          LSnu0 = Start$ls.nu0,
          LStheta0 = rep(Start$ls.theta0, GPar$nBatch),
          sumByCellAll = sum.bycell.bio,
          sumByGeneAll = sum.bygene.bio,
          StoreAdapt = as.numeric(ArgsDef$StoreAdapt),
          EndAdapt = ArgsDef$StopAdapt,
          PrintProgress = as.numeric(ArgsDef$PrintProgress),
          Constrain = ArgsDef$Constrain,
          Index = ArgsDef$Index,
          RefGene = ArgsDef$RefGene,
          RefGenes = ArgsDef$RefGenes,
          ConstrainGene = ArgsDef$ConstrainGene,
          NotConstrainGene = ArgsDef$NotConstrainGene,
          ConstrainType = ArgsDef$ConstrainType,
          StochasticRef = as.numeric(ArgsDef$StochasticRef),
          geneExponent = geneExponent,
          cellExponent = cellExponent,
          mintol_mu = ArgsDef$mintol_mu,
          mintol_delta = ArgsDef$mintol_delta,
          mintol_nu = ArgsDef$mintol_nu,
          mintol_theta = ArgsDef$mintol_theta
        )
      )
    }
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
  if (Verbose) {
    message("-------------------------------------------------------------\n",
            "MCMC running time \n",
            "-------------------------------------------------------------\n",
            "user: ", round(Time['user.self'], 3), "\n",
            "system: ", round(Time['sys.self'], 3), "\n",
            "elapsed: ", round(Time['elapsed'], 3), "\n")
    message("-------------------------------------------------------------\n",
            "Output \n",
            "-------------------------------------------------------------\n")
  }
  # Convert output into a `BASiCS_Chain` object
  indParam <- !(grepl("ls.", names(Chain)) | names(Chain) == "ml")

  ChainClass <- newBASiCS_Chain(parameters = Chain[indParam])

  # Store chain and/or adaptive variances
  HiddenBASiCS_MCMC_OutputStore(ChainClass,
                                Chain,
                                ArgsDef$StoreChains,
                                ArgsDef$StoreAdapt,
                                ArgsDef$StoreDir,
                                ArgsDef$RunName)

  # Store reference gene information (no spikes case only)
  if (!WithSpikes) {
    if (ArgsDef$StoreChains)
      HiddenBASiCS_MCMC_RefFreqStore(Data,
                                     Chain,
                                     ArgsDef$RefGene,
                                     ArgsDef$RefGenes,
                                     ArgsDef$ConstrainType,
                                     ArgsDef$StoreDir,
                                     ArgsDef$RunName)

    if (Verbose) {
      message("-------------------------------------------------------------\n",
              "BASiCS version ", packageVersion("BASiCS"), " : \n",
              "vertical integration (spikes case) \n",
              "-------------------------------------------------------------\n")
    }
  }
  # Return `BASiCS_MCMC` object
  return(ChainClass)
}
