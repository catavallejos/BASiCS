#' Run divide and conquer MCMC with BASiCS
#' 
#' Performs MCMC inference on batches of data. \code{Data} is divided
#' into \code{NSubsets} batches, and \code{BASiCS_MCMC} is run on each
#' batch separately.
#' 
#' Subsets are chosen such that the average library size (when partitioning
#' by cells) or average count (when partitioning by genes) is not significantly
#' different between batches, at a significance level \code{Alpha}.
#' 
#' @param Data SingleCellExperiemnt object
#' @param NSubsets The number of batches to create and perform MCMC inference
#' with.
#' @param SubsetBy A character value specifying whether batches should consist
#' of a subset of the cells in \code{Data} (when \code{SubsetBy="cell"})
#' or a subset of the genes in \code{Data} (when \code{SubsetBy="gene"}).
#' @param Alpha A numeric value specifying the statistical significance level
#' used to determine whether the average library size or average count
#' are significantly different between batches.
#' @param WithSpikes,Regression,PriorParam See \code{\link{BASiCS_MCMC}}.
#' @param BPPARAM A \code{\link{BiocParallelParam}} instance.
#' @param ... Passed to  \code{\link{BASiCS_MCMC}}. All arguments required by
#' \code{\link{BASiCS_MCMC}} must be supplied here, for example
#' \code{N}, \code{Thin}, \code{Burn}.
#' @return A list of \linkS4class{BASiCS_Chain} objects.
#' @examples
#'  bp <- BiocParallel::SnowParam()
#'  Data <- BASiCS_MockSCE()
#'  BASiCS_DivideAndConquer(
#'    Data, 
#'    NSubsets = 2,
#'    SubsetBy = "gene",
#'    N = 8,
#'    Thin = 2,
#'    Burn = 4,
#'    WithSpikes = TRUE,
#'    Regression = TRUE,
#'    BPPARAM = bp
#'  )
#' @references
#' Simple, Scalable and Accurate Posterior Interval Estimation
#' Cheng Li and Sanvesh Srivastava and David B. Dunson
#' arXiv (2016)
#' @export
BASiCS_DivideAndConquer <- function(
    Data,
    NSubsets = 5,
    SubsetBy = c("cell", "gene"),
    Alpha = 0.05,
    WithSpikes,
    Regression,
    BPPARAM = BiocParallel::bpparam(),
    PriorParam = BASiCS_PriorParam(
      Data,
      PriorMu = "EmpiricalBayes"
    ),
    ...
    ) {

  ## These values should be passed in.
  L <- 12
  eta <- 5
  astwo <- 2
  bstwo <- 2
  start <- .BASiCS_MCMC_Start(
    Data,
    PriorParam = PriorParam,
    Regression = Regression,
    WithSpikes = WithSpikes
  )

  if (Regression) {  
    Locations <- .estimateRBFLocations(
      start$mu0,
      L,
      RBFMinMax = FALSE
    )
    PriorParam$RBFLocations <- Locations
  }

  ## To do:
  ## Automatically determine NSubsets?
  ## Balance subsets for variables/batches?
  SubsetBy <- match.arg(SubsetBy)
  cat("Generating partitions...\n")
  Subsets <- .generateSubsets(
    Data = Data,
    NSubsets = NSubsets,
    SubsetBy = SubsetBy,
    Alpha = Alpha,
    WithSpikes = WithSpikes
  )
  cat("Starting MCMC...\n")
  output <- bplapply(
    Subsets,
    function(Subset) {
      ind <- match(rownames(Subset), rownames(Data))
      PriorParam$mu.mu <- PriorParam$mu.mu[ind]
      ind <- match(colnames(Subset), colnames(Data))
      PriorParam$p.phi <- PriorParam$p.phi[ind]
      BASiCS_MCMC(
        Data = Subset,
        WithSpikes = WithSpikes,
        Regression = Regression,
        PriorParam = PriorParam,
        Threads = 1,
        ...
      )
    },
    BPPARAM = BPPARAM
  )
  if (any(vapply(output, class, character(1)) == "try-error")) {
    warning("Some MCMC runs failed!")
  }
  output
}
