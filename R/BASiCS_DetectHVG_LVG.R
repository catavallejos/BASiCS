#' @name BASiCS_DetectVG
#' @aliases BASiCS_DetectHVG BASiCS_DetectLVG BASiCS_DetectHVG_LVG
#'
#' @title Detection method for highly (HVG) and lowly (LVG) variable genes
#'
#' @description Functions to detect highly and lowly variable genes. If the 
#' BASiCS_Chain object was generated using the regression approach,
#' BASiCS finds the top highly variable genes based on the posteriors of the 
#' epsilon parameters. Otherwise, the old approach is used, which initially 
#' performs a variance decomposition.
#'
#' @param Chain an object of class \code{\linkS4class{BASiCS_Chain}}
#' @param Task Search for highly variable genes (\code{Task="HVG"})
#' or lowly variable genes (\code{Task="LVG"}).
#' @param PercentileThreshold Threshold to detect a percentile of variable genes
#' (must be a positive value, between 0 and 1). 
#' Default: \code{PercentileThreshold = NULL}.
#' @param VarThreshold Variance contribution threshold
#' (must be a positive value, between 0 and 1). This is only used when the 
#' BASiCS non-regression model was used to generate the Chain object.
#' Default: \code{VarThreshold = NULL}.
#' @param ProbThreshold Optional parameter. Posterior probability threshold
#' (must be a positive value, between 0 and 1). If \code{EFDR = NULL}, the 
#' posterior probability threshold for the test will be set to
#' \code{ProbThreshold}.
#' @param EpsilonThreshold Threshold for residual overdispersion above which
#' 
#' @param EFDR Target for expected false discovery rate related
#' to HVG/LVG detection. If \code{EFDR = NULL}, EFDR calibration is
#' not performed and the posterior probability threshold is set equal to
#' \code{ProbThreshold}. Default \code{EFDR = 0.10}.
#' @param OrderVariable Ordering variable for output.
#' Possible values: \code{'GeneIndex'}, \code{'GeneName'} and \code{'Prob'}.
#' Default \code{ProbThreshold = 'Prob'}
#' @param Plot If \code{Plot = TRUE} error control and
#' expression versus HVG/LVG probability plots are generated.
#' @param MinESS The minimum effective sample size for a gene to be included 
#' in the HVG or LVG tests. This helps to remove genes with poor mixing from
#' detection of HVGs/LVGs.
#' Default is 100. If set to NA, genes are
#' not checked for effective sample size the tests are performed.
#' @param ... Graphical parameters (see \code{\link[graphics]{par}}).
#'
#' @return An object of class \code{\link[BASiCS]{BASiCS_ResultVG}}.
#' 
#' @details
#' In some cases, the EFDR calibration step may fail to find probability
#' threshold that controls the EFDR at the chosen level. In cases like 
#' @examples
#' 
#' # Loads short example chain (non-regression implementation)
#' data(ChainSC)
#'
#' # Highly and lowly variable genes detection (within a single group of cells)
#' DetectHVG <- BASiCS_DetectHVG(ChainSC, VarThreshold = 0.60,
#'                               EFDR = 0.10, Plot = TRUE)
#' DetectLVG <- BASiCS_DetectLVG(ChainSC, VarThreshold = 0.40,
#'                               EFDR = 0.10, Plot = TRUE)
#'                               
#' # Loads short example chain (regression implementation)
#' data(ChainSCReg)
#'
#' # Highly and lowly variable genes detection (within a single group of cells)
#' DetectHVG <- BASiCS_DetectHVG(ChainSCReg, PercentileThreshold = 0.90,
#'                               EFDR = 0.10, Plot = TRUE)
#' DetectLVG <- BASiCS_DetectLVG(ChainSCReg, PercentileThreshold = 0.10,
#'                               EFDR = 0.10, Plot = TRUE)
#' 
#' ## Highly and lowly variable genes detection based on residual overdispersion
#' ## threshold
#' DetectHVG <- BASiCS_DetectHVG(ChainSCReg, EpsilonThreshold = log(2), Plot = TRUE)
#' DetectLVG <- BASiCS_DetectLVG(ChainSCReg, EpsilonThreshold = -log(2), Plot = TRUE)
#'
#' @details See vignette
#'
#'
#' @seealso \code{\linkS4class{BASiCS_Chain}}
#'
#' @author Catalina A. Vallejos \email{cnvallej@@uc.cl}
#' @author Nils Eling \email{eling@@ebi.ac.uk}
#'
#' @references
#'
#' Vallejos, Marioni and Richardson (2015). PLoS Computational Biology.
#'
#' @rdname BASiCS_DetectVG
#' @export
BASiCS_DetectVG <- function(
    Chain,
    Task = c("HVG", "LVG"),
    PercentileThreshold = NULL, # 0.1,
    VarThreshold = NULL, # 0.5,
    ProbThreshold = 2/3, # 0.5,
    EpsilonThreshold = NULL,
    EFDR = 0.1,
    OrderVariable = c("Prob", "GeneIndex", "GeneName"),
    Plot = FALSE,
    MinESS = 100,
    ...
  ) {
  
  # Check valid input values
  Task <- match.arg(Task)
  .HeaderDetectHVG_LVG(
    Chain = Chain,
    PercentileThreshold = PercentileThreshold,
    EpsilonThreshold = EpsilonThreshold,
    VarThreshold = VarThreshold,
    Plot = Plot
  )
  OrderVariable <- match.arg(OrderVariable)
  .CheckProbEFDR(ProbThreshold, EFDR)
  
  # Define LVG/HVG criteria
  operator <- if (Task == "HVG") `>` else `<`
  if (!is.null(Chain@parameters$beta) & !is.null(PercentileThreshold)) {
    Method <- "Percentile"
  } else if (!is.null(EpsilonThreshold)) {
    Method <- "Epsilon"
  } else {
    Method <- "Variance"
  }

  # Prepare template for output table
  GeneIndex <- seq_along(Chain@parameters$mu)
  GeneName <- colnames(Chain@parameters$mu)
  Table <- cbind.data.frame(
    GeneIndex = GeneIndex,
    GeneName = GeneName,
    Mu = as.vector(matrixStats::colMedians(Chain@parameters$mu)),
    Delta = as.vector(matrixStats::colMedians(Chain@parameters$delta))
  )
  

  if (Method == "Percentile") {
    # Find the epsilon threshold that correspond to the 'PercentileThreshold'
    Epsilon <- as.vector(matrixStats::colMedians(Chain@parameters$epsilon))
    EpsilonThreshold <- stats::quantile(
      Epsilon,
      PercentileThreshold,
      na.rm = TRUE
    )
    Table <- cbind.data.frame(Table, Epsilon = Epsilon)
    # Auxiliary variable to calculate H/LVG prob for a given epsilon threshold
    ProbAux <- operator(Chain@parameters$epsilon, EpsilonThreshold)
    Threshold <- PercentileThreshold 
    
  } else if (Method == "Epsilon") {
    ## Directly supply epsilon threshold
    Epsilon <- as.vector(matrixStats::colMedians(Chain@parameters$epsilon))
    Table <- cbind.data.frame(Table, Epsilon = Epsilon)
    # Auxiliary variable to calculate H/LVG prob for the given epsilon threshold
    ProbAux <- operator(Chain@parameters$epsilon, EpsilonThreshold)
    Threshold <- EpsilonThreshold
    
  } else if (Method == "Variance") {
    # Variance decomposition
    VarDecomp <- HiddenVarDecomp(Chain)
    Table <- cbind.data.frame(Table, 
      Sigma = as.vector(matrixStats::colMedians(VarDecomp$BioVarGlobal)))
    # Auxiliary variable to calculate H/LVG prob for a given variance threshold
    ProbAux <- operator(VarDecomp$BioVarGlobal, VarThreshold)
    Threshold <- VarThreshold
  }

  # Calculate tail posterior probabilities
  Prob <- as.vector(matrixStats::colMeans2(ProbAux))

  # EFDR calibration
  Aux <- .ThresholdSearch(
    Prob[!is.na(Prob)],
    ProbThreshold,
    EFDR,
    Task = paste(Task, "detection")
  )
   
  # Output preparation
  VG <- Prob > Aux$OptThreshold[1]
  Table <- cbind.data.frame(
    Table,
    Prob = Prob,
    VG = VG,
    stringsAsFactors = FALSE
  )
  colnames(Table)[[which(colnames(Table)=="VG")]] <- Task

  # Re-order the table of results
  orderVar <- switch(
    OrderVariable,
    "GeneIndex" = GeneIndex,
    "GeneName" = GeneName,
    "Prob" = Prob
  )
  Table <- Table[order(orderVar, decreasing = TRUE, na.last = TRUE), ]

  # output object
  out <- new("BASiCS_ResultVG",
    Table = Table,
    Name = Task,
    ProbThreshold = Aux$OptThreshold[[1]],
    ProbThresholds = Aux$ProbThresholds,
    EFDRgrid = Aux$EFDRgrid,
    EFNRgrid = Aux$EFNRgrid,
    EFDR = Aux$OptThreshold[[2]],
    EFNR = Aux$OptThreshold[[3]],
    Method = Method,
    Threshold = Threshold, 
    RowData = DataFrame(GeneName = GeneName)
  )
  
  if (Plot) {
    Plots <- list()
    # EFDR / EFNR plot
    Plots$Grid <- BASiCS_PlotVG(out, Plot = "Grid")
    # Output plot : mean vs prob
    Plots$VG <- BASiCS_PlotVG(out, Plot = "VG")
    # Append to existing output object
    out@Extras = list(Plots = Plots)
  }
  return(out)
}

#' @name BASiCS_DetectHVGLVG
#' @aliases BASiCS_DetectHVG BASiCS_DetectLVG BASiCS_DetectHVG_LVG
#' @rdname BASiCS_DetectVG
#' @export
BASiCS_DetectLVG <- function(Chain, ...) {
  # .Deprecated("BASiCS_DetectVG")
  BASiCS_DetectVG(Chain = Chain, Task = "LVG", ...)
}

#' @name BASiCS_DetectHVGLVG
#' @aliases BASiCS_DetectHVG BASiCS_DetectLVG BASiCS_DetectHVG_LVG
#' @rdname BASiCS_DetectVG
#' @export
BASiCS_DetectHVG <- function(Chain, ...) {
  # .Deprecated("BASiCS_DetectVG")
  BASiCS_DetectVG(Chain = Chain, Task = "HVG", ...)
}
