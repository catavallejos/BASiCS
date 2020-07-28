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
#' (must be a positive value, between 0 and 1). Defaults: 0.9 for HVG (top 10 
#' percent), 0.1 for LVG (bottom 10 percent)
#' @param VarThreshold Variance contribution threshold
#' (must be a positive value, between 0 and 1). This is only used when the 
#' BASiCS non-regression model was used to generate the Chain object.
#' @param ProbThreshold Optional parameter. Posterior probability threshold
#' (must be a positive value, between 0 and 1). If \code{EFDR = NULL}, the 
#' posterior probability threshold for the test will be set to \code{ProbThreshold}
#' @param EFDR Target for expected false discovery rate related
#' to HVG/LVG detection. If \code{EFDR = NULL}, EFDR calibration is
#' not performed and the posterior probability threshold is set equal to
#' \code{ProbThreshold}.Default \code{EFDR = 0.10}.
#' @param OrderVariable Ordering variable for output.
#' Possible values: \code{'GeneIndex'}, \code{'GeneName'} and \code{'Prob'}.
#' @param Plot If \code{Plot = TRUE} error control and
#' expression versus HVG/LVG probability plots are generated
#' @param ... Graphical parameters (see \code{\link[graphics]{par}}).
#'
#' @return \code{BASiCS_DetectHVG} returns a list of 4 elements:
#' \describe{
#'     \item{
#'         \code{Table}
#'     }{
#'         Matrix whose columns can contain
#'     }
#'     \describe{
#'         \item{
#'             \code{GeneIndex}
#'         }{
#'             Vector of length \code{q.bio}.
#'             Gene index as in the order present in the analysed
#'             \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#'         }
#'         \item{
#'             \code{GeneName}
#'         }{
#'             Vector of length \code{q.bio}.
#'             Gene name as in the order present in the analysed
#'              \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#'         }
#'         \item{
#'             \code{Mu}
#'         }{
#'             Vector of length \code{q.bio}. For each biological gene,
#'              posterior median of gene-specific mean expression
#'              parameters \eqn{\mu_i}.
#'          }
#'         \item{
#'             \code{Delta}
#'         }{
#'             Vector of length \code{q.bio}. For each biological
#'             gene, posterior median of gene-specific biological
#'             over-dispersion parameter \eqn{\delta_i}.
#'          }
#'         \item{
#'             \code{Sigma}
#'         }{
#'             Vector of length \code{q.bio}.
#'             For each biological gene, proportion of the total variability
#'             that is due to a biological heterogeneity component.
#'         }
#'         \item{
#'             \code{Epsilon}
#'         }{
#'             Vector of length \code{q.bio}.
#'             For each biological gene, posterior median of gene-specific
#'             residual over-dispersion parameter \eqn{\epsilon_i}.
#'         }
#'         \item{
#'             \code{Prob}
#'         }{
#'             Vector of length \code{q.bio}.
#'             For each biological gene, probability of being highly variable
#'             according to the given thresholds.
#'         }
#'         \item{
#'             \code{HVG}
#'         }{
#'             Vector of length \code{q.bio}.
#'             For each biological gene, indicator of being detected as highly
#'             variable according to the given thresholds.
#'         }
#'         \item{
#'             \code{LVG}
#'         }{
#'             Vector of length \code{q.bio}.
#'             For each biological gene, indicator of being detected as lowly
#'             variable according to the given thresholds.
#'         }
#'     }
#'     \item{
#'         \code{ProbThreshold}
#'     }{
#'         Posterior probability threshold.
#'     }
#'     \item{
#'         \code{EFDR}
#'     }{
#'         Expected false discovery rate for the given thresholds.
#'     }
#'     \item{
#'         \code{EFNR}
#'     }{
#'         Expected false negative rate for the given thresholds.
#'     }
#' }
#'
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
    PercentileThreshold = 0.1,
    VarThreshold = 0.5,
    ProbThreshold = 0.5,
    EFDR = 0.1,
    OrderVariable = c("Prob", "GeneIndex","GeneName"),
    Plot = FALSE,
    ...
  ) {
  
  # Check valid input values
  .HeaderDetectHVG_LVG(Chain,
                       PercentileThreshold,
                       VarThreshold,
                       ProbThreshold,
                       EFDR,
                       Plot)
  
  # Define LVG/HVG criteria
  Task <- match.arg(Task)
  OrderVariable <- match.arg(OrderVariable)
  operator <- if (Task == "HVG") `>` else `<`
  if (!is.null(Chain@parameters$beta) & !is.null(PercentileThreshold)) {
    Method <- "Percentile"
  } else {
    Method <- "Variance"
  }
  
  # Prepare template for output table
  GeneIndex <- seq_along(Chain@parameters$mu)
  GeneName <- colnames(Chain@parameters$mu)
  Table <- cbind.data.frame(
    GeneIndex = GeneIndex,
    GeneName = GeneName,
    Mu = matrixStats::colMedians(Chain@parameters$mu),
    Delta = matrixStats::colMedians(Chain@parameters$mu)
  )

  if (Method == "Percentile") {
    
    # Find the epsilon threshold that correspond to the 'PercentileThreshold'
    Epsilon <- matrixStats::colMedians(Chain@parameters$epsilon)
    EpsilonThreshold <- stats::quantile(
      Epsilon,
      PercentileThreshold,
      na.rm = TRUE
    )
    Table <- cbind.data.frame(Table, Epsilon = Epsilon)
    # HVG probability for a given epsilon threshold
    Variable <- operator(Chain@parameters$epsilon, EpsilonThreshold)
    
  } else {
    
    # Variance decomposition
    VarDecomp <- HiddenVarDecomp(Chain)
    Table <- cbind.data.frame(Table, 
      Sigma = matrixStats::colMedians(VarDecomp$BioVarGlobal))
    # H/LVG probability for a given variance threshold
    Variable <- operator(VarDecomp$BioVarGlobal, VarThreshold)
    
  }
  # Calculate tail posterior probabilities
  Prob <- matrixStats::colMeans2(Variable)
  
  # EFDR calibration
  Aux <- .ThresholdSearch(Prob[!is.na(Prob)], 
                          ProbThreshold, 
                          EFDR, 
                          Task = paste(Task, "detection"))

  # Output preparation
  VG <- Prob > Aux$OptThreshold[1]
  Table <- cbind.data.frame(Table,
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
        EFNR = Aux$OptThreshold[[2]],
        Method = Method,
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
BASiCS_DetectLVG <- function(..., PercentileThreshold =  0.1) {
  # .Deprecated("BASiCS_DetectVG")
  BASiCS_DetectVG(
    ...,
    PercentileThreshold = PercentileThreshold,
    Task = "LVG"
  )
}

#' @name BASiCS_DetectHVGLVG
#' @aliases BASiCS_DetectHVG BASiCS_DetectLVG BASiCS_DetectHVG_LVG
#' @rdname BASiCS_DetectVG
#' @export
BASiCS_DetectHVG <- function(..., PercentileThreshold =  0.9) {
  # .Deprecated("BASiCS_DetectVG")
  BASiCS_DetectVG(
    ...,
    PercentileThreshold = PercentileThreshold,
    Task = "HVG"
  )
}

#' Plots of HVG/LVG search.
#' @param object \linkS4class{BASiCS_ResultVG} object.
#' @param Plot Character scalar specifying the type of plot to be made.
#' Options are "Grid" and "VG".
#' @param ... Optional graphical parameters passed to \code{.VGPlot} 
#' (internal function). 
#' 
#' @return A plot.
#' @examples
#' data(ChainSC)
#'
#' # Highly and lowly variable genes detection (within a single group of cells)
#' DetectHVG <- BASiCS_DetectHVG(ChainSC, VarThreshold = 0.60,
#'                               EFDR = 0.10, Plot = TRUE)
#' BASiCS_PlotVG(DetectHVG)
#' @export
BASiCS_PlotVG <- function(object, Plot = c("Grid", "VG"), ...) {
  Plot <- match.arg(Plot)
  if (Plot == "Grid") {
    .VGGridPlot(
      ProbThresholds = object@ProbThresholds,
      EFDRgrid = object@EFDRgrid,
      EFNRgrid = object@EFNRgrid,
      EFDR = object@EFDR
    )
  } else if (Plot == "VG") {
    .VGPlot(
        Task = object@Name,
        Mu = object@Table$Mu,
        Prob = object@Table$Prob,
        OptThreshold = object@ProbThreshold,
        Hits = object@Table[[object@Name]],
        ...
    )
  }
}






