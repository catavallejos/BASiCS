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
#' (must be a positive value, between 0 and 1)
#' @param EFDR Target for expected false discovery rate related
#' to HVG/LVG detection (default = 0.10)
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
  
  Task <- match.arg(Task)
  .HeaderDetectHVG_LVG(Chain,
                       PercentileThreshold,
                       VarThreshold,
                       ProbThreshold,
                       EFDR,
                       Plot)
  

  OrderVariable <- match.arg(OrderVariable)
  
  if (!is.null(Chain@parameters$beta) & !is.null(PercentileThreshold)) {
    Method <- "Percentile"
  } else {
    Method <- "Variance"
  }
  GeneIndex <- seq_along(Chain@parameters$mu)
  GeneName <- colnames(Chain@parameters$mu)
  Table <- cbind.data.frame(
    GeneIndex = GeneIndex,
    GeneName = GeneName
  )

  if (Method == "Percentile") {
    # Find the epsilon threshold that correspond to the 'PercentileThreshold'
    
    nGenes <- ncol(Chain@parameters$epsilon)
    
    Epsilon <- matrixStats::colMedians(Chain@parameters$epsilon)
    EpsilonThreshold <- stats::quantile(
      Epsilon,
      PercentileThreshold,
      na.rm = TRUE
    )
    
    # HVG probability for a given epsilon threshold
    operator <- if (Task == "HVG") `>` else `<`
    Variable <- operator(Chain@parameters$epsilon, EpsilonThreshold)
    Threshold <- ProbThreshold
    Table <- cbind.data.frame(Table, Epsilon = Epsilon)
  } else {
    # Variance decomposition
    VarDecomp <- HiddenVarDecomp(Chain)

    ## outputs
    Sigma <- matrixStats::colMedians(VarDecomp$BioVarGlobal)
    Mu <- matrixStats::colMedians(Chain@parameters$mu)
    Delta <- matrixStats::colMedians(Chain@parameters$delta)

    # H/LVG probability for a given variance threshold
    operator <- if (Task == "HVG") `>` else `<`
    Variable <- operator(VarDecomp$BioVarGlobal, VarThreshold)
    Threshold <- VarThreshold
    Table <- cbind.data.frame(Table,
      Delta = Delta,
      Sigma = Sigma
    )
  }
  Prob <- matrixStats::colMeans2(Variable)
  
  # Threshold search
  Aux <- .ThresholdSearch(Prob[!is.na(Prob)], 
                          ProbThreshold, 
                          EFDR, 
                          Task = paste(Task, "detection"), 
                          Suffix = "")
  
  EFDRgrid <- Aux$EFDRgrid
  EFNRgrid <- Aux$EFNRgrid
  OptThreshold <- Aux$OptThreshold
  
  # Output preparation
  Mu <- matrixStats::colMedians(Chain@parameters$mu)
  VG <- Prob > OptThreshold[1]
  
  
  Table <- cbind.data.frame(Table,
    Mu = Mu,
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
  
  Plots <- list()
  ProbThresholds <- seq(0.5, 0.9995, by = 0.00025)
  if (Plot) {
    # EFDR / EFNR plot
    Plots[[1]] <- .VGGridPlot(
      ProbThresholds = ProbThresholds,
      EFDRgrid = EFDRgrid,
      EFNRgrid = EFNRgrid,
      EFDR = EFDR
    )
    
    # Output plot : mean vs prob
    Plots <- c(Plots, 
      list(
        .VGPlot(
          Task = Task,
          Mu = Mu,
          Prob = Prob,
          OptThreshold = OptThreshold,
          Hits = VG,
          ...
        )
      )
    )
  }

  new("BASiCS_ResultVG",
    Table = Table,
    Name = Task,
    ProbThreshold = OptThreshold[[1]],
    ProbThresholds = ProbThresholds,
    Threshold = Threshold,
    EFDRgrid = EFDRgrid,
    EFNRgrid = EFNRgrid,
    EFDR = OptThreshold[[2]],
    EFNR = OptThreshold[[2]],
    Method = Method,
    RowData = DataFrame(GeneName = GeneName),
    Extras = list(Plots = Plots)
  )
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
#' @param ... Passed to \code{.VGPlot}.
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

.VGGridPlot <- function(ProbThresholds, EFDRgrid, EFNRgrid, EFDR) {
  ggplot2::ggplot() +
    ggplot2::geom_line(
      ggplot2::aes(ProbThresholds, EFDRgrid, color = "EFDR")
    ) +
    ggplot2::geom_line(
      ggplot2::aes(ProbThresholds, EFNRgrid, color = "EFNR")
    ) +
    ggplot2::geom_hline(
      ggplot2::aes(color = "Target EFDR", yintercept = EFDR)
    ) +
    ggplot2::scale_color_brewer(palette = "Set1", name = NULL) +
    ggplot2::labs(x = "Probability threshold", y = "Error rate") +
    ggplot2::ylim(c(0, 1))
}


.VGPlot <- function(
    Task,
    Mu,
    Prob,
    OptThreshold,
    Hits,
    ylim = c(0, 1),
    xlim = c(min(Mu), max(Mu)),
    cex = 1.5,
    pch = 16,
    col = 8,
    bty = "n",
    xlab = "Mean expression",
    ylab = paste(Task, "probability"),
    title = ""
  ) {

  df <- data.frame(Mu, Prob)
  ggplot2::ggplot(df, ggplot2::aes_string(x = "Mu", y = "Prob")) +
    ggplot2::geom_point(
      ggplot2::aes(color = ifelse(Hits, Task, "Other")), 
      pch = pch, cex = cex) +
    ggplot2::scale_color_brewer(palette = "Set1", name = "") +
    ggplot2::geom_hline(
      yintercept = OptThreshold[[1]], lty = 2, col = "black"
    ) +
    ggplot2::scale_x_log10() +
    ggplot2::labs(
      x = xlab,
      y = ylab,
      title = title
    )
}

.VG <- function(x) {
  x@Table[[x@Name]]
}

.HeaderDetectHVG_LVG <- function(Chain,
                           PercentileThreshold,
                           VarThreshold,
                           ProbThreshold,
                           EFDR,
                           Plot) {

  if (!is(Chain, "BASiCS_Chain")) {
    stop("'Chain' is not a BASiCS_Chain class object.")
  }
  
  # Test if the chain contains beta parameters
  if (is.null(Chain@parameters$beta) & is.null(VarThreshold)) {
    stop(
      "'Chain' was not generated using the BASiCS regression model.
      Please supply a values between 0 and 1 for 'VarThreshold'."
    )
  }
  
  # Add a warning that by default, the variance decomposition threshold is
  # used if the user does not supply the PercentileThreshold parameter
  if (is.null(Chain@parameters$beta) & !is.null(PercentileThreshold)) {
    warning(
      "'Chain' was not generated using the BASiCS regression model.\n", 
      "By default, variable genes are detected by testing against ",
      "a variance threshold of ", 100 * VarThreshold, "%"
    )
  }
  
  # Add a warning that it's better to use the regression trend when estimating 
  # highly variable genes 
  # if (!is.null(Chain@parameters$beta) & !is.null(PercentileThreshold)) {
  #   warning(
  #     "'Chain' was generated using the BASiCS regression model.\n",
  #     "By default, the ", 100 * PercentileThreshold, 
  #     " percentile of variable genes will be returned."
  #   )
  # }
  
  if (!is.null(VarThreshold)){ 
    if (VarThreshold < 0 | VarThreshold > 1 | !is.finite(VarThreshold)) {
      stop("Variance contribution threshold must be in (0,1)")
    }
  }
  
  if (!is.null(PercentileThreshold)){
    if(PercentileThreshold < 0 | PercentileThreshold > 1 | 
       !is.finite(PercentileThreshold)) {
      stop("Percentile threshold must be in (0,1)")
    }
  }
  
  if (!is.logical(Plot) | length(Plot) != 1) {
    stop("Please insert `TRUE` or `FALSE` for `Plot` parameter")
  }
  if (!is.null(ProbThreshold)) {
    if (ProbThreshold < 0 | ProbThreshold > 1 | !is.finite(ProbThreshold)) {
      stop(
        "Posterior probability threshold must be contained in (0,1) \n
        For automatic threshold search use `ProbThreshold = NULL`."
      )
    } else {
      message(
        "Posterior probability threshold has been supplied \n",
        "EFDR will not be calibrated. \n"
      )
    }
  }
  if (is.null(ProbThreshold)) {
    message("Posterior probability threshold set by EFDR = ",
            100 * EFDR, "% (+-2.5% tolerance) ...")
  }
}
