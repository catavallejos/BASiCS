#' @name BASiCS_VarThresholdSearchHVG
#' @aliases BASiCS_VarThresholdSearchHVG BASiCS_VarThresholdSearchHVG_LVG
#'
#' @title Detection method for highly and lowly variable genes using
#' a grid of variance contribution thresholds
#'
#' @description Detection method for highly and lowly variable genes
#' using a grid of variance contribution thresholds. Only used when 
#' HVG/LVG are found based on the variance decomposition. 
#'
#' @param Chain an object of class \code{\linkS4class{BASiCS_Chain}}
#' @param Task See \code{?BASiCS_DetectVG}.
#' @param VarThresholdsGrid Grid of values for the variance contribution
#' threshold (they must be contained in (0,1))
#' @param EFDR Target for expected false discovery rate related to
#' HVG/LVG detection. Default: \code{EFDR = 0.10}.
#' @param Progress If \code{Progress = TRUE}, partial output is
#' printed in the console. Default: \code{Progress = TRUE}.
#' @param ... Passed to methods.
#' @examples
#'
#' data(ChainSC)
#' 
#' BASiCS_VarThresholdSearchHVG(ChainSC,
#'                              VarThresholdsGrid = seq(0.55,0.65,by=0.01),
#'                              EFDR = 0.10)
#' BASiCS_VarThresholdSearchLVG(ChainSC,
#'                              VarThresholdsGrid = seq(0.35,0.45,by=0.01),
#'                              EFDR = 0.10)
#'
#' @details See vignette
#'
#' @return
#' \describe{
#' \item{\code{BASiCS_VarThresholdSearchHVG}}{A table displaying the results of
#'       highly variable genes detection for different variance
#'       contribution thresholds.}
#' \item{\code{BASiCS_VarThresholdSearchLVG}}{A table displaying the results of
#'       lowly variable genes detection for different variance
#'       contribution thresholds.}
#' }
#'
#'
#' @seealso \code{\linkS4class{BASiCS_Chain}}
#'
#' @author Catalina A. Vallejos \email{cnvallej@@uc.cl}
#'
#' @references
#' Vallejos, Marioni and Richardson (2015). PLoS Computational Biology.
#'
#' @rdname BASiCS_VarThresholdSearchHVG_LVG
#' @export
BASiCS_VarThresholdSearchVG <- function(Chain,
                                        Task = c("HVG", "LVG"),
                                        VarThresholdsGrid,
                                        EFDR = 0.1,
                                        Progress = TRUE) {

  if (!is(Chain, "BASiCS_Chain")) {
    stop("'object' is not a BASiCS_Chain class object.")
  }
  if (sum(VarThresholdsGrid < 0) > 0 |
      sum(VarThresholdsGrid > 1) > 0 |
      sum(!is.finite(VarThresholdsGrid)) > 0) {
    stop("Variance contribution thresholds must be contained in (0,1).")
  }

  Table <- as.data.frame(matrix(0, nrow = length(VarThresholdsGrid), ncol = 5))
  colnames(Table) <- c("Var. Threshold (%)", "EFDR (%)", "EFNR (%)",
                       "Optimal evidence thres.", "# Detected genes")

  for (i in seq_along(VarThresholdsGrid)) {
    VarThreshold <- VarThresholdsGrid[i]
    if (Progress) {
      message(
        "Evaluating variance contribution threshold = ",
        100 * VarThreshold, " % ... \n")
    }

    suppressMessages(
      DetectHVG <- BASiCS_DetectVG(
        Chain,
        Task = Task,
        EFDR = EFDR,
        VarThreshold = VarThreshold
      )
    )

    Table[i, ] <- c(
      100 * VarThreshold,
      round(100 * DetectHVG@EFDR, 2),
      round(100 * DetectHVG@EFNR, 2),
      DetectHVG@ProbThreshold,
      sum(.VG(DetectHVG))
    )
  }
  return(Table)
}

#' @name BASiCS_VarThresholdSearchHVG
#' @aliases BASiCS_VarThresholdSearchLVG BASiCS_VarThresholdSearchHVG_LVG
#' @rdname BASiCS_VarThresholdSearchHVG_LVG
#' @export
BASiCS_VarThresholdSearchHVG <- function(...) {
  # .Deprecated("BASiCS_VarThresholdSearchVG")
  BASiCS_VarThresholdSearchVG(..., Task = "HVG")
}
#' @name BASiCS_VarThresholdSearchLVG
#' @aliases BASiCS_VarThresholdSearchLVG BASiCS_VarThresholdSearchHVG_LVG
#' @rdname BASiCS_VarThresholdSearchHVG_LVG
#' @export
BASiCS_VarThresholdSearchLVG <- function(...) {
  # .Deprecated("BASiCS_VarThresholdSearchVG")
  BASiCS_VarThresholdSearchVG(..., Task = "LVG")
}
