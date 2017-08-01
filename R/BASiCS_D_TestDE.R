#' @name BASiCS_D_TestDE
#' 
#' @title Detection of genes with changes in expression.
#' 
#' @description Function to assess changes in expression (mean and over-dispersion). This function is no longer in use and will be removed in future releases. Please use \code{\link[BASiCS]{BASiCS_TestDE}} instead.
#' 
#' @param Data an object of class \code{\link[BASiCS]{BASiCS_D_Data-class}}
#' @param object an object of class \code{\link[BASiCS]{BASiCS_D_Chain-class}}
#' @param GeneNames Vector containing gene names to be used in results table (argument to be removed as 'GeneNames' will be an slot of `BASiCS_D_Data` object)
#' @param EpsilonM Minimum fold change tolerance threshold for detecting changes in overall expression (must be a positive real number)
#' @param EpsilonD Minimum fold change tolerance threshold for detecting changes in cell-to-cell biological over dispersion (must be a positive real number)
#' @param EviThresholdM Optional parameter. Evidence threshold for detecting changes in overall expression (must be a positive value, between 0 and 1)
#' @param EviThresholdD Optional parameter. Evidence threshold for detecting changes in cell-to-cell biological over dispersion (must be a positive value, between 0 and 1)
#' @param OrderVariable Ordering variable for output. Must take values in \code{c("GeneIndex", "GeneNames", "ProbDiffExp", "ProbDiffOverDisp")}.
#' @param GroupLabelRef Label assigned to reference group. Default: \code{GroupLabelRef = "Ref"}
#' @param GroupLabelTest Label assigned to reference group. Default: \code{GroupLabelRef = "Test"}
#' @param Plot If \code{Plot = T}, error rates control rates and volcano plots are generated.  
#' @param OffSet Optional argument to remove a fix offset effect (if not previously removed from the MCMC chains). This argument will be removed shorly, once offset removal is built as an internal step. 
#' @param EFDR_M Target for expected false discovery rate related to the comparison of means (default = 0.05)
#' @param EFDR_D Target for expected false discovery rate related to the comparison of dispersions (default = 0.05)
#' @param GenesSelect Optional argument to provide a user-defined list of genes to be considered for the comparison (default = NULL). When used, this argument must be a vector of 'TRUE' (include gene) / 'FALSE' (exclude gene) indicator, with the same length as the number of intrinsic genes and following the same order as how genes are displayed in the table of counts.  This argument is necessary in order to have a meaningful EFDR calibration when the user decides to exclude some genes from the comparison. 
#' @param ... Graphical parameters (see \code{\link[graphics]{par}}).
#' 
#' @return This function is no longer in use and will be removed in future releases. Please use \code{\link[BASiCS]{BASiCS_TestDE}} instead.
#'  
#' @examples
#' 
#' # See vignette
#' 
#' @details See vignette
#' 
#' 
#' @seealso \code{\link[BASiCS]{BASiCS_D_Chain-class}} 
#' 
#' @author Catalina A. Vallejos \email{cnvallej@uc.cl}
#' 
#' @rdname BASiCS_D_TestDE
BASiCS_D_TestDE <- function(Data = NULL, 
                            object = NULL,
                            GeneNames = NULL,
                            EpsilonM = 0.10,
                            EpsilonD = 0.10,
                            EviThresholdM = NULL,
                            EviThresholdD = NULL,
                            OrderVariable = "ProbDiffExp",
                            GroupLabelRef = "Ref",
                            GroupLabelTest = "Test",
                            Plot = FALSE, 
                            OffSet = FALSE, 
                            EFDR_M = 0.05,
                            EFDR_D = 0.05,
                            GenesSelect = NULL, 
                            ...)
{
  message("This function is no longer in use and will be removed in future releases. \n
          Please use `BASiCS_TestDE` instead.")
}