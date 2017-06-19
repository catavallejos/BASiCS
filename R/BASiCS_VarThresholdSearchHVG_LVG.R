#' @name BASiCS_VarThresholdSearchHVG
#' @aliases BASiCS_VarThresholdSearchHVG BASiCS_VarThresholdSearchHVG_LVG
#'
#' @title Detection method for highly and lowly variable genes using a grid of variance contribution thresholds
#'
#' @description Detection method for highly and lowly variable genes using a grid of variance contribution thresholds
#'
#' @param Data an object of class \code{\link[BASiCS]{BASiCS_Data-class}}
#' @param object an object of class \code{\link[BASiCS]{BASiCS_Chain-class}}
#' @param VarThresholdsGrid Grid of values for the variance contribution threshold (they must be contained in (0,1))
#' @param PrintProgress If \code{PrintProgress = TRUE}, partial output is printed in the console.
#'
#' @examples
#'
#' # See
#' help(BASiCS_MCMC)
#'
#' @details See vignette
#'
#' @return
#' \describe{
#' \item{\code{BASiCS_VarThresholdSearchHVG}}{A table displaying the results of highly variable genes detecting for different variance contribution thresholds.}
#' \item{\code{BASiCS_VarThresholdSearchLVG}}{A table displaying the results of lowly variable genes detecting for different variance contribution thresholds.oo}
#' }
#'
#'
#' @seealso \code{\link[BASiCS]{BASiCS_Chain-class}}
#'
#' @author Catalina A. Vallejos \email{cnvallej@@uc.cl}
#'
#' @rdname BASiCS_VarThresholdSearchHVG_LVG

# Change Data class to SE class

BASiCS_VarThresholdSearchHVG=function(
  Data,
  object,
  VarThresholdsGrid, #
  PrintProgress = FALSE)
{
  
  if(!is(object,"BASiCS_Chain")) stop("'object' is not a BASiCS_Chain class object.")
  if(sum(VarThresholdsGrid<0)>0 | sum(VarThresholdsGrid>1)>0 | sum(!is.finite(VarThresholdsGrid))>0 )
    stop("Variance contribution thresholds for HVG and LVG detection must be contained in (0,1).")
  
  Table=matrix(0,nrow=length(VarThresholdsGrid),ncol=5)
  colnames(Table)=c("Var. Threshold (%)","EFDR (%)", "EFNR (%)","Optimal evidence thres.","# Detected genes")
  
  for(i in 1:length(VarThresholdsGrid))
  {
    VarThreshold=VarThresholdsGrid[i]
    
    if(PrintProgress) {print(paste0("Evaluating variance contribution threshold = ",100*VarThreshold," % ...")); cat("\n")}
    
    DetectHVG <- BASiCS_DetectHVG(Data, object, VarThreshold = VarThreshold)
    
    Table[i,]=c(100*VarThreshold, round(100*DetectHVG$EFDR,2),
                round(100*DetectHVG$EFNR,2),
                DetectHVG$EviThreshold,sum(as.numeric(DetectHVG$Table[,7])))
  }
  return(Table)
}

#' @name BASiCS_VarThresholdSearchLVG
#' @aliases BASiCS_VarThresholdSearchLVG BASiCS_VarThresholdSearchHVG_LVG
#' @rdname BASiCS_VarThresholdSearchHVG_LVG
BASiCS_VarThresholdSearchLVG=function(
  Data,
  object,
  VarThresholdsGrid, # Range of values for the variance contribution threshold (they must be contained in (0,1))
  PrintProgress = FALSE)
{
  
  if(!is(object,"BASiCS_Chain")) stop("'object' is not a BASiCS_Chain class object.")
  if(sum(VarThresholdsGrid<0)>0 | sum(VarThresholdsGrid>1)>0 | sum(!is.finite(VarThresholdsGrid))>0 )
    stop("Variance contribution thresholds for HVG and LVG detection must be contained in (0,1).")
  
  Table=matrix(0,nrow=length(VarThresholdsGrid),ncol=5)
  colnames(Table)=c("VarThres (%)","EFDR (%)", "EFNR (%)","Optimal evi thres","# Detected genes")
  
  for(i in 1:length(VarThresholdsGrid))
  {
    VarThreshold=VarThresholdsGrid[i]
    
    if(PrintProgress) {print(paste0("Evaluating variance contribution threshold = ",100*VarThreshold," % ...")); cat("\n")}
    
    DetectLVG <- BASiCS_DetectLVG(Data, object, VarThreshold = VarThreshold)
    
    Table[i,]=c(100*VarThreshold, round(100*DetectLVG$EFDR,2),round(100*DetectLVG$EFNR,2),
                DetectLVG$EviThreshold,sum(as.numeric(DetectLVG$Table[,7])))
    
  }
  return(Table)
}