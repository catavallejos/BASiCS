#' @name BASiCS_VarThresholdSearchHVG
#' @aliases BASiCS_VarThresholdSearchHVG BASiCS_VarThresholdSearchHVG_LVG
#'
#' @title Detection method for highly and lowly variable genes using a grid of variance contribution thresholds
#'
#' @description Detection method for highly and lowly variable genes using a grid of variance contribution thresholds
#'
#' @param Chain an object of class \code{\link[BASiCS]{BASiCS_Chain-class}}
#' @param VarThresholdsGrid Grid of values for the variance contribution threshold (they must be contained in (0,1))
#' @param EFDR Target for expected false discovery rate related to HVG/LVG detection. Default: \code{EFDR = 0.10}.
#' @param Progress If \code{Progress = TRUE}, partial output is printed in the console.
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
#' \item{\code{BASiCS_VarThresholdSearchHVG}}{A table displaying the results of highly variable genes detection for different variance contribution thresholds.}
#' \item{\code{BASiCS_VarThresholdSearchLVG}}{A table displaying the results of lowly variable genes detection for different variance contribution thresholds.}
#' }
#'
#'
#' @seealso \code{\link[BASiCS]{BASiCS_Chain-class}}
#'
#' @author Catalina A. Vallejos \email{cnvallej@@uc.cl}
#' 
#' @references 
#' Vallejos, Marioni and Richardson (2015). Bayesian Analysis of Single-Cell Sequencing data. PLoS Computational Biology. 
#'
#' @rdname BASiCS_VarThresholdSearchHVG_LVG
BASiCS_VarThresholdSearchHVG=function(
  Chain,
  VarThresholdsGrid, #
  EFDR = 0.10,
  Progress = TRUE)
{
  if(!is(Chain,"BASiCS_Chain")) stop("'object' is not a BASiCS_Chain class object.")
  if(sum(VarThresholdsGrid<0)>0 | sum(VarThresholdsGrid>1)>0 | sum(!is.finite(VarThresholdsGrid))>0 )
    stop("Variance contribution thresholds for HVG and LVG detection must be contained in (0,1).")
  
  Table=matrix(0,nrow=length(VarThresholdsGrid),ncol=5)
  colnames(Table)=c("Var. Threshold (%)","EFDR (%)", "EFNR (%)","Optimal evidence thres.","# Detected genes")
  
  for(i in 1:length(VarThresholdsGrid))
  {
    VarThreshold=VarThresholdsGrid[i]
    
    if(Progress) {message(paste0("Evaluating variance contribution threshold = ",100*VarThreshold," % ... \n"))}
    
    suppressMessages(DetectHVG <- BASiCS_DetectHVG(Chain, EFDR = EFDR, VarThreshold = VarThreshold))
    
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
  Chain,
  VarThresholdsGrid, # Range of values for the variance contribution threshold (they must be contained in (0,1))
  EFDR = 0.10,
  Progress = TRUE)
{
  
  if(!is(Chain,"BASiCS_Chain")) stop("'object' is not a BASiCS_Chain class object.")
  if(sum(VarThresholdsGrid<0)>0 | sum(VarThresholdsGrid>1)>0 | sum(!is.finite(VarThresholdsGrid))>0 )
    stop("Variance contribution thresholds for HVG and LVG detection must be contained in (0,1).")
  
  Table=matrix(0,nrow=length(VarThresholdsGrid),ncol=5)
  colnames(Table)=c("VarThres (%)","EFDR (%)", "EFNR (%)","Optimal evi thres","# Detected genes")
  
  for(i in 1:length(VarThresholdsGrid))
  {
    VarThreshold=VarThresholdsGrid[i]
    
    if(Progress) {message(paste0("Evaluating variance contribution threshold = ",100*VarThreshold," % ... \n"))}
    
    suppressMessages(DetectLVG <- BASiCS_DetectLVG(Chain, EFDR = EFDR, VarThreshold = VarThreshold))
    
    Table[i,]=c(100*VarThreshold, round(100*DetectLVG$EFDR,2),round(100*DetectLVG$EFNR,2),
                DetectLVG$EviThreshold,sum(as.numeric(DetectLVG$Table[,7])))
    
  }
  return(Table)
}