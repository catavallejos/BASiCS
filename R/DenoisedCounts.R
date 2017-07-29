#' @name DenoisedCounts
#'
#' @title Calculates normalised and denoised expression expression counts
#'
#' @description Calculates normalised and denoised expression counts, by removing the effect of technical variation.
#'
#' @param Data an object of class \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#' @param Chain an object of class \code{\link[BASiCS]{BASiCS_Chain-class}}
#'
#' @examples
#'
#' # See
#' help(BASiCS_MCMC)
#'
#' @details See vignette
#'
#' @return A matrix of normalised and denoised expression counts.
#'
#' @seealso \code{\link[BASiCS]{BASiCS_Chain-class}}
#'
#' @author Catalina A. Vallejos \email{cnvallej@@uc.cl}
#'
#' @rdname DenoisedCounts
BASiCS_DenoisedCounts=function(
  Data,
  Chain)
{
  if(!is(Data,"SummarizedExperiment")) stop("'Data' is not a SummarizedExperiment class object.")
  if(!is(Chain,"BASiCS_Chain")) stop("'Chain' is not a BASiCS_Chain class object.")
  
  N=dim(Chain@delta)[1]; q.bio=dim(Chain@delta)[2]; n=dim(Chain@phi)[2]
  
  Phi = apply(Chain@phi,2,median)
  Nu = apply(Chain@nu,2,median)
  
  out1 = t( t(assay(Data)[!rowData(Data)$Tech,]) / (Phi * Nu) )
  out2 = t( t(assay(Data)[rowData(Data)$Tech,]) / Nu )
  out = rbind(out1, out2)
  
  return(out)
}
