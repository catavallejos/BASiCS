#' @name BASiCS_CorrectOffset
#'
#' @title Remove global mean expression offset
#'
#' @description Remove global offset in mean expression between two 
#' \code{BASiCS_Chain} objects. 
#' 
#' @param Chain a `BASiCS_MCMC` object to which the offset correction should 
#' be applied (with respect to `ChainRef`).
#' @param ChainRef a `BASiCS_MCMC` object to be used as the reference in the
#' offset correction procedure.
#' @param min.mean Minimum mean expression threshold required for inclusion in
#' offset calculation. Similar to `min.mean` in `scran::computeSumFactors`. 
#'
#' @examples
#'
#' # Loading two 'BASiCS_Chain' objects (obtained using 'BASiCS_MCMC')
#' data(ChainSC)
#' data(ChainRNA)
#' 
#' A <- BASiCS_CorrectOffset(ChainSC, ChainRNA)
#' 
#' # Offset corrected versions for ChainSC (with respect to ChainRNA). 
#' A$Chain
#' A$Offset
#' 
#' @return A list whose first element is an offset corrected version of `Chain` 
#' (using `ChainRef` as a reference), whose second element is the point estimate
#' for the offset and whose third element contains iteration-specific offsets.  
#'
#' @author Catalina A. Vallejos \email{cnvallej@@uc.cl}
#' @author Nils Eling \email{eling@@ebi.ac.uk}
#' @author Alan O'Callaghan
#' 
#' @export
BASiCS_CorrectOffset <- function(Chain, 
                                 ChainRef,
                                 min.mean = 1) {

  # Extract MCMC chains for mean parameters
  mu1 <- Chain@parameters$mu
  mu2 <- ChainRef@parameters$mu

  # Lowly expressed genes are excluded from offset calculation
  # This is similar to what is done in scran:::.rescale_clusters
  # Rough offset estimate applied for this purpose
  # This is based on medians to be more robust (rowMeans2 used before)
  OffsetChain0 <- matrixStats::rowMedians(mu1) / matrixStats::rowMedians(mu2)
  OffsetEst0 <- median(OffsetChain0)
  OffsetRatio <- (colMedians(mu1) / OffsetEst0 + colMedians(mu2)) / 2
  include <- which(OffsetRatio >= min.mean)
  
  # Calculating iteration-specific offset 
  OffsetChain <- matrixStats::rowMedians(mu1[, include]) / 
    matrixStats::rowMedians(mu2[, include])
  # Offset point estimate
  OffsetEst <- median(OffsetChain)
  
  # Application of offset
  Chain_offset <- Chain
  Chain_offset@parameters$mu <- Chain@parameters$mu / OffsetEst
  if("phi" %in% names(Chain@parameters)) {
    Chain_offset@parameters$phi <- Chain@parameters$phi * OffsetEst  
  } else {
    Chain_offset@parameters$s <- Chain@parameters$s * OffsetEst
  }
  
  list(
    Chain = Chain_offset, 
    Offset = OffsetEst, 
    OffsetChain = OffsetChain
  )
}
