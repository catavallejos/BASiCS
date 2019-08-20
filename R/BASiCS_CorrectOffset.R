#' @name BASiCS_TestDE
#'
#' @title Remove global mean expression offset
#'
#' @description Remove global offset in mean expression between two 
#' \code{BASiCS_Chain} objects. 
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
#' @return A list whose first element is an offset corrected version of `Chain1` 
#' (using `Chain2` as a reference), whose second element is the point estimate
#' for the offset and whose third element contains iteration-specific offsets.  
#'
#' @author Catalina A. Vallejos \email{cnvallej@@uc.cl}
#' @author Nils Eling \email{eling@@ebi.ac.uk}
#' @author Alan O'Callaghan \email{a.b.o'callaghan@sms.ed.ac.uk}
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
  include <- which((colMedians(mu1) / OffsetEst0 + colMedians(mu2))/2 >= min.mean)
  
  # Calculating iteration-specific offset 
  OffsetChain <- matrixStats::rowMedians(mu1[,include]) / 
    matrixStats::rowMedians(mu2[,include])
  # Offset point estimate
  OffsetEst <- median(OffsetChain)
  
  # Application of offset
  Chain_offset <- Chain
  Chain_offset@parameters$mu <- Chain@parameters$mu / OffsetEst
  Chain_offset@parameters$phi <- Chain@parameters$phi * OffsetEst
  
  list("Chain" = Chain_offset, 
       "Offset" = OffsetEst, 
       "OffsetChain" = OffsetChain)

}