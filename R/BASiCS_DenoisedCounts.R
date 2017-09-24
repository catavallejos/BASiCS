#' @name BASiCS_DenoisedCounts
#'
#' @title Calculates denoised expression expression counts
#'
#' @description Calculates denoised expression counts by removing 
#' cell-specific technical variation.
#'
#' @param Data an object of class \code{\linkS4class{SingleCellExperiment}}
#' @param Chain an object of class \code{\linkS4class{BASiCS_Chain}} 
#'
#' @examples
#'
#' # See
#' help(BASiCS_MCMC)
#'
#' @details See vignette
#'
#' @return A matrix of denoised expression counts. In line with global scaling
#' normalisation strategies, these are defined as \eqn{X_{ij}/(\phi_j \nu_j)} 
#' for biological genes and \eqn{X_{ij}/(\nu_j)} for spike-in genes. For this
#' calculation \eqn{\phi_j} \eqn{\nu_j} are estimated by their corresponding
#' posterior medians.
#'
#' @seealso \code{\linkS4class{BASiCS_Chain}} 
#'
#' @author Catalina A. Vallejos \email{cnvallej@@uc.cl}
#' 
#' @references 
#' 
#' Vallejos, Marioni and Richardson (2015). PLoS Computational Biology. 
#'
#' @rdname BASiCS_DenoisedCounts
BASiCS_DenoisedCounts = function(Data, Chain) 
{
    if (!is(Data, "SingleCellExperiment")) 
        stop("'Data' is not a SingleCellExperiment class object.")
    if (!is(Chain, "BASiCS_Chain")) 
        stop("'Chain' is not a BASiCS_Chain class object.")
    
    N <- nrow(Chain@delta)
    q.bio <- ncol(Chain@delta)
    n <- nrow(Chain@phi)
    
    Phi <- matrixStats::colMedians(Chain@phi)
    Nu <- matrixStats::colMedians(Chain@nu)
    
    out1 <- t(t(assay(Data)[!isSpike(Data), ])/(Phi * Nu))
    out2 <- t(t(assay(Data)[isSpike(Data), ])/Nu)
    out <- rbind(out1, out2)
    
    return(out)
}
