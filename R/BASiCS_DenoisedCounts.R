#' @name BASiCS_DenoisedCounts
#'
#' @title Calculates denoised expression expression counts
#'
#' @description Calculates denoised expression counts by removing 
#' cell-specific technical variation.
#'
#' @param Data an object of class 
#' \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#' @param Chain an object of class \code{\link[BASiCS]{BASiCS_Chain}}
#'
#' @examples
#'
#' # See
#' help(BASiCS_MCMC)
#'
#' @details See vignette
#'
#' @return A matrix of denoised expression counts. These are defined as 
#' \eqn{X_{ij}/(\phi_j \nu_j)} for biological genes and 
#' \eqn{X_{ij}/(\nu_j)} for spike-in genes.
#'
#' @seealso \code{\link[BASiCS]{BASiCS_Chain}}
#'
#' @author Catalina A. Vallejos \email{cnvallej@@uc.cl}
#' 
#' @references 
#' 
#' Vallejos, Marioni and Richardson (2015). PLoS Computational Biology. 
#'
#' @rdname BASiCS_DenoisedCounts
BASiCS_DenoisedCounts = function(Data, Chain) {
    if (!is(Data, "SingleCellExperiment")) 
        stop("'Data' is not a SingleCellExperiment class object.")
    if (!is(Chain, "BASiCS_Chain")) 
        stop("'Chain' is not a BASiCS_Chain class object.")
    
    N = dim(Chain@delta)[1]
    q.bio = dim(Chain@delta)[2]
    n = dim(Chain@phi)[2]
    
    Phi = apply(Chain@phi, 2, median)
    Nu = apply(Chain@nu, 2, median)
    
    out1 = t(t(assay(Data)[!isSpike(Data), ])/(Phi * Nu))
    out2 = t(t(assay(Data)[isSpike(Data), ])/Nu)
    out = rbind(out1, out2)
    
    return(out)
}
