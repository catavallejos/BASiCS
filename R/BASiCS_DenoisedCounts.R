#' @name BASiCS_DenoisedCounts
#'
#' @title Calculates denoised expression expression counts
#'
#' @description Calculates denoised expression counts by removing 
#' cell-specific technical variation. The latter includes global-scaling
#' normalisation and therefore no further normalisation is required.
#'
#' @param Data an object of class \code{\linkS4class{SingleCellExperiment}}
#' @param Chain an object of class \code{\linkS4class{BASiCS_Chain}} 
#'
#' @examples
#'
#' Data <- makeExampleBASiCS_Data(WithSpikes = TRUE)
#' ## The N and Burn parameters used here are optimised for speed
#' ## and should not be used in regular use.
#' ## For more useful parameters,
#' ## see the vignette (\code{browseVignettes("roxygen2")})
#' Chain <- BASiCS_MCMC(Data, N = 1000, Thin = 10, Burn = 500, 
#'                      Regression = FALSE, PrintProgress = FALSE)
#'
#' DC <- BASiCS_DenoisedCounts(Data, Chain)
#'
#' @details See vignette \code{browseVignettes("roxygen2")}
#'
#' @return A matrix of denoised expression counts. In line with global scaling
#' normalisation strategies, these are defined as \eqn{X_{ij}/(\phi_j \nu_j)} 
#' for biological genes and \eqn{X_{ij}/(\nu_j)} for spike-in genes. For this
#' calculation \eqn{\phi_j} \eqn{\nu_j} are estimated by their corresponding
#' posterior medians. If spike-ins are not used, \eqn{\phi_j} is set equal to 1. 
#'
#' @seealso \code{\linkS4class{BASiCS_Chain}} 
#'
#' @author Catalina A. Vallejos \email{cnvallej@@uc.cl}
#' @author Nils Eling \email{eling@@ebi.ac.uk}
#'
#' @rdname BASiCS_DenoisedCounts
#' @export
BASiCS_DenoisedCounts <- function(Data, Chain) 
{
    if (!is(Data, "SingleCellExperiment")) 
        stop("'Data' is not a SingleCellExperiment class object.")
    if (!is(Chain, "BASiCS_Chain")) 
        stop("'Chain' is not a BASiCS_Chain class object.")
    
    Nu <- matrixStats::colMedians(Chain@parameters$nu)
    if("phi" %in% names(Chain@parameters))
    {
      # Spikes case
      Phi <- matrixStats::colMedians(Chain@parameters$phi)
      out1 <- t(t(counts(Data)[!isSpike(Data), ])/(Phi * Nu))
      out2 <- t(t(counts(Data)[isSpike(Data), ])/Nu) 
      out <- rbind(out1, out2) 
    }
    else
    {
      # No spikes case
      out <- t(t(counts(Data))/Nu)
    }

    rownames(out) <- rownames(Data)
    colnames(out) <- colnames(Data)
    
    return(out)
}
