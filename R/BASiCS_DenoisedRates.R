#' @name BASiCS_DenoisedRates
#'
#' @title Calculates denoised expression rates
#'
#' @description Calculates normalised and denoised expression rates, by 
#' removing the effect of technical variation.
#'
#' @param Data an object of class \code{\linkS4class{SingleCellExperiment}} 
#' @param Chain an object of class \code{\linkS4class{BASiCS_Chain}}
#' @param PrintProgress If \code{TRUE}, partial progress 
#' information is printed in the console. Default: \code{PrintProgress = FALSE}.
#' @param Propensities If \code{TRUE}, returns underlying 
#' expression propensitites \eqn{\rho_{ij}}. 
#' Otherwise, denoised rates \eqn{\mu_i \rho_{ij}} are returned.
#' Default: \code{Propensities = FALSE}. 
#'
#' @examples
#'
#' # See
#' help(BASiCS_MCMC)
#'
#' @details See vignette
#'
#' @return A matrix of denoised expression rates (biological genes only)
#'
#' @seealso \code{\linkS4class{BASiCS_Chain}}
#'
#' @author Catalina A. Vallejos \email{cnvallej@@uc.cl}
#'
#' @references 
#' 
#' Vallejos, Marioni and Richardson (2015). PLoS Computational Biology. 
#'
#' @rdname BASiCS_DenoisedRates
BASiCS_DenoisedRates = function(Data, 
                                Chain, 
                                PrintProgress = FALSE, 
                                Propensities = FALSE) 
{
  if (!is(Data, "SingleCellExperiment")) 
    stop("'Data' is not a SingleCellExperiment class object.")
  if (!is(Chain, "BASiCS_Chain")) 
    stop("'Chain' is not a BASiCS_Chain class object.")
    
  N <- nrow(Chain@delta)
  q.bio <- ncol(Chain@delta)
  n <- ncol(Chain@phi)

  # Message no longer required as calculations are faster
#  message("This calculation requires a loop across the", N, "MCMC iterations\n",
#          "Please be patient ... \n",
#          "To see a progress report use PrintProgress = TRUE. \n")
  
  CountsBio <- assay(Data)[!isSpike(Data),]
  Rho <- HiddenBASiCS_DenoisedRates(CountsBio, 
                                    Chain@mu, Chain@delta,
                                    Chain@phi, Chain@nu,
                                    N, q.bio, n)

  if (Propensities) { out <- Rho } 
  else { out <- Rho * matrixStats::colMedians(Chain@mu) }
    
  return(out)
    
}
