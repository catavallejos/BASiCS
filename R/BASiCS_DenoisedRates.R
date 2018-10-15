#' @name BASiCS_DenoisedRates
#'
#' @title Calculates denoised expression rates
#'
#' @description Calculates normalised and denoised expression rates, by 
#' removing the effect of technical variation.
#'
#' @param Data an object of class \code{\linkS4class{SingleCellExperiment}} 
#' @param Chain an object of class \code{\linkS4class{BASiCS_Chain}}
#' @param Propensities If \code{TRUE}, returns underlying 
#' expression propensitites \eqn{\rho_{ij}}. 
#' Otherwise, denoised rates \eqn{\mu_i \rho_{ij}} are returned.
#' Default: \code{Propensities = FALSE}. 
#'
#' @examples
#'
#' Data <- makeExampleBASiCS_Data(WithSpikes = TRUE)
#' ## The N and Burn parameters used here are optimised for speed
#' ## and should not be used in regular use.
#' ## For more useful parameters,
#' ## see the vignette (\code{browseVignettes("BASiCS")})
#' Chain <- BASiCS_MCMC(Data, N = 1000, Thin = 10, Burn = 500, 
#'                      Regression = FALSE, PrintProgress = FALSE)
#'
#' DR <- BASiCS_DenoisedRates(Data, Chain)
#'
#' @details See vignette
#'
#' @return A matrix of denoised expression rates (biological genes only)
#'
#' @seealso \code{\linkS4class{BASiCS_Chain}}
#'
#' @author Catalina A. Vallejos \email{cnvallej@@uc.cl}
#' @author Nils Eling \email{eling@@ebi.ac.uk}
#'
#' @rdname BASiCS_DenoisedRates
#' @export
BASiCS_DenoisedRates <- function(Data, Chain, Propensities = FALSE) 
{
  if (!is(Data, "SingleCellExperiment")) 
    stop("'Data' is not a SingleCellExperiment class object.")
  if (!is(Chain, "BASiCS_Chain")) 
    stop("'Chain' is not a BASiCS_Chain class object.")
    
  N <- nrow(Chain@parameters$delta)
  q.bio <- ncol(Chain@parameters$delta)
  n <- ncol(Chain@parameters$s)

  if("phi" %in% names(Chain@parameters))
  {
    # Spikes case
    CountsBio <- counts(Data)[!SingleCellExperiment::isSpike(Data),]
    Rho <- HiddenBASiCS_DenoisedRates(CountsBio, 
                                      Chain@parameters$mu, 
                                      t(1/Chain@parameters$delta),
                                      Chain@parameters$phi*Chain@parameters$nu,
                                      N, q.bio, n)    
  }
  else
  {
    # No spikes case
    CountsBio <- counts(Data)
    Rho <- HiddenBASiCS_DenoisedRates(CountsBio, 
                                      Chain@parameters$mu, 
                                      t(1/Chain@parameters$delta),
                                      Chain@parameters$nu,
                                      N, q.bio, n)        
  }


  if (Propensities) { out <- Rho } 
  else { out <- Rho * matrixStats::colMedians(Chain@parameters$mu) }
  
  if("phi" %in% names(Chain@parameters)){
    rownames(out) <- rownames(Data)[!SingleCellExperiment::isSpike(Data)]
  }
  else{
    rownames(out) <- rownames(Data)
  }
  colnames(out) <- colnames(Data)
    
  return(out)
    
}
