#' @name BASiCS_DenoisedRates
#'
#' @title Calculates denoised expression rates
#'
#' @description Calculates normalised and denoised expression rates, by 
#' removing the effect of technical variation.
#'
#' @param Data an object of class 
#' \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#' @param Chain an object of class \code{\link[BASiCS]{BASiCS_Chain}}
#' @param PrintProgress If \code{TRUE}, partial progress 
#' information is printed in the console.
#' @param Propensities If \code{TRUE}, returns underlying 
#' expression propensitites \eqn{\rho_{ij}}. 
#' Otherwise, denoised rates \eqn{\mu_i \rho_{ij}} are returned.
#'
#' @examples
#'
#' # See
#' help(BASiCS_MCMC)
#'
#' @details See vignette
#'
#' @return A matrix of denoised expression rates.
#'
#' @seealso \code{\link[BASiCS]{BASiCS_Chain}}
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
                                Propensities = FALSE) {
    if (!is(Data, "SingleCellExperiment")) 
        stop("'Data' is not a SingleCellExperiment class object.")
    if (!is(Chain, "BASiCS_Chain")) 
        stop("'Chain' is not a BASiCS_Chain class object.")
    
    N = dim(Chain@delta)[1]
    q.bio = dim(Chain@delta)[2]
    n = dim(Chain@phi)[2]
    
    print(paste("This calculation requires a loop across the", N, 
                "MCMC iterations"))
    print("Please be patient ... ")
    cat("\n")
    print("To see a progress report use PrintProgress = TRUE")
    cat("\n")
    
    rho = matrix(0, ncol = n, nrow = q.bio)
    for (m in 1:N) {
        if (PrintProgress) {
            print(paste("Iteration", m, "has been completed."))
        }
        aux1 = assay(Data)[1:q.bio, ] + 1/Chain@delta[m, ]
        aux2 = t(tcrossprod(Chain@phi[m, ] * Chain@nu[m, ], 
                            Chain@mu[m, 1:q.bio])) + 1/Chain@delta[m, ]
        rho = rho + aux1/aux2
    }
    rho = rho/N
    
    if (Propensities) {
        out = rho
    } else {
        out = rho * apply(Chain@mu[, 1:q.bio], 2, median)
    }
    
    return(out)
    
}
