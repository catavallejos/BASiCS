#' @title Simulates expression counts according to the model implemented in BASiCS
#'
#' @description \code{BASiCS_Sim} creates a simulated dataset from the model implemented in BASiCS.
#' This function is used in order to illustrate the performance of the \code{BASiCS} library.
#'
#' @param mu Gene-specific expression levels \eqn{\mu_i} for all biological genes 
#' (vector of length \code{q.bio}, all elements must be positive numbers)
#' @param mu_spikes \eqn{\mu_i} for all technical genes defined as true input molecules (SpikeInfo)
#' (vector of length \code{q-q.bio}, all elements must be positive numbers)
#' @param delta Gene-specific biological cell-to-cell over-dispersion parameters \eqn{\delta_i}, biological genes only
#' (vector of length \code{q.bio}, all elements must be positive numbers)
#' @param phi Cell-specific mRNA content normalising constants \eqn{\phi_j}
#' (vector of length \code{n}, all elements must be positive numbers and the sum of its elements must be equal to \code{n})
#' @param s Cell-specific capture efficiency (or amplification biases if not using UMI based counts) normalising constants \eqn{s_j}
#' (vector of length \code{n}, all elements must be positive numbers)
#' @param theta Technical variability hyper-parameter \eqn{\theta} (must be positive)
#'
#' @return An object of class \code{\link[SingleCellExperiment]{SingleCellExperiment}}, simulated from the model implemented in BASiCS.
#'
#' @examples
#'
#' # Simulated parameter values for 10 genes
#' # (7 biogical and 3 spike-in) measured in 5 cells
#' Mu =  c(8.36, 10.65, 4.88, 6.29, 21.72, 12.93, 30.19)
#' Mu_spike =  c(1010.72, 7.90, 31.59)
#' Delta = c(1.29, 0.88, 1.51, 1.49, 0.54, 0.40, 0.85)
#' Phi = c(1.00, 1.06, 1.09, 1.05, 0.80)
#' S = c(0.38, 0.40, 0.38, 0.39, 0.34)
#' Theta = 0.39
#'
#' Data = BASiCS_Sim(Mu, Mu_spike, Delta, Phi, S, Theta)
#' head(assay(Data))
#' dim(assay(Data))
#' metadata(Data)$SpikeInput
#' isSpike(Data)
#'
#' @author Catalina A. Vallejos \email{cnvallej@@uc.cl}, Nils Eling
#'
#' @references 
#' 
#' Vallejos, Marioni and Richardson (2015). Bayesian Analysis of Single-Cell Sequencing data. PLoS Computational Biology. 
BASiCS_Sim <- function(mu, mu_spikes, delta, phi, s, theta) {
    # Number of cells
    n = length(phi)
    # Total number of genes, including biological and technical ones
    q = length(mu) + length(mu_spikes)
    # Number of biological genes
    q.bio = length(delta)
    
    # Arguments checking
    if (!(is.vector(mu) & is.numeric(mu) & all(mu > 0))) 
        stop("Invalid argument value for 'mu'.")
    if (!(is.vector(mu_spikes) & is.numeric(mu_spikes) & all(mu_spikes > 0))) 
        stop("Invalid argument value for 'mu_spikes'.")
    if (!(is.vector(delta) & is.numeric(delta) & all(delta >= 0))) 
        stop("Invalid argument value for 'delta'.")
    if (!(is.vector(phi) & is.numeric(phi) & all(phi > 0) & isTRUE(all.equal(sum(phi), n, tolerance = 0.01)))) 
        stop("Invalid argument value for 'phi'.")
    if (!(is.vector(s) & is.numeric(s) & all(s > 0) & length(s) == n)) 
        stop("Invalid argument value for 's'.")
    if (!(is.numeric(theta) & length(theta) == 1 & theta >= 0)) 
        stop("Invalid argument value for 'theta'.")
    
    if (!all(c(n, q, q.bio, q - q.bio) > 0)) 
        stop("Arguments' dimensions are not compatible")
    
    # Merge biological and technical genes
    mu = c(mu, mu_spikes)
    
    # Matrix where simulated counts will be stored
    Counts.sim <- matrix(0, ncol = n, nrow = q)
    # Matrix where gene-cell specific simulated random effects will be stored
    rho <- matrix(1, ncol = n, nrow = q.bio)
    # Simulated cell-specific random effects
    if (theta > 0) {
        nu <- rgamma(n, shape = 1/theta, rate = 1/(s * theta))
    } else {
        nu <- s
    }
    # Simulated counts data
    for (i in 1:q) {
        # Biological genes
        if (i <= q.bio) {
            if (delta[i] > 0) {
                rho[i, ] <- rgamma(n, shape = 1/delta[i], rate = 1/delta[i])
            }
            Counts.sim[i, ] <- rpois(n, lambda = phi * nu * rho[i, ] * mu[i])
        } else {
            # Technical genes
            Counts.sim[i, ] <- rpois(n, lambda = nu * mu[i])
        }
    }
    
    rownames(Counts.sim) <- paste0("Gene", 1:q)
    SpikeInfo = data.frame(paste0("Gene", (q.bio + 1):q), mu[(q.bio + 1):q])
    Data = newBASiCS_Data(Counts = Counts.sim, Tech = ifelse(1:q > q.bio, TRUE, FALSE), SpikeInfo = SpikeInfo)
    
    return(Data)
}
