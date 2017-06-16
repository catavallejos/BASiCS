#' @name BASiCS_Chain-class
#' @aliases BASiCS_Chain BASiCS_Chain,BASiCS_Chain-class
#' 
#' @title The BASiCS_Chain class
#' 
#' @description Container of an MCMC sample of the BASiCS' model parameters (see Vallejos et al, 2015) as generated
#' by the function \code{\link[BASiCS]{BASiCS_MCMC}}. 
#' 
#' @slot mu MCMC chain for gene-specific expression levels \eqn{\mu[i]}, defined as true input molecules in case of technical genes 
#' (matrix with \code{q} columns, technical genes located at the end of the matrix, all elements must be positive numbers)
#' @slot delta MCMC chain for gene-specific biological cell-to-cell heterogeneity hyper-parameters \eqn{\delta[i]}, biological genes only 
#' (matrix with \code{q.bio} columns, all elements must be positive numbers)
#' @slot phi MCMC chain for cell-specific mRNA content normalising constants \eqn{\phi[j]}
#' (matrix with \code{n} columns, all elements must be positive numbers and the sum of its elements must be equal to \code{n})
#' @slot s MCMC chain for cell-specific capture efficiency (or amplification biases if not using UMI based counts) normalising constants \eqn{s[j]}
#' (matrix with \code{n} columns, all elements must be positive numbers)
#' @slot nu MCMC chain for cell-specific random effects \eqn{\nu[j]}
#' (matrix with \code{n} columns, all elements must be positive numbers)
#' @slot theta MCMC chain for technical variability hyper-parameter(s) \eqn{\theta} (matrix, all elements must be positive, each colum 
#' represents 1 batch)
#'   
#' @examples
#' 
#' # A BASiCS_Chain object created by the BASiCS_MCMC function.
#' Data = makeExampleBASiCS_Data()
#' MCMC_Output <- BASiCS_MCMC(Data, N = 100, Thin = 2, Burn = 2)
#' 
#' @author Catalina A. Vallejos \email{catalina.vallejos@@mrc-bsu.cam.ac.uk}
#' 
#' @references Vallejos, Marioni and Richardson (2015). Bayesian Analysis of Single-Cell Sequencing data. PLoS Computational Biology.
setClass("BASiCS_Chain",
         representation = representation(
           mu = "matrix",
           delta = "matrix",
           phi = "matrix",
           s = "matrix",
           nu = "matrix",
           theta = "matrix"),
         validity = function(object){
           errors <- character()
           
           if(length(object@mu)==0 | length(object@delta)==0 | 
              length(object@phi)==0 | length(object@s)==0 | 
              length(object@nu)==0 | length(object@theta)==0) {errors <-c(errors,"One or more slots are missing"); stop(errors)}   
           
           N = nrow(object@mu)
           n = ncol(object@phi)
           if(nrow(object@delta) != N |
              nrow(object@phi) != N | nrow(object@s) != N |
              nrow(object@nu) != N | nrow(object@theta) != N |
              ncol(object@mu) != ncol(object@delta) |
              ncol(object@s) != n | ncol(object@nu) != n) {errors <-c(errors,"Slots' dimensions are not compatible")}
           
           if(sum(!is.finite(object@mu)) + sum(!is.finite(object@delta)) +
              sum(!is.finite(object@phi)) + sum(!is.finite(object@s)) +
              sum(!is.finite(object@nu)) + sum(!is.finite(object@theta))) {errors <-c(errors,"One or more of the slots contains NAs or Infinite values")}
           
           if (length(errors) == 0) TRUE else errors
         }
)


#' @name BASiCS_Summary-class
#' @aliases BASiCS_Summary,BASiCS_Summary-class
#' 
#' @title The BASiCS_Summary class
#' 
#' @slot mu Posterior medians (first column), lower (second column) and upper (third column) limits of gene-specific expression levels \eqn{\mu[i]}.
#' @slot delta Posterior medians (first column), lower (second column) and upper (third column) limits of gene-specific biological cell-to-cell heterogeneity hyper-parameters \eqn{\delta[i]}, biological genes only 
#' @slot phi Posterior medians (first column), lower (second column) and upper (third column) limits of cell-specific mRNA content normalising constants \eqn{\phi[j]}
#' @slot s Posterior medians (first column), lower (second column) and upper (third column) limits of cell-specific capture efficiency (or amplification biases if not using UMI based counts) normalising constants \eqn{s[j]}
#' @slot nu Posterior medians (first column), lower (second column) and upper (third column) limits of cell-specific random effects \eqn{\nu[j]}
#' @slot theta Posterior median (first column), lower (second column) and upper (third column) limits of technical variability hyper-parameter \eqn{\theta} (each row represents one batch)
#' 
#' @examples
#' 
#' # A BASiCS_Summary object created by the Summary method.
#' Data = makeExampleBASiCS_Data()
#' MCMC_Output <- BASiCS_MCMC(Data, N = 100, Thin = 2, Burn = 2)
#' MCMC_Summary <- Summary(MCMC_Output)
#' 
#' @description Container of a summary of a \code{\link[BASiCS]{BASiCS_Chain-class}} object.  
#' In each slot, first column contains posterior medians, second column contains the lower limits of an high posterior
#' density interval and third column contains the upper limits of high posterior density intervals.
setClass("BASiCS_Summary",
         representation = representation(
           mu = "matrix",
           delta = "matrix",
           phi = "matrix",
           s = "matrix",
           nu = "matrix",
           theta = "matrix"),
         validity = function(object){
           if(sum(!is.finite(object@mu))>0 | sum(!is.finite(object@delta))>0 | 
                sum(!is.finite(object@phi))>0 | sum(!is.finite(object@s))>0 |
                sum(!is.finite(object@nu))>0 | sum(!is.finite(object@theta))>0) stop("Invalid slots") 
           else {TRUE} 
         }
)

#' @name BASiCS_D_Chain-class
#' @aliases BASiCS_D_Chain BASiCS_D_Chain,BASiCS_D_Chain-class
#' 
#' @title The BASiCS_D_Chain class
#' 
#' @description Container of an MCMC sample of the BASiCS model parameters when comparing 2 populations of cells. 
#' For each population, these are independently generated by the function \code{\link[BASiCS]{BASiCS_MCMC}} and 
#' afterwards combined through the function \code{\link[BASiCS]{CombineBASiCS_Chain}} . 
#' 
#' @slot muTest MCMC chain for gene-specific expression levels \eqn{\mu[i]} in test group, defined as true input molecules in case of technical genes 
#' (matrix with \code{q.bio} columns, all elements must be positive numbers)
#' @slot muRef MCMC chain for gene-specific expression levels \eqn{\mu[i]} in reference group, defined as true input molecules in case of technical genes 
#' (matrix with \code{q.bio} columns, all elements must be positive numbers)
#' @slot deltaTest MCMC chain for gene-specific biological cell-to-cell heterogeneity hyper-parameters \eqn{\delta[i]} 
#' in test group, biological genes only (matrix with \code{q} columns, all elements must be positive numbers)
#' @slot deltaRef MCMC chain for gene-specific biological cell-to-cell heterogeneity hyper-parameters \eqn{\delta[i]} in reference group, biological genes only 
#' (matrix with \code{q} columns, all elements must be positive numbers)
#' @slot phi MCMC chain for cell-specific mRNA content normalising constants \eqn{\phi[j]}
#' (matrix with \code{n} columns, all elements must be positive numbers and the sum of its elements must be equal to \code{n})
#' @slot s MCMC chain for cell-specific capture efficiency (or amplification biases if not using UMI based counts) normalising constants \eqn{s[j]}
#' (matrix with \code{n} columns, all elements must be positive numbers)
#' @slot nu MCMC chain for cell-specific random effects \eqn{\nu[j]}
#' (matrix with \code{n} columns, all elements must be positive numbers)
#' @slot thetaTest MCMC chain for technical variability hyper-parameter \eqn{\theta[test]} in the test sample (vector, all elements must be positive)
#' @slot thetaRef MCMC chain for technical variability hyper-parameter \eqn{\theta[ref]} in the reference sample (vector, all elements must be positive)
#' @slot offset Offset value to capture global changes in mRNA content. Default value = 1 (when no offset correction has been performed)
#'   
#' @examples
#' 
#' # See vignette
#' 
#' @author Catalina A. Vallejos \email{catalina.vallejos@@mrc-bsu.cam.ac.uk}
setClass("BASiCS_D_Chain",
         representation = representation(
           muTest = "matrix",
           muRef = "matrix",
           deltaTest = "matrix",
           deltaRef = "matrix",
           phi = "matrix",
           s = "matrix",
           nu = "matrix",
           thetaTest = "matrix",
           thetaRef = "matrix",
           offset = "numeric"),
         validity = function(object){
           errors <- character()
           
           if(length(object@muTest)==0 | length(object@muRef)==0 | 
              length(object@deltaTest)==0 | length(object@deltaRef)==0 | 
              length(object@phi)==0 | length(object@s)==0 | length(object@nu)==0 | 
              length(object@thetaTest)==0 | length(object@thetaRef)==0 | 
              is.null(offset)) {errors <-c(errors,"One or more slots are missing"); stop(errors)}   
           
           N = nrow(object@muTest)
           q.bio = ncol(object@muTest)
           n = ncol(object@phi)
           if(nrow(object@muRef) != N | 
              nrow(object@deltaTest) != N | nrow(object@deltaRef) != N |
              nrow(object@phi) != N | nrow(object@s) != N |
              nrow(object@nu) != N | nrow(object@thetaTest) != N | nrow(object@thetaRef) != N |
              ncol(object@muRef) != q.bio | 
              ncol(object@deltaTest) != q.bio | ncol(object@deltaRef) != q.bio |
              ncol(object@s) != n | ncol(object@nu) != n) {errors <-c(errors,"Slots' dimensions are not compatible")}
           
           if(sum(!is.finite(object@muTest)) + sum(!is.finite(object@muRef)) + 
              sum(!is.finite(object@deltaTest)) + sum(!is.finite(object@deltaRef)) +
              sum(!is.finite(object@phi)) + sum(!is.finite(object@s)) + sum(!is.finite(object@nu)) + 
              sum(!is.finite(object@thetaTest)) + 
              sum(!is.finite(object@thetaRef)) +
              !is.finite(object@offset) > 0) {errors <-c(errors,"One or more of the slots contains NAs or Infinite values")}
           
           if (length(errors) == 0) TRUE else errors
         }
)

#' @name BASiCS_D_Summary-class
#' @aliases BASiCS_D_Summary,BASiCS_D_Summary-class
#' 
#' @title The BASiCS_D_Summary class contains posterior medians and the limits of 
#' the highest posterior density (HPD) interval for all BASiCS model parameters
#' 
#' @slot muTest Posterior medians (first column), lower (second column) and upper (third column) 
#' limits of HPD interval for gene-specific expression levels \eqn{\mu[i]} in test group.
#' @slot muRef Posterior medians (first column), lower (second column) and upper (third column) 
#' limits of HPD interval for gene-specific expression levels \eqn{\mu[i]} in reference group.
#' @slot deltaTest Posterior medians (first column), lower (second column) and upper (third column) 
#' limits of HPD interval for gene-specific biological cell-to-cell heterogeneity hyper-parameters \eqn{\delta[i]} in test group
#' @slot deltaRef Posterior medians (first column), lower (second column) and upper (third column) 
#' limits of HPD interval for gene-specific biological cell-to-cell heterogeneity hyper-parameters \eqn{\delta[i]} in reference group
#' @slot phi Posterior medians (first column), lower (second column) and upper (third column) 
#' limits of HPD interval for cell-specific mRNA content normalising constants \eqn{\phi[j]}
#' @slot s Posterior medians (first column), lower (second column) and upper (third column) 
#' limits of HPD interval for cell-specific capture efficiency (or amplification biases if not using UMI based counts) normalising constants \eqn{s[j]}
#' @slot of HPD interval for Posterior medians (first column), lower (second column) and upper (third column) 
#' limits of HPD interval for cell-specific random effects \eqn{\nu[j]}
#' @slot thetaTest Posterior median (first column), lower (second column) and upper (third column) 
#' limits of HPD interval for technical variability hyper-parameter \eqn{\theta[test]} for test sample (vector, all elements must be positive)
#' @slot thetaRef Posterior median (first column), lower (second column) and upper (third column) 
#' limits of HPD interval for technical variability hyper-parameter \eqn{\theta[ref]} for reference sample (vector, all elements must be positive)
#' @slot offset Offset value to capture global changes in mRNA content. Default value = 1 (when no offset correction has been performed)
#' @slot probHPD Posterior probability contained in the HPD interval. Must be a value between 0 and 1. 
#'  
#' @examples
#' 
#' # see help(BASiCS_D_MCMC)
#' 
#' @description Container of a summary of a \code{\link[BASiCS]{BASiCS_D_Chain-class}} object.  
setClass("BASiCS_D_Summary",
         representation = representation(
           muTest = "matrix",
           muRef = "matrix",
           deltaTest = "matrix",
           deltaRef = "matrix",
           phi = "matrix",
           s = "matrix",
           nu = "matrix",
           thetaTest = "matrix",
           thetaRef = "matrix",
           offset = "numeric",
           probHPD = "numeric"),
         validity = function(object){
           if(sum(!is.finite(object@muTest))>0 | sum(!is.finite(object@muRef))>0 | 
                sum(!is.finite(object@deltaTest))>0 | sum(!is.finite(object@deltaRef))>0 | 
                sum(!is.finite(object@phi))>0 | sum(!is.finite(object@s))>0 |
                sum(!is.finite(object@nu))>0 | 
                sum(!is.finite(object@thetaTest))>0 | sum(!is.finite(object@thetaRef))>0 |
                is.null(object@offset)  | !is.finite(object@offset) | is.na(object@offset) |
                is.null(object@probHPD) | is.na(object@probHPD) | !is.finite(object@probHPD) | 
                object@probHPD < 0 | object@probHPD > 1) stop("Invalid slots") 
           else {TRUE} 
         }
)
