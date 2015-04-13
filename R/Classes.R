#' @name BASiCS_Data-class
#' @aliases BASiCS_Data,BASiCS_Data-class
#' 
#' @title The BASiCS_Data class 
#' 
#' @description Container of expression counts from single-cell 
#' sequencing experiments in the format required by BASiCS (see Vallejos et al 2015).
#' 
#' @slot Counts Matrix of dimensions \code{q} times \code{n} whose elements corresponds to the simulated expression counts. 
#' First \code{q.bio} rows correspond to biological genes. Last \code{q-q.bio} rows correspond to technical spike-in genes. 
#' @slot Tech Logical vector of length \code{q}. If \code{Tech = F} the gene is biological; otherwise the gene is spike-in.
#' @slot SpikeInput Vector of length \code{q-q.bio} whose elements indicate the simulated input concentrations for the spike-in genes. 
#' 
#' @examples
#'
#' CountsAux = matrix(rpois(10*5, 1), ncol = 5)
#' TechAux = c(rep(FALSE,7),rep(TRUE,3))
#' SpikeInputAux = rgamma(3,1,1)
#' Data = new('BASiCS_Data', Counts = CountsAux, Tech = TechAux, SpikeInput = SpikeInputAux)
#' head(counts(Data))
#' dim(counts(Data, type="biological"))
#' dim(counts(Data, type="technical"))
#' displayTechIndicator(Data)
#' displaySpikeInput(Data)
#' 
#' # For usage, see documentation of functions: BASiCS_Sim, BASiCS_MCMC_Start, BASiCS_MCMC.
#' 
#' @author Catalina A. Vallejos \email{catalina.vallejos@@mrc-bsu.cam.ac.uk}
#' 
#' @seealso \code{\link[BASiCS]{makeExampleBASiCS_Data}}, \code{\link[BASiCS]{BASiCS_Data-methods}}, \code{\link[BASiCS]{BASiCS_MCMC_Start}}, \code{\link[BASiCS]{BASiCS_MCMC}}.
#' 
#' @references Vallejos, Marioni and Richardson (2015). Bayesian Analysis of Single-Cell Sequencing data.
setClass("BASiCS_Data",
         representation = representation(
           Counts = "matrix",
           Tech = "vector",
           SpikeInput = "vector"),
         validity = function(object){
           errors <- character()
           
           if(!(is.numeric(object@Counts) & all(object@Counts>=0) & 
                  sum(!is.finite(object@Counts))==0 )) errors <- c(errors, "Invalid value for Counts")
           
           if(!(is.logical(object@Tech))) errors <- c(errors, "Invalid value for Tech")
           
           if(!(is.numeric(object@SpikeInput) & all(object@SpikeInput>0) & 
                  sum(!is.finite(object@SpikeInput))==0 )) errors <- c(errors, "Invalid value for SpikeInput.")
           
           q = nrow(object@Counts)
           q.bio = q - length(object@SpikeInput)
           n = ncol(object@Counts)
           
           if(!( length(object@Tech) == q & sum(!object@Tech) == q.bio )) errors <- c(errors, "Argument's dimensions are not compatible.")
           
           if(!( sum(object@Tech[1:q.bio]) == 0 & sum(object@Tech[(q.bio+1):q])==q-q.bio )) errors <- c(errors, "Expression counts are not in the right format (spike-in genes must be at the bottom of the matrix).")
                        
           if (length(errors) == 0) TRUE else errors
         }
)

#' @name BASiCS_Chain-class
#' @aliases BASiCS_Chain,BASiCS_Chain-class
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
#' @slot theta MCMC chain for technical variability hyper-parameter \eqn{\theta} (vector, all elements must be positive)
#'   
#' @examples
#' 
#' # A BASiCS_Chain object created by the BASiCS_MCMC function.
#' Data = makeExampleBASiCS_Data()
#' Result <- BASiCS_MCMC(Data, N = 100, Thin = 2, Burn = 2)
#' Result
#' 
#' @author Catalina A. Vallejos \email{catalina.vallejos@@mrc-bsu.cam.ac.uk}
#' 
#' @references Vallejos, Marioni and Richardson (2015). Bayesian Analysis of Single-Cell Sequencing data.
setClass("BASiCS_Chain",
         representation = representation(
           mu = "matrix",
           delta = "matrix",
           phi = "matrix",
           s = "matrix",
           nu = "matrix",
           theta = "vector"),
         validity = function(object){
           errors <- character()
           
           if(length(object@mu)==0 | length(object@delta)==0 | 
                length(object@phi)==0 | length(object@s)==0 | 
                length(object@nu)==0 | length(object@theta)==0) {errors <-c(errors,"One or more slots are missing"); stop(errors)}   
           
           N = nrow(object@mu)
           n = ncol(object@phi)
           if(nrow(object@delta) != N |
              nrow(object@phi) != N | nrow(object@s) != N |
              nrow(object@nu) != N | length(object@theta) != N |
              ncol(object@mu) <= ncol(object@delta) |
              ncol(object@s) != n | ncol(object@nu) != n) {errors <-c(errors,"Slots' dimensions are not compatible")}
           
           if(sum(!is.finite(object@mu)) + sum(!is.finite(object@delta)) +
              sum(!is.finite(object@phi)) + sum(!is.finite(object@s)) +
              sum(!is.finite(object@nu)) + sum(!is.finite(object@theta))) {errors <-c(errors,"One or more of the slots contains NAs or Infinite values")}
                    
           if (length(errors) == 0) TRUE else errors
         }
         )


# Redefine this class. For each slot, it will contain posterior medians and HPD95% limits for the corresponding parameter. 
#' @name BASiCS_Summary-class
#' @aliases BASiCS_Summary,BASiCS_Summary-class
#' 
#' @title The BASiCS_Summary class
#' 
#' @description This class contains valid values for BASiCS model parameters (see Vallejos et al, 2015)
setClass("BASiCS_Summary",
         representation = representation(
           mu = "matrix",
           delta = "matrix",
           phi = "matrix",
           s = "matrix",
           nu = "matrix",
           theta = "matrix"),
         validity = function(object){
           
#           if(!is.list(object@PostMedian) | !is.list(object@HPD_UpLimit) | !is.list(object@HPD_LowLimit)) stop("Invalid slots")  
           TRUE 
         }
)