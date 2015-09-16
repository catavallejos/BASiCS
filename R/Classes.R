#' @name BASiCS_Data-class
#' @aliases BASiCS_Data,BASiCS_Data-class
#' 
#' @title The BASiCS_Data class 
#' 
#' @description Container of expression counts from single-cell 
#' sequencing experiments in the format required by BASiCS (see Vallejos et al 2015).
#' 
#' @slot Counts Matrix of dimensions \code{q} times \code{n} whose elements corresponds to the raw expression counts. 
#' First \code{q.bio} rows correspond to biological genes. Last \code{q-q.bio} rows correspond to technical spike-in genes. 
#' @slot Tech Logical vector of length \code{q}. If \code{Tech = F} the gene is biological; otherwise the gene is spike-in.
#' @slot SpikeInput Vector of length \code{q-q.bio} whose elements indicate the input number of molecules for the spike-in genes (amount per cell). 
#' @slot GeneNames Vector of length \code{q} containing gene names. Default value: \code{GeneNames = paste("Gene", 1:q)}, 
#' with numbering order as in the input dataset. 
#' @slot BatchInfo Vector of lenght \code{n} containing batch indicators. 
#' 
#' @examples
#'
#' Counts = matrix(rpois(10*5, 1), ncol = 5)
#' Tech = c(rep(FALSE,7),rep(TRUE,3))
#' SpikeInput = rgamma(3,1,1)
#' Data = newBASiCS_Data(Counts, Tech, SpikeInput)
#' 
#' head(counts(Data))
#' dim(counts(Data, type="biological"))
#' dim(counts(Data, type="technical"))
#' displayTechIndicator(Data)
#' displaySpikeInput(Data)
#' displayBatchInfo(Data)
#' 
#' @author Catalina A. Vallejos \email{catalina.vallejos@@mrc-bsu.cam.ac.uk}
#' 
#' @seealso \code{\link[BASiCS]{BASiCS_Data-methods}}, \code{\link[BASiCS]{BASiCS_MCMC}}.
#' 
#' @references Vallejos, Marioni and Richardson (2015). Bayesian Analysis of Single-Cell Sequencing data.
setClass("BASiCS_Data",
         representation = representation(
           Counts = "matrix",
           Tech = "vector",
           SpikeInput = "vector", 
           GeneNames = "vector", 
           BatchInfo = "vector"),
         validity = function(object){
           errors <- character()
           
           if(!(is.numeric(object@Counts) & all(object@Counts>=0) & 
                  sum(!is.finite(object@Counts))==0 )) errors <- c(errors, "Invalid value for Counts")
           
           if(sum(object@Counts %% 1) > 0) errors <- c(errors, "Invalid value for Counts (entries must be positive integers)")          
           
           if(!(is.logical(object@Tech))) errors <- c(errors, "Invalid value for Tech")
           
           if(!(is.numeric(object@SpikeInput) & all(object@SpikeInput>0) & 
                  sum(!is.finite(object@SpikeInput))==0 )) errors <- c(errors, "Invalid value for SpikeInput.")
           
           q = nrow(object@Counts)
           q.bio = q - length(object@SpikeInput)
           n = ncol(object@Counts)
          
           if(!( length(object@Tech) == q & sum(!object@Tech) == q.bio )) 
             errors <- c(errors, "Argument's dimensions are not compatible.")
           
           if(length(object@GeneNames) != q)
             errors <- c(errors, "Incorrect length of the vector stored in the GeneNames slot.")
           
           if(!( sum(object@Tech[1:q.bio]) == 0 & sum(object@Tech[(q.bio+1):q])==q-q.bio )) 
             errors <- c(errors, "Expression counts are not in the right format (spike-in genes must be at the bottom of the matrix).")
           
           if(sum(apply(object@Counts[ object@Tech,],2,sum) == 0) > 0) 
             errors <- c(errors, "Some cells have zero reads mapping back to the spike-in genes. Please remove them before creating the BASiCS_Data object.")
           
           if(sum(apply(object@Counts[!object@Tech,],2,sum) == 0) > 0) 
             errors <- c(errors, "Some cells have zero reads mapping back to the intrinsic genes. Please remove them before creating the BASiCS_Data object.")
           
           if(sum(apply(object@Counts,1,sum) == 0) > 0) 
             errors <- c(errors, "Some genes have zero counts across all cells. Please remove them before creating the BASiCS_Data object.")
           
           if(sum(apply(object@Counts,1,sum) > 0) == 1) 
             errors <- c(errors, "Some genes have non-zero counts only in 1 cell. Please remove them before creating the BASiCS_Data object.")
           
           if(length(object@BatchInfo) != n) 
             errors <- c(errors, "BatchInfo slot is not compatible with the number of cells contained in Counts slot.")
           
           if (length(errors) == 0) TRUE else errors
         }
)

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
#' @references Vallejos, Marioni and Richardson (2015). Bayesian Analysis of Single-Cell Sequencing data.
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