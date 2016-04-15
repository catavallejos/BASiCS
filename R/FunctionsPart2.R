#' @title Creates a BASiCS_D_Data object from expression counts matrices and experimental information about spike-in genes
#' 
#' @description \code{newBASiCS_D_Data} creates a \code{\link[BASiCS]{BASiCS_D_Data-class}} object from 
#' expression counts matrices and experimental information about spike-in genes.
#' 
#' @param CountsTest Matrix of dimensions \code{q} times \code{n_test} whose elements corresponds to the raw expression counts in the test sample. 
#' First \code{q.bio} rows correspond to biological genes. Last \code{q-q.bio} rows correspond to technical spike-in genes. 
#' @param CountsRef Matrix of dimensions \code{q} times \code{n_ref} whose elements corresponds to the raw expression counts in the reference sample. 
#' First \code{q.bio} rows correspond to biological genes. Last \code{q-q.bio} rows correspond to technical spike-in genes.
#' @param Tech Logical vector of length \code{q}. If \code{Tech = F} the gene is biological; otherwise the gene is spike-in.
#' @param SpikeInputTest Vector of length \code{q-q.bio} whose elements indicate the input number of molecules for the spike-in genes in the test sample (amount per cell). 
#' @param SpikeInputRef Vector of length \code{q-q.bio} whose elements indicate the input number of molecules for the spike-in genes in the reference sample (amount per cell). 
#' @param BatchInfoTest Vector of length \code{n_test} indicating experimental batches for test group. 
#' @param BatchInfoRef Vector of length \code{n_ref} indicating experimental batches for test group. 
#' @param GeneNames Vector of gene ids. 
#' 
#' @return An object of class \code{\link[BASiCS]{BASiCS_D_Data-class}}.
#' 
#' @examples
#' 
#' # Rather than using this function, we recomend to use 'CombineBASiCS_Data'
#' 
#' @seealso \code{\link[BASiCS]{BASiCS_D_Data-class}}
#' 
#' @author Catalina A. Vallejos \email{catalina.vallejos@@mrc-bsu.cam.ac.uk}
newBASiCS_D_Data <- function(CountsTest, 
                              CountsRef, 
                              Tech, 
                              SpikeInputTest, 
                              SpikeInputRef,
                              BatchInfoTest = NULL,
                              BatchInfoRef = NULL,
                              GeneNames)
{
  if (is.null(BatchInfoTest)) { BatchInfoTest = rep(1, times = ncol(CountsTest))}
  if (is.null(BatchInfoRef)) { BatchInfoRef = rep(1, times = ncol(CountsRef))}
  
  Data <- new("BASiCS_D_Data", CountsTest = CountsTest, CountsRef = CountsRef, Tech = Tech, 
              SpikeInputTest = SpikeInputTest, SpikeInputRef = SpikeInputRef, 
              BatchInfoTest = BatchInfoTest, BatchInfoRef = BatchInfoRef, GeneNames = GeneNames)
#  show(Data)
  cat('\n')
  cat('NOTICE: BASiCS requires a pre-filtered dataset \n')
  cat('    - You must remove poor quality cells before creating the BASiCS data object \n')
  cat('    - We recommend to pre-filter very lowly expressed transcripts before creating the object. \n')
  cat('      Inclusion criteria may vary for each data. For example, remove transcripts \n')
  cat('          - with very low total counts across of all cells \n')
  cat('          - that are only expressed in few cells \n')
  cat('            (by default genes expressed in only 1 cell are not accepted) \n')
  cat('          - with very low total counts across the cells where the transcript is expressed \n')
  cat('\n')
  cat(' BASiCS_Filter can be used for this purpose \n \n')
  
  return(Data)  
}

#' @title Creates a BASiCS_D_Data object based on 2 independent `BASiCS_Data` objects
#' 
#' @description Creates a \code{\link[BASiCS]{BASiCS_D_Data-class}} object from 
#' based on 2 independent  \code{\link[BASiCS]{BASiCS_Data-class}} objects
#' 
#' @param DataTest A \code{\link[BASiCS]{BASiCS_Data-class}} object (test group)
#' @param DataRef A \code{\link[BASiCS]{BASiCS_Data-class}} object (reference group)
#' 
#' @return An object of class \code{\link[BASiCS]{BASiCS_D_Data-class}}.
#' 
#' @examples
#' 
#' DataTest = makeExampleBASiCS_Data()
#' DataRef = makeExampleBASiCS_Data()
#' Data = CombineBASiCS_Data(DataTest, DataRef)
#' 
#' @seealso \code{\link[BASiCS]{BASiCS_D_Data-class}}
#' 
#' @author Catalina A. Vallejos \email{catalina.vallejos@@mrc-bsu.cam.ac.uk}
CombineBASiCS_Data <- function(DataTest, DataRef)
{
  if(sum(DataTest@GeneNames != DataRef@GeneNames) > 0) stop("DataTest and DataRef do not contain the same genes and/or they are in a different order")
  Data <- newBASiCS_D_Data(CountsTest = DataTest@Counts,
                            CountsRef = DataRef@Counts, 
                            Tech = DataTest@Tech,
                            SpikeInputTest = DataTest@SpikeInput,
                            SpikeInputRef = DataRef@SpikeInput,
                            BatchInfoTest = DataTest@BatchInfo, 
                            BatchInfoRef = DataRef@BatchInfo,
                            GeneNames = DataTest@GeneNames)
  show(Data)
  
  Data
}

#' @name makeExampleBASiCS_D_Data
#' 
#' @title Create a simple example of a BASiCS_D_Data object with random data
#' 
#' @description A simple \code{\link[BASiCS]{BASiCS_D_Data-class}} object is generated by simulating a dataset from the
#' model in BASiCS. This is used for the examples in the package and the vignette. 
#' 
#' @return An object of class \code{\link[BASiCS]{BASiCS_D_Data-class}}, simulated from the model implemented in BASiCS (for 2 populations of cells). 
#' It contains 120 genes (100 biological and 20 spike-in) and 20 cells (10 in each sample). 
#' 
#' @examples 
#' Data = makeExampleBASiCS_D_Data()
#' is(Data, "BASiCS_D_Data")
#' 
#' @author Catalina A. Vallejos \email{catalina.vallejos@@mrc-bsu.cam.ac.uk}
makeExampleBASiCS_D_Data <- function()
{
  mu =  c( 8.36,  10.65,   4.88,   6.29,  12.93,   3.89,   6.34,  12.45,   8.08,   7.31,  
           15.56,  15.91,  12.24,  15.96,  19.14,   4.20,   6.75,  27.74,   8.88,  21.21,  
           19.89,   7.14,  11.09,   7.19,  20.64,  73.90,   9.05,   6.13,  16.40,   6.08,  
           17.89,   6.98,  10.25,  14.05,   8.14,   5.67,   6.95,  11.16,  11.83,   7.56, 
           159.05,  16.41,   4.58,  15.46,  10.96,  25.05,  21.72,  30.19,  83.92,  57.87,  
           18.54,   8.50,   4.05,   5.37,   4.63,   4.08,   3.75,   5.02,  27.74,  10.28,   
           3.91,  13.10,   8.23,   3.64,  80.77,  20.36,   3.47,  20.12,  13.29,   7.92,  
           25.51, 173.01,  27.71,   4.89,  33.10,   3.42,   6.19,   4.29,   5.19,   8.36,  
           10.27,   6.86,   5.72,  37.25,   3.82,  23.97,   5.80,  14.30,  29.07,   5.30,   
           7.47,   8.44,   4.24,  16.15,  23.39, 120.22,   8.92,  97.15,   9.75,  10.07,
           4.74, 176.19,   8.57,   6.93,  29.28,   9.75,   9.42,  14.41, 152.29,  80.91,  
           26.54,  34.98,   7.07,  41.33,  13.51,   4.54,   8.97,  21.72,   4.78,   6.07,   
           6.36,  10.44, 212.32,  15.26,  55.39,  14.86,   5.31,   7.70,  12.54,  11.91,
           6.02, 107.51,  17.09,   9.50,  10.73,  36.97,   5.05,   7.11,   5.69,   0.05, 
           1010.72,   7.90,  31.59,  63.17,   1.97, 252.68,  31.59,  31.59,  31.59,  63.17, 
           4042.89,4042.89,   3.95, 126.34, 252.68,   7.90,  15.79, 126.34,   3.95,  15.79,
           126.34,   3.95,  15.79,   0.99,1010.72,  15.79,2021.45,4042.89,  31.59,  63.17)
  tau = c(rep(2, times = 5), rep(1.5, times = 5), rep(0, times = 130))
  delta = c(1.29, 0.88, 1.51, 1.49, 0.27, 0.53, 1.31, 0.81, 0.72, 0.70, 
            0.96, 0.58, 1.15, 0.82, 0.25, 5.32, 1.13, 0.31, 0.66, 0.27, 
            0.76, 1.39, 1.18, 1.57, 0.55, 0.17, 1.40, 1.47, 0.57, 2.55, 
            0.62, 0.77, 1.47, 0.91, 1.53, 2.89, 1.43, 0.77, 1.37, 0.57, 
            0.15, 0.33, 3.99, 0.47, 0.87, 0.86, 0.54, 0.40, 0.85, 0.26,
            0.97, 1.25, 2.20, 2.19, 1.26, 1.89, 1.70, 1.89, 0.69, 1.63,
            2.83, 0.29, 1.21, 2.06, 0.20, 0.34, 0.71, 0.61, 0.71, 1.20, 
            0.88, 0.17, 0.25, 1.48, 0.31, 2.49, 2.75, 1.43, 2.65, 1.97,
            0.84, 0.81, 2.75, 0.53, 2.23, 0.45, 1.87, 0.74, 0.53, 0.58, 
            0.94, 0.72, 2.61, 1.56, 0.37, 0.07, 0.90, 0.12, 0.76, 1.45, 
            1.04, 0.16, 2.16, 1.90, 0.28, 1.01, 1.23, 0.46, 0.13, 0.20, 
            0.16, 0.44, 0.84, 0.08, 0.69, 2.62, 1.00, 0.26, 1.50, 2.74, 
            0.92, 0.71, 0.11, 0.66, 0.34, 0.22, 1.32, 0.99, 0.98, 0.79, 
            1.68, 0.08, 0.39, 0.84, 1.10, 0.24, 0.71, 1.00, 2.06, 0.05) 
  omega = c(2, 1.5, 1, 0, 0, 2, 1.5, 1, 0, 0, 2, 1.5, 1, 0, 0, rep(0, times = 125)) 
  phi = c(1.09, 1.16, 1.19, 1.14, 0.87, 1.10, 0.48, 1.06, 0.94, 0.97, 1.09, 1.16, 1.19, 1.14, 0.87, 1.10, 0.48, 1.06, 0.94, 0.97,
          1.09, 1.16, 1.19, 1.14, 0.87, 1.10, 0.48, 1.06, 0.94, 0.97, 1.09, 1.16, 1.19, 1.14, 0.87, 1.10, 0.48, 1.06, 0.94, 0.97,
          1.09, 1.16, 1.19, 1.14, 0.87, 1.10, 0.48, 1.06, 0.94, 0.97, 1.09, 1.16, 1.19, 1.14, 0.87, 1.10, 0.48, 1.06, 0.94, 0.97,
          1.09, 1.16, 1.19, 1.14, 0.87, 1.10, 0.48, 1.06, 0.94, 0.97, 1.09, 1.16, 1.19, 1.14, 0.87, 1.10, 0.48, 1.06, 0.94, 0.97,
          1.09, 1.16, 1.19, 1.14, 0.87, 1.10, 0.48, 1.06, 0.94, 0.97, 1.09, 1.16, 1.19, 1.14, 0.87, 1.10, 0.48, 1.06, 0.94, 0.97)
  s = c(0.38, 0.40, 0.38, 0.39, 0.34, 0.39, 0.31, 0.39, 0.40, 0.37, 0.38, 0.40, 0.38, 0.39, 0.34, 0.39, 0.31, 0.39, 0.40, 0.37,
        0.38, 0.40, 0.38, 0.39, 0.34, 0.39, 0.31, 0.39, 0.40, 0.37, 0.38, 0.40, 0.38, 0.39, 0.34, 0.39, 0.31, 0.39, 0.40, 0.37,
        0.38, 0.40, 0.38, 0.39, 0.34, 0.39, 0.31, 0.39, 0.40, 0.37, 0.38, 0.40, 0.38, 0.39, 0.34, 0.39, 0.31, 0.39, 0.40, 0.37,
        0.38, 0.40, 0.38, 0.39, 0.34, 0.39, 0.31, 0.39, 0.40, 0.37, 0.38, 0.40, 0.38, 0.39, 0.34, 0.39, 0.31, 0.39, 0.40, 0.37,
        0.38, 0.40, 0.38, 0.39, 0.34, 0.39, 0.31, 0.39, 0.40, 0.37, 0.38, 0.40, 0.38, 0.39, 0.34, 0.39, 0.31, 0.39, 0.40, 0.37)
  theta_test = 0.5
  theta_ref = 0.5
  
  n_test = 50; n_ref = 50
  
  q = length(mu); q.bio = length(delta); n = length(phi)
  
  # Matrix where simulated counts will be stored
  Counts.sim<-matrix(0,ncol=n,nrow=q)
  # Matrix where gene-cell specific simulated random effects will be stored 
  rho<-matrix(1,ncol=n,nrow=q.bio)
  # Simulated cell-specific random effects
  test.ind=c(rep(T,n_test),rep(F,n-n_test))
  theta=test.ind*theta_test+(1-test.ind)*theta_ref
  if(all(theta>0)) {set.seed(1000); nu<-rgamma(n,shape=1/theta,rate=1/(s*theta))}
  else {nu<-s}
  # Simulated counts data
  for(i in 1:q)
  {
    # Biological genes
    if(i<=q.bio)
    {
      if(delta[i]>0)
      {
        set.seed(i); 
        rho[i,]<-rgamma(n, shape = exp(-(test.ind*omega[i]+(1-test.ind)*0))/delta[i], rate = exp(-(test.ind*omega[i]+(1-test.ind)*0))/delta[i])
      }
      set.seed(i+50000); Counts.sim[i,]<-rpois(n,lambda=phi*nu*mu[i]*exp(test.ind*tau[i]+(1-test.ind)*0)*rho[i,])
    }
    # Technical genes
    else {set.seed(i+20000); Counts.sim[i,]<-rpois(n,lambda=nu*mu[i])}
  }
  
  Data = newBASiCS_D_Data(CountsTest = Counts.sim[,1:n_test], 
                           CountsRef = Counts.sim[,(n_test+1):n], 
                           Tech = ifelse(1:q > q.bio, T, F), 
                           SpikeInputTest = mu[(q.bio+1):q], 
                           SpikeInputRef = mu[(q.bio+1):q],
                           GeneNames = paste0("Gene",1:q))    
  return(Data)  
}

HiddenOffSetCorrection <- function(MCMC_Output) # BASiCS_D_Chain object
{ 
  median(rowSums(MCMC_Output@muTest)/rowSums(MCMC_Output@muRef)) 
}

#' @title Creates a BASiCS_D_Chain object from pre-computed MCMC chains
#' 
#' @description \code{BASiCS_D_Chain} creates a \code{\link[BASiCS]{BASiCS_D_Chain-class}} object from pre-computed MCMC chains.
#' 
#' @param muTest MCMC chain for gene-specific expression levels \eqn{\mu[i]} (test group), defined as true input molecules in case of technical genes 
#' (matrix with \code{q} columns, technical genes located at the end of the matrix, all elements must be positive numbers)
#' @param muRef MCMC chain for gene-specific log-fold changes \eqn{\mu[i]} (reference group), defined as (log) difference in input molecules in case of technical genes 
#' (matrix with \code{q} columns, technical genes located at the end of the matrix)
#' @param deltaTest MCMC chain for gene-specific biological cell-to-cell heterogeneity hyper-parameters \eqn{\delta[i]} (test group), biological genes only 
#' (matrix with \code{q.bio} columns, all elements must be positive numbers)
#' @param deltaRef MCMC chain for gene-specific log-fold change in biological cell-to-cell heterogeneity hyper-parameters \eqn{\delta[i]} (reference group), biological genes only 
#' (matrix with \code{q.bio} columns, all elements must be real numbers)
#' @param phi MCMC chain for cell-specific mRNA content normalising constants \eqn{\phi[j]}
#' (matrix with \code{n} columns, all elements must be positive numbers and the sum of its elements must be equal to \code{n})
#' @param s MCMC chain for cell-specific capture efficiency (or amplification biases if not using UMI based counts) normalising constants \eqn{s[j]}
#' (matrix with \code{n} columns, all elements must be positive numbers)
#' @param nu MCMC chain for cell-specific random effects \eqn{\nu[j]}
#' (matrix with \code{n} columns, all elements must be positive numbers)
#' @param thetaTest MCMC chain for technical variability hyper-parameter \eqn{\theta_test} in the test sample (vector, all elements must be positive)
#' @param thetaRef MCMC chain for technical variability hyper-parameter \eqn{\theta_ref} in the reference sample (vector, all elements must be positive)
#' @param offset Offset value to be corrected (default = \code{NULL} to be internally calculated)
#' #@param offsetCorrect \code{TRUE}/\code{FALSE} value to indicate if offset correction is required (default value = \code{TRUE})
#' 
#' @return An object of class \code{\link[BASiCS]{BASiCS_D_Chain-class}}. 
#' 
#' @examples
#' 
#' # Data = makeExampleBASiCS_D_Data()
#' # MCMC_Output <- BASiCS_D_MCMC(Data, N = 50, Thin = 5, Burn = 5, 
#' #                StoreChains = TRUE, StoreDir = getwd(), RunName = "Example")
#' 
#' # ChainMuTest = as.matrix(read.table("chain_muTest_Example.txt"))
#' # ChainMuRef = as.matrix(read.table("chain_muRef_Example.txt"))
#' # ChainDeltaTest = as.matrix(read.table("chain_deltaTest_Example.txt"))
#' # ChainDeltaRef = as.matrix(read.table("chain_deltaRef_Example.txt"))
#' # ChainPhi = as.matrix(read.table("chain_phi_Example.txt"))
#' # ChainS = as.matrix(read.table("chain_s_Example.txt"))
#' # ChainNu = as.matrix(read.table("chain_nu_Example.txt"))
#' # ChainThetaTest = read.table("chain_thetaTest_Example.txt")[,1]
#' # ChainThetaRef = read.table("chain_thetaRef_Example.txt")[,1]
#' 
#' # MCMC_Output_Load <- newBASiCS_D_Chain(muTest = ChainMuTest, 
#' #                                        muRef = ChainMuRef, 
#' #                                        deltaTest = ChainDeltaTest, 
#' #                                        deltaRef = ChainDeltaRef, 
#' #                                        phi = ChainPhi, 
#' #                                        s = ChainS, 
#' #                                        nu = ChainNu, 
#' #                                        thetaTest = ChainThetaTest, 
#' #                                        thetaRef = ChainThetaRef)
#' 
#' @seealso \code{\link[BASiCS]{BASiCS_D_Chain-class}}
#' 
#' @author Catalina A. Vallejos \email{catalina.vallejos@@mrc-bsu.cam.ac.uk}
newBASiCS_D_Chain <- function(muTest, muRef, deltaTest, deltaRef, 
                               phi, s, nu, thetaTest, thetaRef, offset = 1)
{
  Chain <- new("BASiCS_D_Chain", muTest = muTest, muRef = muRef, 
               deltaTest = deltaTest, deltaRef = deltaRef, 
               phi = phi, s = s, nu = nu, 
               thetaTest = thetaTest, thetaRef = thetaRef, offset = offset)
  show(Chain)
  return(Chain)  
}

#' @title Creates a BASiCS_D_Chain object based on 2 independent `BASiCS_Chain` objects
#' 
#' @description Creates a \code{\link[BASiCS]{BASiCS_D_Chain-class}} object from 
#' based on 2 independent  \code{\link[BASiCS]{BASiCS_Chain-class}} objects
#' 
#' @param ChainTest A \code{\link[BASiCS]{BASiCS_Chain-class}} object (test group)
#' @param ChainRef A \code{\link[BASiCS]{BASiCS_Chain-class}} object (reference group)
#' 
#' @return An object of class \code{\link[BASiCS]{BASiCS_D_Chain-class}}.
#' 
#' @examples
#' 
#' # See vignette
#' 
#' @seealso \code{\link[BASiCS]{BASiCS_D_Data-class}}
#' 
#' @author Catalina A. Vallejos \email{catalina.vallejos@@mrc-bsu.cam.ac.uk}
CombineBASiCS_Chain <- function(ChainTest, ChainRef)
{
  Chain <- newBASiCS_D_Chain(muTest = ChainTest@mu,
                              muRef = ChainRef@mu,
                              deltaTest = ChainTest@delta,
                              deltaRef = ChainRef@delta,
                              phi = cbind(ChainTest@phi, ChainRef@phi),
                              s = cbind(ChainTest@s, ChainRef@s),
                              nu = cbind(ChainTest@nu, ChainRef@nu),
                              thetaTest = ChainTest@theta,
                              thetaRef = ChainRef@theta)
  Chain
}



########################################################################
## Required functions for differential expression test: Bayes Factors
########################################################################

HiddenTailProbUpDV<-function(chain,threshold){return(sum(chain>threshold)/length(chain))}
HiddenTailProbLowDV<-function(chain,threshold){return(sum(chain<threshold)/length(chain))}

HiddenProbDE<-function(
  chain, # MCMC chain for log-fold change parameter (in absolute value)
  tol) # Minimum tolerance difference for log-fold change
{  
  return(apply(chain, 2, HiddenTailProbUpDV, threshold = tol))
}

HiddenEFDRDV<-function(
  EviThreshold, # Evidence threshold (it must be contained in (0,1))
  Prob) # Probability of changes in expression (for a given tolerance)
{  
  return(sum((1-Prob)*I(Prob > EviThreshold))/sum(I(Prob > EviThreshold)))
}

HiddenEFNRDV<-function(
  EviThreshold, # Evidence threshold (it must be contained in (0,1))
  Prob) # Probability of changes in expression (for a given tolerance) 
{
  return(sum(Prob*I(Prob <= EviThreshold))/sum(I(Prob <= EviThreshold)))
}

#' @name BASiCS_D_TestDE
#' 
#' @title Detection of genes with changes in expression
#' 
#' @description Function to assess changes in expression (mean and over-dispersion)
#' 
#' @param Data an object of class \code{\link[BASiCS]{BASiCS_D_Data-class}}
#' @param object an object of class \code{\link[BASiCS]{BASiCS_D_Chain-class}}
#' @param GeneNames Vector containing gene names to be used in results table (argument to be removed as 'GeneNames' will be an slot of `BASiCS_D_Data` object)
#' @param EpsilonM Minimum fold change tolerance threshold for detecting changes in overall expression (must be a positive real number)
#' @param EpsilonD Minimum fold change tolerance threshold for detecting changes in cell-to-cell biological over dispersion (must be a positive real number)
#' @param EviThresholdM Optional parameter. Evidence threshold for detecting changes in overall expression (must be a positive value, between 0 and 1)
#' @param EviThresholdD Optional parameter. Evidence threshold for detecting changes in cell-to-cell biological over dispersion (must be a positive value, between 0 and 1)
#' @param OrderVariable Ordering variable for output. Must take values in \code{c("GeneIndex", "GeneNames", "ProbDiffExp", "ProbDiffOverDisp")}.
#' @param GroupLabelRef Label assigned to reference group. Default: \code{GroupLabelRef = "Ref"}
#' @param GroupLabelTest Label assigned to reference group. Default: \code{GroupLabelRef = "Test"}
#' @param Plot If \code{Plot = T}, error rates control rates and volcano plots are generated.  
#' @param OffSet Optional argument to remove a fix offset effect (if not previously removed from the MCMC chains). This argument will be removed shorly, once offset removal is built as an internal step. 
#' @param EFDR_M Target for expected false discovery rate related to the comparison of means (default = 0.05)
#' @param EFDR_D Target for expected false discovery rate related to the comparison of dispersions (default = 0.05)
#' @param GenesSelect Optional argument to provide a user-defined list of genes to be considered for the comparison (default = NULL). When used, this argument must be a vector of 'TRUE' (include gene) / 'FALSE' (exclude gene) indicator, with the same length as the number of intrinsic genes and following the same order as how genes are displayed in the table of counts.  This argument is necessary in order to have a meaningful EFDR calibration when the user decides to exclude some genes from the comparison. 
#' @param ... Graphical parameters (see \code{\link[graphics]{par}}).
#' 
#' @return \code{BASiCS_D_TestDE} returns a list of 4 elements:
#' \describe{
#' \item{\code{Table}}{A \code{\link[base]{data.frame}} containing the results of the differential expression test}
#'    \describe{
#'    \item{\code{Mu}}{Vector of length \code{q.bio}. For each biological gene, posterior median of gene-specific expression levels \eqn{\mu[i]}}
#'    \item{\code{Delta}}{Vector of length \code{q.bio}. For each biological gene, posterior median of gene-specific biological cell-to-cell heterogeneity hyper-parameter \eqn{\delta[i]}}
#'    \item{\code{Sigma}}{Vector of length \code{q.bio}. For each biological gene, proportion of the total variability that is due to a cell-to-cell biological heterogeneity component. }
#'    \item{\code{Prob}}{Vector of length \code{q.bio}. For each biological gene, probability of being highly variable according to the given thresholds.}
#'    \item{\code{HVG}}{Vector of length \code{q.bio}. For each biological gene, indicator of being detected as highly variable according to the given thresholds. }
#'    }
#' \item{\code{EviThreshold}}{Evidence thresholds.}
#' \item{\code{EFDR}}{Expected false discovery rate for the given thresholds.}
#' \item{\code{EFNR}}{Expected false negative rate for the given thresholds.}
#' }
#'  
#' @examples
#' 
#' # See vignette
#' 
#' @details See vignette
#' 
#' 
#' @seealso \code{\link[BASiCS]{BASiCS_D_Chain-class}} 
#' 
#' @author Catalina A. Vallejos \email{catalina.vallejos@@mrc-bsu.cam.ac.uk}
#' 
#' @rdname BASiCS_D_TestDE
BASiCS_D_TestDE <- function(Data, 
                             object,
                             GeneNames,
                             EpsilonM = 0.10,
                             EpsilonD = 0.10,
                             EviThresholdM = NULL,
                             EviThresholdD = NULL,
                             OrderVariable = "ProbDiffExp",
                             GroupLabelRef = "Ref",
                             GroupLabelTest = "Test",
                             Plot = FALSE, 
                             OffSet = FALSE, 
                             EFDR_M = 0.05,
                             EFDR_D = 0.05,
                             GenesSelect = NULL, 
                             ...)
{
  # Checking validity of input arguments
  if(!is(Data,"BASiCS_D_Data")) stop("'Data' is not a 'BASiCS_D_Data' class object.")    
  if(!is(object,"BASiCS_D_Chain")) stop("'object' is not a 'BASiCS_D_Chain' class object.")    
  if(EpsilonM < 0 | !is.finite(EpsilonM)) stop("Minimum tolerance of fold-change 'EpsilonM' must be a positive real value")
  if(EpsilonD < 0 | !is.finite(EpsilonD)) stop("Minimum tolerance of fold-change 'EpsilonD' must be a positive real value")
  if(!is.logical(Plot) | length(Plot)!=1) stop("Please insert T or F for Plot parameter")
  if(!is.null(EviThresholdM) | !is.null(EviThresholdD))
  {
    if(EviThresholdM < 0 | EviThresholdM > 1 | !is.finite(EviThresholdM)) 
      stop("Evidence threshold 'EviThresholdM' must be contained in (0,1) \n For automatic threshold search use EviThresholdM = NULL.")    
    if(EviThresholdD < 0 | EviThresholdD > 1 | !is.finite(EviThresholdD)) 
      stop("Evidence threshold 'EviThresholdD' must be contained in (0,1) \n For automatic threshold search use EviThresholdD = NULL.")    
  }
  if(!(OrderVariable %in% c("GeneIndex", "GeneNames", "ProbDiffExp", "ProbDiffOverDisp"))) stop("Invalid 'OrderVariable' value")
  if(!is.character(GroupLabelRef) | length(GroupLabelRef) > 1) stop("Invalid value for 'GroupLabelRef'")
  if(!is.character(GroupLabelTest) | length(GroupLabelTest) > 1) stop("Invalid value for 'GroupLabelTest'")
  if(!is.null(GenesSelect) & (length(GenesSelect) != length(GeneNames))) stop("Invalid value for 'GenesSelect'")
  if(!is.null(GenesSelect) & !is.logical(GenesSelect)) stop("Invalid value for 'GenesSelect'")
 # if(!is.null(GenesSelect) & (length(GenesSelect) == sum(!GenesSelect))) stop("Invalid value for 'GenesSelect'")
  
  nTest = ncol(Data@CountsTest)
  nRef = ncol(Data@CountsRef)
  n = nTest + nRef
  # Changes in overall expression
  if(OffSet)
  {
    # With offset correction
    ChainMuRefOffSet = object@muRef / rowSums(object@muRef)
    ChainMuTestOffSet = object@muTest / rowSums(object@muTest)
    MedianMuRefOffSet = apply(ChainMuRefOffSet, 2, median)
    MedianMuTestOffSet = apply(ChainMuTestOffSet, 2, median)
    
    ChainTau = log(ChainMuTestOffSet / ChainMuRefOffSet )
    MedianTau = apply(ChainTau, 2, median)
  }
  else
  {
    ChainTau = log(object@muTest / object@muRef)
    MedianTau = apply(ChainTau, 2, median)  
  }
  
  if(EpsilonM > 0) {ProbM = HiddenProbDE(chain = abs(ChainTau), tol = EpsilonM)}
  else 
  {
    ProbM_aux = HiddenProbDE(chain = ChainTau, tol = 0)
    ProbM = 2*pmax(ProbM_aux, 1-ProbM_aux) - 1
  }
  if(is.null(EviThresholdM))
  {
    EviThresholds <- seq(0.5,0.9995,by=0.00025)
    
    if(is.null(GenesSelect))
    {
      EFDR <- sapply(EviThresholds, HiddenEFDRDV, Prob = ProbM)
      EFNR <- sapply(EviThresholds, HiddenEFNRDV, Prob = ProbM)      
    }
    else
    {
      EFDR <- sapply(EviThresholds, HiddenEFDRDV, Prob = ProbM[GenesSelect])
      EFNR <- sapply(EviThresholds, HiddenEFNRDV, Prob = ProbM[GenesSelect])        
    }

    optimal = round(median(which(abs(EFDR - EFDR_M) == min(abs(EFDR - EFDR_M)))))
    OptThresholdM <- c(EviThresholds[optimal], EFDR[optimal], EFNR[optimal])

    EviThresholdM = OptThresholdM[1]
    if(is.na(EviThresholdM)) 
    { 
      cat("EFDR calibration failed. Probability threshold automatically set equal to 0.90 \n")
      EviThresholdM = 0.90
    }
  }
  else
  {
    if(is.null(GenesSelect))
    {
      EFDR = HiddenEFDRDV(EviThresholdM, ProbM)
      EFNR = HiddenEFNRDV(EviThresholdM, ProbM)
    }
    else
    {
      EFDR = HiddenEFDRDV(EviThresholdM, ProbM[GenesSelect])
      EFNR = HiddenEFNRDV(EviThresholdM, ProbM[GenesSelect])      
    }
    OptThresholdM <- c(EviThresholdM, EFDR, EFNR)
  }  
  
  # Changes in cell-to-cell biological over dispersion
  ChainOmega = log(object@deltaTest / object@deltaRef)
  MedianOmega = apply(ChainOmega, 2, median)
  if(EpsilonD > 0) {ProbD = HiddenProbDE(chain = abs(ChainOmega), tol = EpsilonD)}
  else 
  {
    ProbD_aux = HiddenProbDE(chain = ChainOmega, tol = 0)
    ProbD = 2*pmax(ProbD_aux, 1-ProbD_aux) - 1
  }
  if(is.null(EviThresholdD))
  {
    EviThresholds <- seq(0.5,0.9995,by=0.00025)
    
    if(is.null(GenesSelect))
    {
      EFDR <- sapply(EviThresholds, HiddenEFDRDV, Prob = ProbD)
      EFNR <- sapply(EviThresholds, HiddenEFNRDV, Prob = ProbD)
    }
    else
    {
      EFDR <- sapply(EviThresholds, HiddenEFDRDV, Prob = ProbD[GenesSelect])
      EFNR <- sapply(EviThresholds, HiddenEFNRDV, Prob = ProbD[GenesSelect])
    }
    
    optimal = round(median(which(abs(EFDR - EFDR_D) == min(abs(EFDR - EFDR_D)))))
    OptThresholdD <- c(EviThresholds[optimal], EFDR[optimal], EFNR[optimal])
    
    EviThresholdD = OptThresholdD[1]
    if(is.na(EviThresholdD)) 
    { 
      cat("EFDR calibration failed. Probability threshold automatically set equal to 0.90 \n")
      EviThresholdD = 0.90
    }
  }
  else
  {
    if(is.null(GenesSelect))
    {
      EFDR = HiddenEFDRDV(EviThresholdD, ProbD)
      EFNR = HiddenEFNRDV(EviThresholdD, ProbD)
    }
    else
    {
      EFDR = HiddenEFDRDV(EviThresholdD, ProbD[GenesSelect])
      EFNR = HiddenEFNRDV(EviThresholdD, ProbD[GenesSelect])      
    }
    OptThresholdD <- c(EviThresholdD, EFDR, EFNR)
  } 
  
  # Differential expression results
  TestPlusM = which(ProbM > EviThresholdM & MedianTau > 0)
  RefPlusM = which(ProbM > EviThresholdM & MedianTau < 0)  
  ResultDiffExp = rep("NoDiff", length(MedianTau))
  ResultDiffExp[TestPlusM] = paste0(GroupLabelTest,"+")
  ResultDiffExp[RefPlusM] = paste0(GroupLabelRef,"+") 
  if(!is.null(GenesSelect)) {ResultDiffExp[!GenesSelect] = "ExcludedByUser"}
  
  # Differential cell-to-cell biological over dispersion results
  TestPlusD = which(ProbD > EviThresholdD & MedianOmega > 0)
  RefPlusD = which(ProbD > EviThresholdD & MedianOmega < 0)
  ResultDiffOverDisp = rep("NoDiff", length(MedianTau))
  ResultDiffOverDisp[TestPlusD] = paste0(GroupLabelTest,"+")
  ResultDiffOverDisp[RefPlusD] = paste0(GroupLabelRef,"+")
  if(!is.null(GenesSelect)) {ResultDiffOverDisp[!GenesSelect] = "ExcludedByUser"}
  
  # Gene-specific parameters' summaries
  MedianMuTest = apply(object@muTest, 2, median)
  MedianMuRef = apply(object@muRef, 2, median)  
  MedianDeltaTest = apply(object@deltaTest, 2, median)
  MedianDeltaRef = apply(object@deltaRef, 2, median)
  MuBase=(MedianMuRef * nRef + MedianMuTest * nTest)/n
  DeltaBase=(MedianDeltaRef * nRef + MedianDeltaTest * nTest)/n
  
  FullList = cbind.data.frame("GeneNames" = GeneNames,
                              "ExpOverall"= round(as.numeric(MuBase),3),
                              "ExpTest"= round(as.numeric(MedianMuTest),3),
                              "ExpRef"= round(as.numeric(MedianMuRef),3),
                              "ExpFC"= round(as.numeric(exp(MedianTau)),3),
                              "ExpLogFC"= round(as.numeric(MedianTau),3),
                              "ProbDiffExp"= round(as.numeric(ProbM),3),
                              "ResultDiffExp" = ResultDiffExp,
                              "OverDispOverall"= round(as.numeric(DeltaBase),3),
                              "OverDispTest"= round(as.numeric(MedianDeltaTest),3),
                              "OverDispRef"= round(as.numeric(MedianDeltaRef),3),  
                              "OverDispFC"= round(as.numeric(exp(MedianOmega)),3), 
                              "OverDispLogFC"= round(as.numeric(MedianOmega),3),
                              "ProbDiffOverDisp"= round(as.numeric(ProbD),3),
                              "ResultDiffOverDisp" = ResultDiffOverDisp,
                              stringsAsFactors = FALSE)
  
  GeneIndex = 1:length(MuBase)
  
  if(OrderVariable == "GeneIndex") orderVar = GeneIndex
  if(OrderVariable == "GeneNames") orderVar = GeneNames
  if(OrderVariable == "ProbDiffExp") orderVar = ProbM
  if(OrderVariable == "ProbDiffOverDisp") orderVar = ProbD
  
  FullList = FullList[order(orderVar, decreasing = TRUE),]
  
  if(!is.null(GenesSelect))
  {
    cat("--------------------------------------------------------------------- \n")
    cat(paste("The user excluded ", sum(!GenesSelect), " genes from the comparison. \n"))
    cat("These genes are marked as 'ExcludedByUser' in the results table and excluded from EFDR calibration. \n")
    cat("--------------------------------------------------------------------- \n")    
  }

  
  cat("--------------------------------------------------------------------- \n")
  cat(paste(length(TestPlusM) + length(RefPlusM), " genes with a change on the overall expression:  \n"))
  cat(paste("- Higher expression in ",GroupLabelTest,"group:", length(TestPlusM), "\n"))
  cat(paste("- Higher expression in ",GroupLabelRef,"group:", length(RefPlusM), "\n"))
  cat(paste("- Fold change tolerance = ", round(100*EpsilonM,2), "% \n"))    
  cat(paste("- Evidence threshold = ", OptThresholdM[1], "\n"))
  cat(paste("- EFDR = ", round(100*OptThresholdM[2],2), "% \n"))  
  cat(paste("- EFNR = ", round(100*OptThresholdM[3],2), "% \n"))  
  cat("--------------------------------------------------------------------- \n")
  cat("\n")
  cat("--------------------------------------------------------------------- \n")
  cat(paste(length(TestPlusD) + length(RefPlusD), " genes with a change on the cell-to-cell biological over dispersion:  \n"))
  cat(paste("- Higher over dispersion in ",GroupLabelTest,"group:", length(TestPlusD), "\n"))
  cat(paste("- Higher over dispersion in ",GroupLabelRef,"group:", length(RefPlusD), "\n"))
  cat(paste("- Fold change tolerance = ", round(100*EpsilonD,2), "% \n"))    
  cat(paste("- Evidence threshold = ", OptThresholdD[1], "\n"))
  cat(paste("- EFDR = ", round(100*OptThresholdD[2],2), "% \n"))  
  cat(paste("- EFNR = ", round(100*OptThresholdD[3],2), "% \n"))    
  cat("--------------------------------------------------------------------- \n")
  
  if(Plot)
  {    
    args <- list(...)
    
#    if(Search)
#    {      
#      par(ask=T)
#      
#      plot(EviThresholds, EFDR, type = "l", lty = 1, bty = "n", ylab = "Error rate", xlab = "Evidence threshold", ylim = c(0,1))
#      lines(EviThresholds, EFNR, lty = 2)      
#      legend('topleft', c("EFDR", "EFNR"), lty = 1:2, bty = "n")
#    }
    
#    if("ylim" %in% names(args)) {ylim = args$ylim} else{ylim = c(0, 1)}
#    if("xlim" %in% names(args)) {xlim = args$xlim} else{xlim = c(min(Mu),max(Mu))}
#    cex = ifelse("cex" %in% names(args),args$cex, 1.5)
#    pch = ifelse("pch" %in% names(args),args$pch, 16)
#    col = ifelse("col" %in% names(args),args$col, 8)
#    bty = ifelse("bty" %in% names(args),args$bty, "n")
#    cex.lab = ifelse("cex.lab" %in% names(args),args$cex.lab, 1)
#    cex.axis = ifelse("cex.axis" %in% names(args),args$cex.axis, 1)
#    cex.main = ifelse("cex.main" %in% names(args),args$cex.main, 1) 
#    xlab = ifelse("xlab" %in% names(args),args$xlab, expression(mu[i]))
#    ylab = ifelse("ylab" %in% names(args),args$ylab, "HVG probability")
#    main = ifelse("main" %in% names(args),args$main, "") 
    
#    plot(Mu, Prob, log="x", pch = pch, ylim = ylim, xlim = xlim, col = col, cex = cex,
#         bty = bty, cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main,
#         xlab = xlab, ylab = ylab, main = main)
#    abline(h = OptThreshold[1], lty = 2, col = "black")
#    points(Mu[HVG], Prob[HVG], pch = pch, col = "red", cex = cex)
    
#    par(ask=F)
  }
  
  list("Table" = FullList, 
       "DiffExpSummary" = list("EviThreshold" = OptThresholdM[1],
                               "EFDR" = OptThresholdM[2], 
                               "EFNR" = OptThresholdM[3]),
       "DiffOverDispSummary" = list("EviThreshold" = OptThresholdD[1],
                                    "EFDR" = OptThresholdD[2], 
                                    "EFNR" = OptThresholdD[3]))
}