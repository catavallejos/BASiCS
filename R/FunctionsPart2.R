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
#' @author Catalina A. Vallejos \email{cnvallej@uc.cl}
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
#' @author Catalina A. Vallejos \email{cnvallej@uc.cl}
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
#' @author Catalina A. Vallejos \email{cnvallej@uc.cl}
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
#' @author Catalina A. Vallejos \email{cnvallej@uc.cl}
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
#' @author Catalina A. Vallejos \email{cnvallej@uc.cl}
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
