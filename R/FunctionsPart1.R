#' @title Filter for input datasets
#' 
#' @description \code{BASiCS_Filter} indicates which transcripts and cells pass a pre-defined inclusion criteria.
#' 
#' @param Counts Matrix of dimensions \code{q} times \code{n} whose elements corresponds to the simulated expression counts. 
#' First \code{q.bio} rows correspond to biological genes. Last \code{q-q.bio} rows correspond to technical spike-in genes. 
#' @param Tech Logical vector of length \code{q}. If \code{Tech = F} the gene is biological; otherwise the gene is spike-in.
#' @param SpikeInput Vector of length \code{q-q.bio} whose elements indicate the simulated input concentrations for the spike-in genes. 
#' @param MinTotalCountsPerCell Minimum value of total expression counts required per cell (biological and technical)
#' @param MinTotalCountsPerGene Minimum value of total expression counts required per transcript (biological and technical)
#' @param MinCellsWithExpression Minimum number of cells where expression must be detected (positive count). Criteria applied to each transcript.
#' @param MinAvCountsPerCellsWithExpression Minimum average number of counts per cells where expression is detected. Criteria applied to each transcript.
#' 
#' @return A list of 2 elements
#' \describe{
#' \item{\code{Counts}}{Filtered matrix of expression counts}
#' \item{\code{Tech}}{Filtered vector of spike-in indicators}
#' \item{\code{SpikeInput}}{Filtered vector of spike-in genes input molecules}
#' \item{\code{IncludeGenes}}{Inclusion indicators for transcripts}
#' \item{\code{IncludeCells}}{Inclusion indicators for cells}
#' }
#' 
#' @examples
#' 
#' Counts = matrix(rpois(50*10, 2), ncol = 10)
#' Tech = c(rep(FALSE,40),rep(TRUE,10))
#' SpikeInput = rgamma(10,1,1)
#' Filter = BASiCS_Filter(Counts, Tech, SpikeInput)
#' 
#' FilterData = newBASiCS_Data(Filter$Counts, Filter$Tech, Filter$SpikeInput)
#' 
#' @seealso \code{\link[BASiCS]{BASiCS_Data-class}}
#' 
#' @author Catalina A. Vallejos \email{catalina.vallejos@@mrc-bsu.cam.ac.uk}
#' 
#' @references Vallejos, Marioni and Richardson (2015). Bayesian Analysis of Single-Cell Sequencing data.
BASiCS_Filter <- function(Counts, Tech, SpikeInput, 
                          MinTotalCountsPerCell = 2, MinTotalCountsPerGene = 2, 
                          MinCellsWithExpression = 2, MinAvCountsPerCellsWithExpression = 2)
{
  q = length(Tech)
  q.bio = q - sum(SpikeInput)
  n = ncol(Counts)
  
  CellIndex = 1:n
  GeneIndex = 1:q
  
  # Remove cells with zero counts in either biological or technical genes
  IncludeCells = ifelse(apply(Counts[!Tech,],2,sum)>0 & apply(Counts[Tech,],2,sum)>0, T, F)
  if(sum(IncludeCells) == 0) stop('All cells have zero biological or technical counts (across all transcripts) \n')
  IncludeCells = ifelse(apply(Counts,2,sum) >= MinTotalCountsPerCell, IncludeCells, F)
  Counts1 = Counts[,IncludeCells]
  
  # Remove transcripts with zero counts across all cells
  IncludeBio = ifelse(apply(Counts1[!Tech,],1,sum) >= MinTotalCountsPerGene, T, F)
  IncludeTech = ifelse(apply(Counts1[Tech,],1,sum) >= MinTotalCountsPerGene, T, F)
  
  # Remove transcripts expressed in less than 'MinExpressedCells' cells
  NonZero <- I(Counts1>0)
  IncludeBio = ifelse(apply(NonZero[!Tech,], 1, sum) >= MinCellsWithExpression, IncludeBio, F)
  IncludeTech = ifelse(apply(NonZero[Tech,], 1, sum) >= MinCellsWithExpression, IncludeTech, F)
  
  # Remove transcripts with low counts in the cells where they are expressed  
  IncludeBio = ifelse(apply(Counts1[!Tech,], 1, sum)/apply(NonZero[!Tech,], 1, sum) >= MinAvCountsPerCellsWithExpression, IncludeBio, F)
  IncludeTech = ifelse(apply(Counts1[Tech,], 1, sum)/apply(NonZero[Tech,], 1, sum) >= MinAvCountsPerCellsWithExpression, IncludeTech, F)
    
  list("Counts" = Counts1[c(IncludeBio,IncludeTech),], "Tech" = Tech[c(IncludeBio,IncludeTech)],
       "SpikeInput" = SpikeInput[IncludeTech],
       "IncludeGenes" = c(IncludeBio,IncludeTech), "IncludeCells" = IncludeCells)
}

#' @title Creates a BASiCS_Data object from a matrix of expression counts and experimental information about spike-in genes
#' 
#' @description \code{newBASiCS_Data} creates a \code{\link[BASiCS]{BASiCS_Data-class}} object from 
#' a matrix of expression counts and experimental information about spike-in genes.
#' 
#' @param Counts Matrix of dimensions \code{q} times \code{n} whose elements corresponds to the simulated expression counts. 
#' First \code{q.bio} rows correspond to biological genes. Last \code{q-q.bio} rows correspond to technical spike-in genes. 
#' @param Tech Logical vector of length \code{q}. If \code{Tech = F} the gene is biological; otherwise the gene is spike-in.
#' @param SpikeInput Vector of length \code{q-q.bio} whose elements indicate the simulated input concentrations for the spike-in genes. 
#' @param BatchInfo Vector of length \code{n} whose elements indicate batch information. 
#' 
#' @return An object of class \code{\link[BASiCS]{BASiCS_Data-class}}.
#' 
#' @examples
#' 
#' set.seed(1)
#' Counts = matrix(rpois(10*5, 1), ncol = 5)
#' Tech = c(rep(FALSE,7),rep(TRUE,3))
#' set.seed(2)
#' SpikeInput = rgamma(3,1,1)
#' Data = newBASiCS_Data(Counts, Tech, SpikeInput)
#' 
#' @seealso \code{\link[BASiCS]{BASiCS_Data-class}}
#' 
#' @author Catalina A. Vallejos \email{catalina.vallejos@@mrc-bsu.cam.ac.uk}
#' 
#' @references Vallejos, Marioni and Richardson (2015). Bayesian Analysis of Single-Cell Sequencing data.
newBASiCS_Data <- function(Counts, Tech, SpikeInput, BatchInfo = NULL)
{
  if(is.null(BatchInfo)) {BatchInfo = rep(1, times = ncol(Counts))}
  Data <- new("BASiCS_Data", Counts = Counts, Tech = Tech, SpikeInput = SpikeInput, BatchInfo = BatchInfo)
  show(Data)
  cat('\n')
  cat('NOTICE: BASiCS requires a pre-filtered dataset \n')
  cat('    - You must remove poor quality cells before creating the BASiCS data object \n')
  cat('    - We recommend to pre-filter very lowly expressed transcripts before creating the object. \n')
  cat('      Inclusion criteria may vary for each data. For example, remove transcripts \n')
  cat('          - with very low total counts across of all samples \n')
  cat('          - that are only expressed in few cells \n')
  cat('            (by default genes expressed in only 1 cell are not accepted) \n')
  cat('          - with very low total counts across the samples where the transcript is expressed \n')
  cat('\n')
  cat(' BASiCS_Filter can be used for this purpose \n')
  return(Data)  
}

#' @title Simulates expression counts according to the model implemented in BASiCS 
#' 
#' @description \code{BASiCS_Sim} creates a simulated dataset from the model implemented in BASiCS. 
#' This function is used in order to illustrate the performance of the \code{BASiCS} library.
#' 
#' @param mu Gene-specific expression levels \eqn{\mu[i]}, defined as true input molecules in case of technical genes 
#' (vector of length \code{q}, technical genes located at the end of the vector, all elements must be positive numbers)
#' @param delta Gene-specific biological cell-to-cell heterogeneity hyper-parameters \eqn{\delta[i]}, biological genes only 
#' (vector of length \code{q.bio}, all elements must be positive numbers)
#' @param phi Cell-specific mRNA content normalising constants \eqn{\phi[j]}
#' (vector of length \code{n}, all elements must be positive numbers and the sum of its elements must be equal to \code{n})
#' @param s Cell-specific capture efficiency (or amplification biases if not using UMI based counts) normalising constants \eqn{s[j]}
#' (vector of length \code{n}, all elements must be positive numbers)
#' @param theta Technical variability hyper-parameter \eqn{\theta} (must be positive)
#' 
#' @return An object of class \code{\link[BASiCS]{BASiCS_Data-class}}, simulated from the model implemented in BASiCS.
#'  
#' @examples
#' 
#' # Simulated parameter values for 10 genes 
#' # (7 biogical and 5 spike-in) measured in 5 cells
#' Mu =  c(8.36, 10.65, 4.88, 6.29, 21.72, 12.93, 30.19, 1010.72, 7.90, 31.59)
#' Delta = c(1.29, 0.88, 1.51, 1.49, 0.54, 0.40, 0.85)
#' Phi = c(1.00, 1.06, 1.09, 1.05, 0.80)
#' S = c(0.38, 0.40, 0.38, 0.39, 0.34)
#' Theta = 0.39
#' 
#' Data = BASiCS_Sim(Mu, Delta, Phi, S, Theta)
#' head(counts(Data))
#' dim(counts(Data, type="biological"))
#' dim(counts(Data, type="technical"))
#' displayTechIndicator(Data)
#' displaySpikeInput(Data)
#' 
#' @seealso \code{\link[BASiCS]{BASiCS_Data-class}}
#' 
#' @author Catalina A. Vallejos \email{catalina.vallejos@@mrc-bsu.cam.ac.uk}
#' 
#' @references Vallejos, Marioni and Richardson (2015). Bayesian Analysis of Single-Cell Sequencing data.
BASiCS_Sim<-function(
  mu, 
  delta, 
  phi, 
  s, 
  theta) 
{
  # Number of cells 
  n = length(phi)
  # Total number of genes, including biological and technical ones
  q = length(mu)
  # Number of biological genes
  q.bio = length(delta)
  
  # Arguments checking
  if(!(is.vector(mu) & is.numeric(mu) & all(mu>0))) stop("Invalid argument value for 'mu'.")
  if(!(is.vector(delta) & is.numeric(delta) & all(delta>=0))) stop("Invalid argument value for 'delta'.")
  if(!(is.vector(phi) & is.numeric(phi) & all(phi>0) & all.equal(sum(phi),n))) stop("Invalid argument value for 'phi'.")
  if(!(is.vector(s) & is.numeric(s) & all(s>0) & length(s)==n)) stop("Invalid argument value for 's'.")
  if(!(is.numeric(theta) & length(theta)==1 & theta>=0)) stop("Invalid argument value for 'theta'.")
  
  if(!all(c(n,q,q.bio,q-q.bio)>0)) stop("Arguments' dimensions are not compatible")
  
  # Matrix where simulated counts will be stored
  Counts.sim<-matrix(0,ncol=n,nrow=q)
  # Matrix where gene-cell specific simulated random effects will be stored 
  rho<-matrix(1,ncol=n,nrow=q.bio)
  # Simulated cell-specific random effects
  if(theta>0) {nu<-rgamma(n,shape=1/theta,rate=1/(s*theta))}
  else {nu<-s}
  # Simulated counts data
  for(i in 1:q)
  {
    # Biological genes
    if(i<=q.bio)
    {
      if(delta[i]>0){rho[i,]<-rgamma(n,shape=1/delta[i],rate=1/delta[i])}
      Counts.sim[i,]<-rpois(n,lambda=phi*nu*rho[i,]*mu[i])
    }
    # Technical genes
    else {Counts.sim[i,]<-rpois(n,lambda=nu*mu[i])}
  }  
  
  Data = newBASiCS_Data(Counts = Counts.sim, Tech = ifelse(1:q > q.bio, T, F), SpikeInput = mu[(q.bio+1):q])
  
  return(Data)
}

#' @name makeExampleBASiCS_Data
#' 
#' @title Create a simple example of a BASiCS_Data object with random data
#' 
#' @description A simple \code{\link[BASiCS]{BASiCS_Data-class}} object is generated by simulating a dataset from the
#' model in BASiCS (Vallejos et al 2015). This is used for the examples in the package and the vignette. 
#' 
#' @param WithBatch If true, 2 batches are generated (each of them containing 10 cells)
#' 
#' @return An object of class \code{\link[BASiCS]{BASiCS_Data-class}}, simulated from the model implemented in BASiCS. 
#' It contains 120 genes (100 biological and 20 spike-in) and 20 cells. 
#' 
#' @examples 
#' Data = makeExampleBASiCS_Data()
#' is(Data, "BASiCS_Data")
#' 
#' @author Catalina A. Vallejos \email{catalina.vallejos@@mrc-bsu.cam.ac.uk}
#' 
#' @references Vallejos, Marioni and Richardson (2015). Bayesian Analysis of Single-Cell Sequencing data.  
makeExampleBASiCS_Data <- function(WithBatch = FALSE)
{
  Mu =  c( 8.36,  10.65,   4.88,   6.29,  21.72,  12.93,  30.19,  83.92,   3.89,   6.34,  
           57.87,  12.45,   8.08,   7.31,  15.56,  15.91,  12.24,  15.96,  19.14,   4.20,   
           6.75,  27.74,   8.88,  21.21,  19.89,   7.14,  11.09,   7.19,  20.64,  73.90,   
           9.05,   6.13,  16.40,   6.08,  17.89,   6.98,  10.25,  14.05,   8.14,   5.67,   
           6.95,  11.16,  11.83,   7.56, 159.05,  16.41,   4.58,  15.46,  10.96,  25.05,  
           18.54,   8.50,   4.05,   5.37,   4.63,   4.08,   3.75,   5.02,  27.74,  10.28,   
           3.91,  13.10,   8.23,   3.64,  80.77,  20.36,   3.47,  20.12,  13.29,   7.92,  
           25.51, 173.01,  27.71,   4.89,  33.10,   3.42,   6.19,   4.29,   5.19,   8.36,  
           10.27,   6.86,   5.72,  37.25,   3.82,  23.97,   5.80,  14.30,  29.07,   5.30,   
           7.47,   8.44,   4.24,  16.15,  23.39, 120.22,   8.92,  97.15,   9.75,  10.07, 
           1010.72,   7.90,  31.59,  63.17,   1.97, 252.68,  31.59,  31.59,  31.59,  63.17, 
           4042.89,4042.89,   3.95, 126.34, 252.68,   7.90,  15.79, 126.34,   3.95,  15.79)
  Delta = c(1.29, 0.88, 1.51, 1.49, 0.54, 0.40, 0.85, 0.27, 0.53, 1.31, 
            0.26, 0.81, 0.72, 0.70, 0.96, 0.58, 1.15, 0.82, 0.25, 5.32, 
            1.13, 0.31, 0.66, 0.27, 0.76, 1.39, 1.18, 1.57, 0.55, 0.17, 
            1.40, 1.47, 0.57, 2.55, 0.62, 0.77, 1.47, 0.91, 1.53, 2.89, 
            1.43, 0.77, 1.37, 0.57, 0.15, 0.33, 3.99, 0.47, 0.87, 0.86,
            0.97, 1.25, 2.20, 2.19, 1.26, 1.89, 1.70, 1.89, 0.69, 1.63, 
            2.83, 0.29, 1.21, 2.06, 0.20, 0.34, 0.71, 0.61, 0.71, 1.20, 
            0.88, 0.17, 0.25, 1.48, 0.31, 2.49, 2.75, 1.43, 2.65, 1.97,
            0.84, 0.81, 2.75, 0.53, 2.23, 0.45, 1.87, 0.74, 0.53, 0.58, 
            0.94, 0.72, 2.61, 1.56, 0.37, 0.07, 0.90, 0.12, 0.76, 1.45)
  Phi = c(1.09, 1.16, 1.19, 1.14, 0.87, 1.10, 0.48, 1.06, 0.94, 0.97,
          1.09, 1.16, 1.19, 1.14, 0.87, 1.10, 0.48, 1.06, 0.94, 0.97)
  S = c(0.38, 0.40, 0.38, 0.39, 0.34, 0.39, 0.31, 0.39, 0.40, 0.37,
        0.38, 0.40, 0.38, 0.39, 0.34, 0.39, 0.31, 0.39, 0.40, 0.37)
  
  q = length(Mu); q.bio = length(Delta); n = length(Phi)
  
  if(!WithBatch) {Theta = 0.5}
  else
  {
    # 2 batches, 10 cells each
    Theta1 = 0.5; Theta2 = 0.75; 
    Theta = ifelse(1:n <= 10, Theta1, Theta2)
  }
  
  # Matrix where simulated counts will be stored
  Counts.sim<-matrix(0,ncol=n,nrow=q)
  # Matrix where gene-cell specific simulated random effects will be stored 
  Rho<-matrix(1,ncol=n,nrow=q.bio)
  # Simulated cell-specific random effects
  if(all(Theta>0)) {set.seed(1000); Nu<-rgamma(n,shape=1/Theta,rate=1/(S*Theta))}
  else {Nu<-S}
  # Simulated counts data
  for(i in 1:q)
  {
    # Biological genes
    if(i<=q.bio)
    {
      if(Delta[i]>0){set.seed(i); Rho[i,]<-rgamma(n,shape=1/Delta[i],rate=1/Delta[i])}
      set.seed(i+10000); Counts.sim[i,]<-rpois(n,lambda=Phi*Nu*Rho[i,]*Mu[i])
    }
    # Technical genes
    else {set.seed(i+20000); Counts.sim[i,]<-rpois(n,lambda=Nu*Mu[i])}
  }  
  
  if(!WithBatch)
  {
    Data = newBASiCS_Data(Counts = Counts.sim, Tech = ifelse(1:q > q.bio, T, F), SpikeInput = Mu[(q.bio+1):q])  
  }
  else
  {
    Data = newBASiCS_Data(Counts = Counts.sim, Tech = ifelse(1:q > q.bio, T, F), 
                          SpikeInput = Mu[(q.bio+1):q], BatchInfo = c(rep(1,10), rep(2,10)))
  }
  
  
  return(Data)  
}


HiddenBASiCS_MCMC_Start<-function(
  Data,
  ...)
{
  if(!is(Data,"BASiCS_Data")) stop("'Data' is not a BASiCS_Data class object.")
  
  # Number of instrinsic genes
  q.bio<-sum(Data@Tech==F)
  # Number of cells
  n <- dim(Data@Counts)[2]
  
  # Initialize phi as a vector of ones
  phi0 = rep(1, times = n)
  # Initialize s as the empirical capture efficiency rates
  s0=colSums(Data@Counts[Data@Tech,])/sum(Data@SpikeInput); nu0=s0 
  # Initialize mu using average 'normalised counts' across cells 
  # (correcting by empirical capture efficiency rates)
  # and true input values for spike-in genes
  nCountsBio <- t( t(Data@Counts[!Data@Tech,]) / s0 )
  meansBio <- rowMeans( nCountsBio )
  mu0<-c(meansBio,Data@SpikeInput)
  
  # Random stating value for delta
  delta0=rgamma(q.bio,1,1)
  
  # Random stating value for theta (within typically observed range)
  theta0=runif(1, min = 0.2, max = 1)
  
  # If given, load default values for adaptive proposal variances
  args <- list(...)
  ls.mu0 = ifelse("ls.mu0" %in% names(args),args$ls.mu0,-4)
  ls.delta0 = ifelse("ls.delta0" %in% names(args),args$ls.delta0,-3)
  ls.phi0 = ifelse("ls.phi0" %in% names(args),args$ls.phi0,0)
  ls.nu0 = ifelse("ls.nu0" %in% names(args),args$ls.nu0,-10)
  ls.theta0 = ifelse("ls.theta0" %in% names(args),args$ls.theta0,-6)
  
  # Starting values for the proposal variances
  ls.mu0 =  pmax(2 * log (0.02 * abs(log(mu0))),ls.mu0)
  ls.delta0 =  pmax(2 * log (0.02 * abs(log(delta0))),ls.delta0)
  ls.phi0 = ifelse(n<200, pmax(2*log(n),ls.phi0), 11) # 6
  ls.nu0 =  pmax(2 * log (0.02 * abs(log(nu0))),ls.nu0)
  ls.theta0 =  pmax(2 * log (0.02 * abs(log(theta0))),ls.theta0)
  
  return(list("mu0"=mu0, "delta0"=delta0, "phi0"=phi0, "s0"=s0, "nu0"=nu0, "theta0"=theta0,
              "ls.mu0"=ls.mu0, "ls.delta0"=ls.delta0, "ls.phi0"=ls.phi0, "ls.nu0"=ls.nu0, "ls.theta0"=ls.theta0))
}

#' @title BASiCS MCMC sampler
#' 
#' @description MCMC sampler to perform Bayesian inference for single-cell mRNA sequencing datasets using the model described in Vallejos et al (2015). 
#' 
#' @param Data A \code{\link[BASiCS]{BASiCS_Data-class}} object.
#' @param N Total number of iterations for the MCMC sampler. Use \code{N>=max(4,Thin)}, \code{N} being a multiple of \code{Thin}. 
#' @param Thin Thining period for the MCMC sampler. Use \code{Thin>=2}.
#' @param Burn Burn-in period for the MCMC sampler. Use \code{Burn>=1}, \code{Burn<N}, \code{Burn} being a multiple of \code{Thin}.
#' @param ... Optional parameters.
#' \describe{
#' \item{\code{PriorParam}}{List of 7 elements, containing the hyper-parameter values required for the adopted prior (see Vallejos et al, 2015). All elements must be positive real numbers. 
#' \describe{
#'   \item{\code{a.delta}}{Shape hyper-parameter for the Gamma(\code{a.delta},\code{b.delta}) prior that is shared by all gene-specific biological cell-to-cell heterogeneity hyper-parameters \eqn{\delta[i]}. 
#'   Default: \code{a.delta = 1}.}
#'   \item{\code{b.delta}}{Rate hyper-parameter for the Gamma(\code{a.delta},\code{b.delta}) prior that is shared by all gene-specific biological cell-to-cell heterogeneity hyper-parameters \eqn{\delta[i]}.
#'   Default: \code{b.delta = 1}.}
#'   \item{\code{p.phi}}{Dirichlet hyper-parameter for the joint of all (scaled by \code{n}) cell-specific mRNA content normalising constants \eqn{\phi[j] / n}.
#'   Default: \code{p.phi = rep(1, n)}.}
#'   \item{\code{a.s}}{Shape hyper-parameter for the Gamma(\code{a.s},\code{b.s}) prior that is shared by all cell-specific capture efficiency normalising constants \eqn{s[j]}.
#'   Default: \code{a.s = 1}.}
#'   \item{\code{b.s}}{Rate hyper-parameter for the Gamma(\code{a.s},\code{b.s}) prior that is shared by all cell-specific capture efficiency normalising constants \eqn{s[j]}.
#'   Default: \code{b.s = 1}.}
#'   \item{\code{a.theta}}{Shape hyper-parameter for the Gamma(\code{a.theta},\code{b.theta}) prior for technical noise hyper-parameter \eqn{\theta}.
#'   Default: \code{a.theta = 1}.}
#'   \item{\code{b.theta}}{Rate hyper-parameter for the Gamma(\code{a.theta},\code{b.theta}) prior for technical noise hyper-parameter \eqn{\theta}.
#'   Default: \code{b.theta = 1}.}
#' }}
#' \item{\code{AR}}{Optimal acceptance rate for adaptive Metropolis Hastings updates. It must be a positive number between 0 and 1. Default (and recommended): \code{ar = 0.44}}.
#' 
#' \item{\code{StopAdapt}}{Iteration at which adaptive proposals are not longer adapted. Use \code{stopAdapt>=1}. Default: \code{StopAdapt = Burn}.}
#' 
#' \item{\code{StoreChains}}{If \code{StoreChains = T}, MCMC chains of each parameter are stored in separate .txt files. (\code{RunName} argument used for file names). Default: \code{StoreChains = F}.} 
#' \item{\code{StoreAdapt}}{If \code{StoreAdapt = T}, trajectory of adaptive proposal variances (log scale) for each parameter are stored in separate .txt files. (\code{RunName} argument used for file names). Default: \code{StoreAdapt = F}.}
#' \item{\code{StoreDir}}{Directory where MCMC chain will be stored (only required if \code{Store = T}). Default: \code{StoreDir = getwd()}.}
#' \item{\code{RunName}}{Run-name to be used when storing chains and/or adaptive proposal variances in .txt files.} 
#' \item{\code{PrintProgress}}{If \code{PrintProgress = T}, intermediate output is displayed in the console. }
#' }
#' 
#' @return An object of class \code{\link[BASiCS]{BASiCS_Chain-class}}. 
#' 
#' @examples
#' 
#' # Built-in simulated dataset
#' Data = makeExampleBASiCS_Data()
#' 
#' # Only a short run of the MCMC algorithm for illustration purposes
#' # Longer runs migth be required to reach convergence
#' MCMC_Output <- BASiCS_MCMC(Data, N = 10000, Thin = 10, Burn = 5000, PrintProgress = FALSE)
#' head(displayChainBASiCS(MCMC_Output, Param = "mu"))
#' head(displayChainBASiCS(MCMC_Output, Param = "delta"))
#' head(displayChainBASiCS(MCMC_Output, Param = "phi"))
#' head(displayChainBASiCS(MCMC_Output, Param = "s"))
#' head(displayChainBASiCS(MCMC_Output, Param = "nu"))
#' head(displayChainBASiCS(MCMC_Output, Param = "theta"))
#' 
#' # Traceplots
#' plot(MCMC_Output, Param = "mu", Gene = 1)
#' plot(MCMC_Output, Param = "delta", Gene = 1)
#' plot(MCMC_Output, Param = "phi", Cell = 1)
#' plot(MCMC_Output, Param = "s", Cell = 1)
#' plot(MCMC_Output, Param = "nu", Cell = 1)
#' plot(MCMC_Output, Param = "theta", Batch = 1)
#' 
#' # Calculating posterior medians and 95% HPD intervals
#' MCMC_Summary <- Summary(MCMC_Output)
#' head(displaySummaryBASiCS(MCMC_Summary, Param = "mu"))
#' head(displaySummaryBASiCS(MCMC_Summary, Param = "delta"))
#' head(displaySummaryBASiCS(MCMC_Summary, Param = "phi"))
#' head(displaySummaryBASiCS(MCMC_Summary, Param = "s"))
#' head(displaySummaryBASiCS(MCMC_Summary, Param = "nu"))
#' head(displaySummaryBASiCS(MCMC_Summary, Param = "theta"))
#' 
#' # Graphical display of posterior medians and 95% HPD intervals
#' plot(MCMC_Summary, Param = "mu", main = "All genes")
#' plot(MCMC_Summary, Param = "mu", Genes = 1:10, main = "First 10 genes")
#' plot(MCMC_Summary, Param = "delta", main = "All genes")
#' plot(MCMC_Summary, Param = "delta", Genes = c(2,5,10,50,100), main = "5 customized genes")
#' plot(MCMC_Summary, Param = "phi", main = "All cells")
#' plot(MCMC_Summary, Param = "phi", Cells = 1:5, main = "First 5 cells")
#' plot(MCMC_Summary, Param = "s", main = "All cells")
#' plot(MCMC_Summary, Param = "s", Cells = 1:5, main = "First 5 cells")
#' plot(MCMC_Summary, Param = "nu", main = "All cells")
#' plot(MCMC_Summary, Param = "nu", Cells = 1:5, main = "First 5 cells")
#' plot(MCMC_Summary, Param = "theta")
#' 
#' # To constrast posterior medians of cell-specific parameters
#' plot(MCMC_Summary, Param = "phi", Param2 = "s")
#' plot(MCMC_Summary, Param = "phi", Param2 = "nu")
#' plot(MCMC_Summary, Param = "s", Param2 = "nu")
#' 
#' # To constrast posterior medians of gene-specific parameters
#' plot(MCMC_Summary, Param = "mu", Param2 = "delta", log = "x")
#' 
#' # Highly and lowly variable genes detection
#' #DetectHVG <- BASiCS_DetectHVG(Data, MCMC_Output, VarThreshold = 0.70, Plot = TRUE)
#' #DetectLVG <- BASiCS_DetectLVG(Data, MCMC_Output, VarThreshold = 0.40, Plot = TRUE)
#' 
#' #plot(MCMC_Summary, Param = "mu", Param2 = "delta", log = "x", col = 8)
#' #points(DetectHVG$Table[DetectHVG$Table$HVG==1,2], DetectHVG$Table[DetectHVG$Table$HVG==1,3], 
#' #       pch = 16, col = "red", cex = 1)
#' #points(DetectLVG$Table[DetectLVG$Table$LVG==1,2], DetectLVG$Table[DetectLVG$Table$LVG==1,3], 
#' #       pch = 16, col = "blue", cex = 1)
#' 
#' # If variance thresholds are not fixed
#' #BASiCS_VarThresholdSearchHVG(Data, MCMC_Output, VarThresholdsGrid = seq(0.70,0.75,by=0.01))
#' #BASiCS_VarThresholdSearchLVG(Data, MCMC_Output, VarThresholdsGrid = seq(0.40,0.45,by=0.01))
#' 
#' @author Catalina A. Vallejos \email{catalina.vallejos@@mrc-bsu.cam.ac.uk}
#' 
#' @references Vallejos, Marioni and Richardson (2015). Bayesian Analysis of Single-Cell Sequencing data.
BASiCS_MCMC <- function(
  Data, 
  N, 
  Thin, 
  Burn,  
  ...)
{
  
  if(!is(Data,"BASiCS_Data")) stop("'Data' is not a BASiCS_Data class object.")
  
  # SOME QUANTITIES USED THROUGHOUT THE MCMC ALGORITHM
  q=length(Data@Tech); q.bio=length(Data@Tech[!Data@Tech]); n=dim(Data@Counts)[2]
  
  args <- list(...)
  if("PriorParam" %in% names(args)) {PriorParam = args$PriorParam}
  else { PriorParam = list(a.delta = 1, b.delta = 1, p.phi = rep(1, times = n), a.s = 1, b.s = 1, a.theta = 1, b.theta = 1)}
  AR = ifelse("AR" %in% names(args),args$AR, 0.44)
  StopAdapt = ifelse("StopAdapt" %in% names(args),args$StopAdapt, Burn)
  StoreChains = ifelse("StoreChains" %in% names(args),args$StoreChains, F)
  StoreAdapt = ifelse("StoreAdapt" %in% names(args),args$StoreAdapt, F)
  StoreDir = ifelse("StoreDir" %in% names(args),args$StoreDir, getwd())  
  RunName = ifelse("RunName" %in% names(args),args$RunName, "")
  PrintProgress = ifelse("PrintProgress" %in% names(args),args$PrintProgress, TRUE)
  
  if(!(length(N) == 1 | length(Thin) == 1 | length(Burn) == 1)) stop("Invalid parameter values.")
  if(!(N%%Thin==0 & N>=max(4,Thin))) stop("Please use an integer value for N. It must also be a multiple of thin (N>=4)).")
  if(!(Thin%%1==0 & Thin>=2)) stop("Please use an integer value for Thin (Thin>=2).") 
  if(!(Burn%%Thin==0 & Burn<N & Burn>=1)) stop("Please use an integer value for Burn. It must also be lower than N and a multiple of thin (Burn>=1).")
  
  if(!(PriorParam$a.delta>0  & length(PriorParam$a.delta) == 1 &
         PriorParam$b.delta>0  & length(PriorParam$b.delta) == 1 &
         all(PriorParam$p.phi>0) & length(PriorParam$p.phi) == n &
         PriorParam$a.s>0      & length(PriorParam$a.s) == 1 &
         PriorParam$b.s>0      & length(PriorParam$b.s) == 1 &
         PriorParam$a.theta>0  & length(PriorParam$a.theta) == 1 &
         PriorParam$b.theta>0) & length(PriorParam$b.theta) == 1) stop("Invalid prior hyper-parameter values.")
  
  if(!(AR>0 & AR<1 & length(AR) == 1)) stop("Invalid AR value. Recommended value: AR = 0.44.")
  if(!(StopAdapt>0)) stop("Invalid StopAdapt value.")
  if(!(is.logical(StoreChains) & length(StoreChains) == 1)) stop("Invalid StoreChains value.")
  if(!(is.logical(StoreAdapt) & length(StoreAdapt) == 1)) stop("Invalid StoreAdapt value.")
  if(!(file.info(StoreDir)["isdir"])) stop("Invalid StoreDir value.")
  
  # SOME SUMS USED THROUGHOUT THE MCMC ALGORITHM
  sum.bycell.all<-apply(Data@Counts,1,sum)
  sum.bycell.bio<-apply(Data@Counts[1:q.bio,],1,sum)
  sum.bygene.all<-apply(Data@Counts,2,sum)
  sum.bygene.bio<-apply(Data@Counts[1:q.bio,],2,sum)
  
  # GENERATING STARTING VALUES 
  Start=HiddenBASiCS_MCMC_Start(Data)
  # Starting values for MCMC chains
  mu0=as.vector(Start$mu0); delta0=as.vector(Start$delta0)
  phi0=as.vector(Start$phi0); s0=as.vector(Start$s0)
  nu0=as.vector(Start$nu0); theta0=as.numeric(Start$theta0)
  # Starting values for adaptive proposal variances
  ls.mu0=as.vector(Start$ls.mu0); ls.delta0=as.vector(Start$ls.delta0)
  ls.phi0=as.numeric(Start$ls.phi0) 
  ls.nu0=as.vector(Start$ls.nu0); ls.theta0=as.numeric(Start$ls.theta0)
  
  StoreAdaptNumber = as.numeric(StoreAdapt)
  nBatch = length(unique(Data@BatchInfo))
  
  if(nBatch > 1)
  {
    BatchDesign = model.matrix(~as.factor(Data@BatchInfo)-1)
    
    # MCMC SAMPLER (FUNCTION IMPLEMENTED IN C++)
    Time = system.time(Chain <- HiddenBASiCS_MCMCcppBatch(
      N,
      Thin,
      Burn,
      as.matrix(Data@Counts),
      BatchDesign,
      mu0, delta0, phi0, s0, nu0, theta0,
      PriorParam$a.delta, PriorParam$b.delta,
      PriorParam$p.phi,
      PriorParam$a.s, PriorParam$b.s,
      PriorParam$a.theta, PriorParam$b.theta,
      AR,
      ls.mu0, ls.delta0, ls.phi0, ls.nu0, ls.theta0,
      sum.bycell.all, sum.bycell.bio, sum.bygene.all, sum.bygene.bio,
      StoreAdaptNumber,StopAdapt,as.numeric(PrintProgress)))  
  }
  else
  {
    # MCMC SAMPLER (FUNCTION IMPLEMENTED IN C++)
    Time = system.time(Chain <- HiddenBASiCS_MCMCcpp(
      N,
      Thin,
      Burn,
      as.matrix(Data@Counts),
      mu0, delta0, phi0, s0, nu0, theta0,
      PriorParam$a.delta, PriorParam$b.delta,
      PriorParam$p.phi,
      PriorParam$a.s, PriorParam$b.s,
      PriorParam$a.theta, PriorParam$b.theta,
      AR,
      ls.mu0, ls.delta0, ls.phi0, ls.nu0, ls.theta0,
      sum.bycell.all, sum.bycell.bio, sum.bygene.all, sum.bygene.bio,
      StoreAdaptNumber,StopAdapt,as.numeric(PrintProgress)))   
  }
  
  cat("--------------------------------------------------------------------"); cat("\n")
  cat("MCMC running time"); cat("\n")
  cat("--------------------------------------------------------------------"); cat("\n")
  print(Time)
  cat("\n")
  
  OldDir = getwd()
  
  if(StoreChains)
  {   
    setwd(StoreDir)
    
    cat("--------------------------------------------------------------------"); cat("\n")
    cat("Storing MCMC chains of model parameters as .txt files in"); cat("\n")
    cat(paste0("'",StoreDir,"' directory ... ")); cat("\n")
    cat("--------------------------------------------------------------------"); cat("\n")
    
    write.table(Chain$mu,paste0("chain_mu_",RunName,".txt"),col.names=F,row.names=F)
    write.table(Chain$delta,paste0("chain_delta_",RunName,".txt"),col.names=F,row.names=F)
    write.table(Chain$phi,paste0("chain_phi_",RunName,".txt"),col.names=F,row.names=F)
    write.table(Chain$s,paste0("chain_s_",RunName,".txt"),col.names=F,row.names=F) 
    write.table(Chain$nu,paste0("chain_nu_",RunName,".txt"),col.names=F,row.names=F)
    write.table(Chain$theta,paste0("chain_theta_",RunName,".txt"),col.names=F,row.names=F)
    
    setwd(OldDir)
  }
  
  if(StoreAdapt)
  {
    setwd(StoreDir)
    
    cat("--------------------------------------------------------------------"); cat("\n")
    cat("Storing trajectories of adaptive proposal variances (log-scale) as .txt files in"); cat("\n")
    cat(paste0("'",StoreDir,"' directory ... ")); cat("\n")
    cat("--------------------------------------------------------------------"); cat("\n")
    
    write.table(Chain$ls.mu,paste0("chain_ls.mu_",RunName,".txt"),col.names=F,row.names=F)
    write.table(Chain$ls.delta,paste0("chain_ls.delta_",RunName,".txt"),col.names=F,row.names=F)      
    write.table(Chain$ls.phi,paste0("chain_ls.phi_",RunName,".txt"),col.names=F,row.names=F)
    write.table(Chain$ls.nu,paste0("chain_ls.nu_",RunName,".txt"),col.names=F,row.names=F)
    write.table(Chain$ls.theta,paste0("chain_ls.theta_",RunName,".txt"),col.names=F,row.names=F)
    
    setwd(OldDir)
  }
  
  cat("--------------------------------------------------------------------"); cat("\n")
  cat("Output"); cat("\n")
  cat("--------------------------------------------------------------------"); cat("\n")
  
  ChainClass <- newBASiCS_Chain(mu = Chain$mu, delta = Chain$delta, phi = Chain$phi, 
                                s = Chain$s, nu = Chain$nu, theta = Chain$theta)
  
  return(ChainClass)
}

#' @title Creates a BASiCS_Chain object from pre-computed MCMC chains
#' 
#' @description \code{BASiCS_Chain} creates a \code{\link[BASiCS]{BASiCS_Chain-class}} object from pre-computed MCMC chains.
#' 
#' @param mu MCMC chain for gene-specific expression levels \eqn{\mu[i]}, defined as true input molecules in case of technical genes 
#' (matrix with \code{q} columns, technical genes located at the end of the matrix, all elements must be positive numbers)
#' @param delta MCMC chain for gene-specific biological cell-to-cell heterogeneity hyper-parameters \eqn{\delta[i]}, biological genes only 
#' (matrix with \code{q.bio} columns, all elements must be positive numbers)
#' @param phi MCMC chain for cell-specific mRNA content normalising constants \eqn{\phi[j]}
#' (matrix with \code{n} columns, all elements must be positive numbers and the sum of its elements must be equal to \code{n})
#' @param s MCMC chain for cell-specific capture efficiency (or amplification biases if not using UMI based counts) normalising constants \eqn{s[j]}
#' (matrix with \code{n} columns, all elements must be positive numbers)
#' @param nu MCMC chain for cell-specific random effects \eqn{\nu[j]}
#' (matrix with \code{n} columns, all elements must be positive numbers)
#' @param theta MCMC chain for technical variability hyper-parameter \eqn{\theta} (vector, all elements must be positive) 
#' 
#' @return An object of class \code{\link[BASiCS]{BASiCS_Chain-class}}. 
#' 
#' @examples
#' 
#' # Data = makeExampleBASiCS_Data()
#' # MCMC_Output <- BASiCS_MCMC(Data, N = 50, Thin = 5, Burn = 5, 
#' #                StoreChains = TRUE, StoreDir = getwd(), RunName = "Test")
#' 
#' # ChainMu = as.matrix(read.table("chain_mu_Test.txt"))
#' # ChainDelta = as.matrix(read.table("chain_delta_Test.txt"))
#' # ChainPhi = as.matrix(read.table("chain_phi_Test.txt"))
#' # ChainS = as.matrix(read.table("chain_s_Test.txt"))
#' # ChainNu = as.matrix(read.table("chain_nu_Test.txt"))#
#' # ChainTheta = read.table("chain_theta_Test.txt")[,1]
#' 
#' # MCMC_Output_Load <- newBASiCS_Chain(mu = ChainMu, delta = ChainDelta, 
#' #   phi = ChainPhi, s = ChainS, nu = ChainNu, theta = ChainTheta)
#' 
#' @seealso \code{\link[BASiCS]{BASiCS_Chain-class}}
#' 
#' @author Catalina A. Vallejos \email{catalina.vallejos@@mrc-bsu.cam.ac.uk}
#' 
#' @references Vallejos, Marioni and Richardson (2015). Bayesian Analysis of Single-Cell Sequencing data.
newBASiCS_Chain <- function(mu, delta, phi, s, nu, theta)
{
  Chain <- new("BASiCS_Chain", mu = mu, delta = delta, phi = phi, s = s, nu = nu, theta = theta)
  show(Chain)
  return(Chain)  
}

##########################################################################################
# Highly variable genes detection based on a BASiCS_Chain object #########################
##########################################################################################

HiddenVarDecomp <- function(Data, object)
{
  
  if(!is(object,"BASiCS_Chain")) stop("'object' is not a BASiCS_Chain class object.")
  
  N = nrow(object@delta); q.bio = ncol(object@delta) 
  UniqueBatch = unique(Data@BatchInfo)
  nBatch = length(UniqueBatch)
    
  if(nBatch > 1) {Theta = apply(object@theta, 1, median)}
  else{ Theta = as.vector(object@theta)}
   
  # To store global values (uses median values across all cells)
  PhiS = apply(object@phi *object@s, 1, median) 
  Aux = (1/(PhiS * object@mu[,1:q.bio])) + object@delta * (Theta+1)
  TechVarGlobal = Theta / ( Aux + Theta )
  BioVarGlobal = (object@delta * (Theta + 1)) / (Aux + Theta)
  
  # To store batch specific values (in arrays)
  TechVarBatch = array(0, dim = c(N, q.bio, nBatch)) # Technical
  BioVarBatch = array(0, dim = c(N, q.bio, nBatch)) # Biological
 
  if(nBatch > 1)
  {
    for(Batch in 1:nBatch)
    {   
      PhiSBatch = apply(object@phi[, Data@BatchInfo == UniqueBatch[Batch]] *
                          object@s[, Data@BatchInfo == UniqueBatch[Batch]], 1, median) 
      Aux = (1/(PhiSBatch * object@mu[,1:q.bio])) + object@delta *(object@theta[,Batch]+1)
      TechVarBatch[,,Batch] = object@theta[,Batch] / ( Aux + object@theta[,Batch] )
      BioVarBatch[,,Batch] = (object@delta * (object@theta[,Batch] + 1)) / (Aux + object@theta[,Batch])
    }     
  }

  if(nBatch > 1)
  {
    list("TechVarGlobal"=TechVarGlobal,
         "BioVarGlobal"=BioVarGlobal,
         "TechVarBatch"=TechVarBatch,
         "BioVarBatch"=BioVarBatch)
  }
  else
  {
    list("TechVarGlobal"=TechVarGlobal,
         "BioVarGlobal"=BioVarGlobal)    
  }
}

HiddenTailProbUp<-function(chain,threshold){return(sum(chain>threshold)/length(chain))}
HiddenTailProbLow<-function(chain,threshold){return(sum(chain<threshold)/length(chain))}

HiddenProbHVG<-function(
  VarThreshold, # Variance contribution threshold for HVG/LVG detection. Value must be between 0 and 1.
  VarDecomp) # Output of the variance decomposition obtained using BASiCS_Variance
{  
  return(apply(VarDecomp$BioVarGlobal, 2, HiddenTailProbUp, threshold = VarThreshold))
}

HiddenProbLVG<-function(
  VarThreshold, # Variance contribution threshold for HVG/LVG detection. Value must be between 0 and 1.
  VarDecomp) # Output of the variance decomposition obtained using BASiCS_Variance
{  
  return(apply(VarDecomp$BioVarGlobal, 2, HiddenTailProbLow, threshold = VarThreshold))
}

HiddenEFDR<-function(
  EviThreshold, # Evidence threshold (it must be contained in (0,1))
  VarThreshold, # Variance contribution threshold choosen by the user (it must be contained in (0,1))
  Prob) # Output of the variance decomposition obtained using BASiCS_Variance
{  
  return(sum((1-Prob)*I(Prob>EviThreshold))/sum(I(Prob>EviThreshold)))
}

HiddenEFNR<-function(
  EviThreshold, # Evidence threshold (it must be contained in (0,1))
  VarThreshold, # Variance contribution threshold choosen by the user (it must be contained in (0,1))
  Prob) # Output of the variance decomposition obtained using BASiCS_Variance  
{
  return(sum(Prob*I(EviThreshold>=Prob))/sum(I(EviThreshold>=Prob)))
}

#' @name BASiCS_VarianceDecomp
#' @aliases BASiCS_VarianceDecomp
#' 
#' @title Decomposition of gene expression variability according to BASiCS
#' 
#' @param Data an object of class \code{\link[BASiCS]{BASiCS_Data-class}}
#' @param object an object of class \code{\link[BASiCS]{BASiCS_Chain-class}}
#' @param GeneNames Vector of characters containing gene names (must be in the same order as in the Counts matrix).
#' @param OrderVariable Ordering variable for output. Must take values in \code{c("GeneNames", "BioVarGlobal", "TechVarGlobal", "ShotNoiseGlobal")}.
#' @param Plot If \code{TRUE}, a barplot of the variance decomposition (global and by batches, if any) is generated
#' @param ... Other arguments to be passed to \code{\link[graphics]{barplot}}
#' 
#' @return A \code{\link[base]{data.frame}} whose first 4 columns correspond to
#' \describe{
#' \item{\code{GeneName}}{Gene name (as indicated by user)}
#' \item{\code{BioVarGlobal}}{Percentage of variance explained by a biological cell-to-cell heterogeneity component (overall across all cells)}
#' \item{\code{TechVarGlobal}}{Percentage of variance explained by the technical cell-to-cell heterogeneity component (overall across all cells)}
#' \item{\code{ShotNoiseGlobal}}{Percentage of variance explained by the shot noise component (baseline, overall across all cells)}
#' }
#' If more than 1 batch of cells are being analysed, the remaining columns contain the corresponding variance decomposition calculated within each batch. 
#'  
#' @examples
#' 
#' # See
#' help(BASiCS_MCMC)
#' 
#' @details See vignette
#' 
#' 
#' @seealso \code{\link[BASiCS]{BASiCS_Chain-class}} 
#' 
#' @author Catalina A. Vallejos \email{catalina.vallejos@@mrc-bsu.cam.ac.uk}
#' 
#' @references Vallejos, Marioni and Richardson (2015). Bayesian Analysis of Single-Cell Sequencing data.
#' 
#' @rdname BASiCS_VarianceDecomp
BASiCS_VarianceDecomp <- function(Data,
                                  object, 
                                  GeneNames = NULL,
                                  OrderVariable = "BioVarGlobal",
                                  Plot = TRUE,
                                  ...)
{  
  if(!(OrderVariable %in% c("GeneNames", "BioVarGlobal", "TechVarGlobal", "ShotNoise"))) stop("Invalid 'OrderVariable' value.")
  
  q.bio = ncol(object@delta) 
  UniqueBatch = unique(Data@BatchInfo)
  nBatch = length(UniqueBatch)
  
  # Calculating variance decomposition
  VarDecomp = HiddenVarDecomp(Data, object)
  
  # Global values 
  BioVarGlobal = apply(VarDecomp$BioVarGlobal, 2, median)
  TechVarGlobal = apply(VarDecomp$TechVarGlobal, 2, median)
  ShotNoiseGlobal = 1-BioVarGlobal-TechVarGlobal
  
  Genes = 1:q.bio
  if(is.null(GeneNames)) {GeneNames = paste("Gene", Genes)}

  if(nBatch > 1)
  {
    VarDecompBatch = NULL
    
    for(Batch in 1:nBatch)
    {
      VarDecompBatch = cbind(VarDecompBatch,
                             apply(VarDecomp$BioVarBatch[,,Batch], 2, median),
                             apply(VarDecomp$TechVarBatch[,,Batch], 2, median),
                             1 - apply(VarDecomp$BioVarBatch[,,Batch], 2, median) - apply(VarDecomp$TechVarBatch[,,Batch], 2, median))
    }
    
    colnames(VarDecompBatch) = paste0(rep(c("BioVarBatch", "TechBatch", "ShotNoiseBatch"),nBatch),rep(1:nBatch,each = 3))

    out = cbind.data.frame("GeneIndex" = Genes,
                           "GeneNames" = GeneNames,
                           "BioVarGlobal" = BioVarGlobal, 
                           "TechVarGlobal" = TechVarGlobal, 
                           "ShotNoiseGlobal" = ShotNoiseGlobal,
                           VarDecompBatch,
                           stringsAsFactors = FALSE)
  }
  else
  {
    out = cbind.data.frame("GeneIndex" = Genes,
                           "GeneNames" = GeneNames,
                           "BioVarGlobal" = BioVarGlobal, 
                           "TechVarGlobal" = TechVarGlobal, 
                           "ShotNoiseGlobal" = ShotNoiseGlobal,
                           stringsAsFactors = FALSE)    
  }
  
  # Re-order before output 
  if(OrderVariable == "GeneNames") orderVar = GeneNames
  if(OrderVariable == "BioVarGlobal") orderVar = BioVarGlobal
  if(OrderVariable == "TechVarGlobal") orderVar = TechVarGlobal
  if(OrderVariable == "ShotNoiseGlobal") orderVar = 1-BioVarGlobal-TechVarGlobal
  out = out[order(orderVar, decreasing = TRUE),]
  
  if(Plot)
  {
    args <- list(...)
    main = ifelse("main" %in% names(args),args$main, "Overall variance decomposition") 
    ylab = ifelse("ylab" %in% names(args),args$ylab, "% of variance")
    beside = ifelse("beside" %in% names(args),args$beside, FALSE)
    if("col" %in% names(args)) {col = args$col} else{col = c("lightblue", "mistyrose", "lightcyan")}
    if("legend" %in% names(args)) {legend = args$legend} else{legend = c("Biological", "Technical", "Baseline")}
    if("args.legend" %in% names(args)) {args.legend = args$args.legend} else{args.legend = list(x = "bottomright")}
    if("names.arg" %in% names(args)) {names.arg = args$names.arg} 
    else
    {
      if(nBatch > 1) {names.arg = c("Overall", paste("Batch ", 1:nBatch))}
      else {names.arg = c("Overall")}
    }
    
    outmat = matrix(apply(out[,-c(1:2)], 2, mean), nrow = 3, byrow = FALSE)
    barplot(outmat, 
            beside = beside, main = main, ylab = ylab, 
            col = col, legend = legend,
            args.legend = args.legend,
            names.arg = names.arg)    
  }
  
  return(out)
}

#' @name BASiCS_DetectHVG
#' @aliases BASiCS_DetectHVG BASiCS_DetectHVG_LVG
#' 
#' @title Detection method for highly and lowly variable genes
#' 
#' @param Data an object of class \code{\link[BASiCS]{BASiCS_Data-class}}
#' @param object an object of class \code{\link[BASiCS]{BASiCS_Chain-class}}
#' @param VarThreshold Variance contribution threshold (must be a positive value, between 0 and 1)
#' @param EviThreshold Optional parameter. Evidence threshold (must be a positive value, between 0 and 1)
#' @param OrderVariable Ordering variable for output. Must take values in \code{c("GeneIndex", "Mu", "Delta", "Sigma", "Prob")}.
#' @param Plot If \code{Plot = T} a plot of the gene specific expression level against HVG or LVG is generated.
#' @param GeneNames Vector of characters containing gene names (must be in the same order as in the Counts matrix).
#' @param ... Graphical parameters (see \code{\link[graphics]{par}}).
#' 
#' @return \code{BASiCS_DetectHVG} returns a list of 4 elements:
#' \describe{
#' \item{\code{Table}}{Matrix whose columns contain}
#'    \describe{
#'    \item{\code{GeneNames}}{Gene names (if names are provided by the user, otherwise returns 'Gene 1', 'Gene 2', etc.)}
#'    \item{\code{Mu}}{Vector of length \code{q.bio}. For each biological gene, posterior median of gene-specific expression levels \eqn{\mu[i]}}
#'    \item{\code{Delta}}{Vector of length \code{q.bio}. For each biological gene, posterior median of gene-specific biological cell-to-cell heterogeneity hyper-parameter \eqn{\delta[i]}}
#'    \item{\code{Sigma}}{Vector of length \code{q.bio}. For each biological gene, proportion of the total variability that is due to a cell-to-cell biological heterogeneity component. }
#'    \item{\code{Prob}}{Vector of length \code{q.bio}. For each biological gene, probability of being highly variable according to the given thresholds.}
#'    \item{\code{HVG}}{Vector of length \code{q.bio}. For each biological gene, indicator of being detected as highly variable according to the given thresholds. }
#'    }
#' \item{\code{EviThreshold}}{Evidence threshold.}
#' \item{\code{EFDR}}{Expected false discovery rate for the given thresholds.}
#' \item{\code{EFNR}}{Expected false negative rate for the given thresholds.}
#' }
#' \code{BASiCS_DetectLVG} produces a similar output, replacing the element \code{HVG} by \code{LVG}, an indicator of a gene being detected as lowly variable according to the given thresholds.  
#'  
#' @examples
#' 
#' # See
#' help(BASiCS_MCMC)
#' 
#' @details See vignette
#' 
#' 
#' @seealso \code{\link[BASiCS]{BASiCS_Chain-class}} 
#' 
#' @author Catalina A. Vallejos \email{catalina.vallejos@@mrc-bsu.cam.ac.uk}
#' 
#' @references Vallejos, Marioni and Richardson (2015). Bayesian Analysis of Single-Cell Sequencing data.
#' 
#' @rdname BASiCS_DetectHVG_LVG
BASiCS_DetectHVG <- function(Data, 
                             object,
                             VarThreshold,
                             EviThreshold = NULL,
                             OrderVariable = "Prob",
                             GeneNames = NULL,
                             Plot = FALSE, 
                             ...)
{
  if(!is(object,"BASiCS_Chain")) stop("'object' is not a BASiCS_Chain class object.")    
  if(VarThreshold<0 | VarThreshold>1 | !is.finite(VarThreshold)) stop("Variance contribution thresholds for HVG/LVG detection must be contained in (0,1)")
  if(!is.logical(Plot) | length(Plot)!=1) stop("Please insert T or F for Plot parameter")
  if(!is.null(EviThreshold))
  {
    if(EviThreshold<0 | EviThreshold>1 | !is.finite(EviThreshold)) 
      stop("Evidence thresholds for HVG and LVG detection must be contained in (0,1) \n For automatic threshold search use EviThreshold = NULL.")    
  }
  if(!(OrderVariable %in% c("GeneNames", "Mu", "Delta", "Sigma", "Prob"))) stop("Invalid 'OrderVariable' value")
  Search = F
  if(is.null(EviThreshold)) Search = T
  
  VarDecomp <- HiddenVarDecomp(Data, object)
  Prob <- HiddenProbHVG(VarThreshold = VarThreshold, VarDecomp = VarDecomp)
  
  if(length(EviThreshold) == 0)
  {
    EviThresholds <- seq(0.5,0.9995,by=0.0005)

    EFDR <- sapply(EviThresholds, HiddenEFDR, VarThreshold = VarThreshold, Prob = Prob)
    EFNR <- sapply(EviThresholds, HiddenEFNR, VarThreshold = VarThreshold, Prob = Prob)
    
    above<-EFDR>EFNR
    optimal<-which(diff(above)!=0)
    EviThreshold = EviThresholds[optimal]
    if(length(optimal)>0){OptThreshold <- c(EviThreshold, EFDR[optimal], EFNR[optimal])}
    else 
    {
      print("It is not possible to find an optimal evidence threshold for the given variance contribution threshold. \n")
      optimal <- round(median(which(abs(EFDR - EFNR) == min(abs(EFDR - EFNR), na.rm = T))))
      if(length(optimal)>0)
      {
        print("Returned value is such that the difference between EFDR and EFNR is minimised.")
        OptThreshold <- c(EviThreshold, EFDR[optimal], EFNR[optimal])
      }
      else
      {
        cat("Numerical issues when computing EFDR and EFNR. Please try a different variance contribution threshold")
        OptThreshold <- rep("Not found",3) 
      }
    }  
  }
  else
  {
     EFDR = HiddenEFDR(EviThreshold, VarThreshold, Prob)
     EFNR = HiddenEFNR(EviThreshold, VarThreshold, Prob)
     OptThreshold <- c(EviThreshold, EFDR, EFNR)
  } 
  
  Sigma <- apply(VarDecomp$BioVarGlobal, 2, median)
  Mu <- apply(object@mu[,1:length(Sigma)], 2, median)
  Delta <- apply(object@delta, 2, median)
  if(OptThreshold[1] == "Not found") {HVG = rep("Not found", length(Sigma))}
  else
  { 
    HVG <- ifelse(Prob > OptThreshold[1], TRUE, FALSE);  
    HVG <- ifelse(Prob >= 0.5, HVG, FALSE); 
  }
    
  qbio = length(Sigma)
  Genes = 1:qbio
  if(is.null(GeneNames)) {GeneNames = paste("Gene", Genes)}
  
  if(Plot)
  {    
    args <- list(...)
    
    if(Search)
    {      
      par(ask=T)
      
      plot(EviThresholds, EFDR, type = "l", lty = 1, bty = "n", ylab = "Error rate", xlab = "Evidence threshold", ylim = c(0,1))
      lines(EviThresholds, EFNR, lty = 2)      
      legend('topleft', c("EFDR", "EFNR"), lty = 1:2, bty = "n")
    }
    
    if("ylim" %in% names(args)) {ylim = args$ylim} else{ylim = c(0, 1)}
    if("xlim" %in% names(args)) {xlim = args$xlim} else{xlim = c(min(Mu),max(Mu))}
    cex = ifelse("cex" %in% names(args),args$cex, 1.5)
    pch = ifelse("pch" %in% names(args),args$pch, 16)
    col = ifelse("col" %in% names(args),args$col, 8)
    bty = ifelse("bty" %in% names(args),args$bty, "n")
    cex.lab = ifelse("cex.lab" %in% names(args),args$cex.lab, 1)
    cex.axis = ifelse("cex.axis" %in% names(args),args$cex.axis, 1)
    cex.main = ifelse("cex.main" %in% names(args),args$cex.main, 1) 
    xlab = ifelse("xlab" %in% names(args),args$xlab, expression(mu[i]))
    ylab = ifelse("ylab" %in% names(args),args$ylab, "HVG probability")
    main = ifelse("main" %in% names(args),args$main, "") 
    
    plot(Mu, Prob, log="x", pch = pch, ylim = ylim, xlim = xlim, col = col, cex = cex,
                    bty = bty, cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main,
                    xlab = xlab, ylab = ylab, main = main)
    abline(h = OptThreshold[1], lty = 2, col = "black")
    points(Mu[HVG], Prob[HVG], pch = pch, col = "red", cex = cex)
    
    par(ask=F)
  }
  
  cat(paste(sum(HVG), " genes classified as highly variable using: \n"))
  cat(paste("- Variance contribution threshold = ", round(100*VarThreshold,2), "% \n"))    
  cat(paste("- Evidence threshold = ", OptThreshold[1], "\n"))
  cat(paste("- EFDR = ", round(100*OptThreshold[2],2), "% \n"))  
  cat(paste("- EFNR = ", round(100*OptThreshold[3],2), "% \n"))  
  
  GeneIndex = 1:length(Mu)
  Table = cbind.data.frame("GeneIndex" = Genes,
                           "GeneNames" = GeneNames,
                           "mu" = Mu,
                           "delta" = Delta,
                           "Sigma" = Sigma,
                           "Prob" = Prob,
                           "HVG" = HVG, stringsAsFactors = FALSE)
  
  if(OrderVariable == "GeneNames") orderVar = GeneNames
  if(OrderVariable == "Mu") orderVar = Mu
  if(OrderVariable == "Delta") orderVar = Delta
  if(OrderVariable == "Sigma") orderVar = Sigma
  if(OrderVariable == "Prob") orderVar = Prob
  Table = Table[order(orderVar, decreasing = TRUE),]
  
  list("Table" = Table, 
       "EviThreshold" = OptThreshold[1], "EFDR" = OptThreshold[2], "EFNR" = OptThreshold[3])
}

#' @name BASiCS_DetectLVG
#' @aliases BASiCS_DetectLVG BASiCS_DetectHVG_LVG
#' @rdname BASiCS_DetectHVG_LVG
BASiCS_DetectLVG <- function(Data, 
                             object,
                             VarThreshold,
                             EviThreshold = NULL,
                             OrderVariable = "Prob",
                             GeneNames = NULL,
                             Plot = FALSE, 
                             ...)
{
  if(!is(object,"BASiCS_Chain")) stop("'object' is not a BASiCS_Chain class object.")    
  if(VarThreshold<0 | VarThreshold>1 | !is.finite(VarThreshold)) stop("Variance contribution thresholds for HVG/LVG detection must be contained in (0,1)")
  if(!is.logical(Plot) | length(Plot)!=1) stop("Please insert T or F for Plot parameter")
  if(!is.null(EviThreshold))
  {
    if(EviThreshold<0 | EviThreshold>1 | !is.finite(EviThreshold)) 
      stop("Evidence thresholds for HVG and LVG detection must be contained in (0,1) \n For automatic threshold search use EviThreshold = NULL.")    
  }
  if(!(OrderVariable %in% c("GeneNames", "GeneIndex", "Mu", "Delta", "Sigma", "Prob"))) stop("Invalid 'OrderVariable' value")
  Search = F
  if(is.null(EviThreshold)) Search = T
  
  VarDecomp <- HiddenVarDecomp(Data, object)
  Prob <- HiddenProbLVG(VarThreshold = VarThreshold, VarDecomp = VarDecomp)
  
  if(length(EviThreshold) == 0)
  {
    EviThresholds <- seq(0.5,0.9995,by=0.0005)
    EFDR <- sapply(EviThresholds, HiddenEFDR, VarThreshold = VarThreshold, Prob = Prob)
    EFNR <- sapply(EviThresholds, HiddenEFNR, VarThreshold = VarThreshold, Prob = Prob)
    
    above<-EFDR>EFNR
    optimal<-which(diff(above)!=0)
    EviThreshold = EviThresholds[optimal]
    if(length(optimal)>0){OptThreshold <- c(EviThreshold, EFDR[optimal], EFNR[optimal])}
    else 
    {
      print("It is not possible to find an optimal evidence threshold for the given variance contribution threshold. \n")
      optimal <- round(median(which(abs(EFDR - EFNR) == min(abs(EFDR - EFNR), na.rm = T))))
      if(length(optimal)>0)
      {
        print("Returned value is such that the difference between EFDR and EFNR is minimised.")
        OptThreshold <- c(EviThreshold, EFDR[optimal], EFNR[optimal])
      }
      else
      {
        cat("Numerical issues when computing EFDR and EFNR. Please try a different variance contribution threshold")
        OptThreshold <- rep("Not found",3) 
      }
    } 
  }
  else
  {
    EFDR = HiddenEFDR(EviThreshold, VarThreshold, Prob)
    EFNR = HiddenEFNR(EviThreshold, VarThreshold, Prob)
    OptThreshold <- c(EviThreshold, EFDR, EFNR)
  }
  
  Sigma <- apply(VarDecomp$BioVarGlobal, 2, median)
  Mu <- apply(object@mu[,1:length(Sigma)], 2, median)
  Delta <- apply(object@delta, 2, median)
  if(OptThreshold[1] == "Not found") {LVG = rep("Not found", length(Sigma))}
  else
  { 
    LVG <- ifelse(Prob > OptThreshold[1], TRUE, FALSE) 
    LVG <- ifelse(Prob >= 0.5, LVG, FALSE); 
  }
  
  qbio = length(Sigma)
  Genes = 1:qbio
  if(is.null(GeneNames)) {GeneNames = paste("Gene", Genes)}
  
  if(Plot)
  {    
    args <- list(...)
        
    if(Search)
    {
      par(ask=T)
      
      plot(EviThresholds, EFDR, type = "l", lty = 1, bty = "n", ylab = "Error rate", xlab = "Evidence threshold", ylim = c(0,1))
      lines(EviThresholds, EFNR, lty = 2)      
      legend('topleft', c("EFDR", "EFNR"), lty = 1:2, bty = "n")
    }
    
    if("ylim" %in% names(args)) {ylim = args$ylim} else{ylim = c(0, 1)}
    if("xlim" %in% names(args)) {xlim = args$xlim} else{xlim = c(min(Mu),max(Mu))}
    cex = ifelse("cex" %in% names(args),args$cex, 1.5)
    pch = ifelse("pch" %in% names(args),args$pch, 16)
    col = ifelse("col" %in% names(args),args$col, 8)
    bty = ifelse("bty" %in% names(args),args$bty, "n")
    cex.lab = ifelse("cex.lab" %in% names(args),args$cex.lab, 1)
    cex.axis = ifelse("cex.axis" %in% names(args),args$cex.axis, 1)
    cex.main = ifelse("cex.main" %in% names(args),args$cex.main, 1) 
    xlab = ifelse("xlab" %in% names(args),args$xlab, expression(mu[i]))
    ylab = ifelse("ylab" %in% names(args),args$ylab, "LVG probability")
    main = ifelse("main" %in% names(args),args$main, "") 
    
    plot(Mu, Prob, log="x", pch = pch, ylim = ylim, xlim = xlim, col = col, cex = cex, 
         bty = bty, cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main,
         xlab = xlab, ylab = ylab, main = main)
    abline(h = OptThreshold[1], lty = 2, col = "black")
    points(Mu[LVG], Prob[LVG], pch = pch, col = "red", cex = cex)
    
    par(ask=F)
    
  }
  
  cat(paste(sum(LVG), " genes classified as lowly variable using: \n"))
  cat(paste("- Variance contribution threshold = ", round(100*VarThreshold,2), "% \n"))    
  cat(paste("- Evidence threshold = ", OptThreshold[1], "\n"))
  cat(paste("- EFDR = ", round(100*OptThreshold[2],2), "% \n"))  
  cat(paste("- EFNR = ", round(100*OptThreshold[3],2), "% \n")) 
   
  GeneIndex = 1:length(Mu)
  Table = cbind.data.frame("GeneIndex" = Genes,
                           "GeneNames" = GeneNames,
                           "mu" = Mu,
                           "delta" = Delta,
                           "Sigma" = Sigma,
                           "Prob" = Prob,
                           "LVG" = LVG, stringsAsFactors = FALSE)
  
  if(OrderVariable == "GeneNames") orderVar = GeneNames
  if(OrderVariable == "Mu") orderVar = Mu
  if(OrderVariable == "Delta") orderVar = Delta
  if(OrderVariable == "Sigma") orderVar = Sigma
  if(OrderVariable == "Prob") orderVar = Prob
  Table = Table[order(orderVar, decreasing = TRUE),]
  
  list("Table" = Table, 
       "EviThreshold" = OptThreshold[1], "EFDR" = OptThreshold[2], "EFNR" = OptThreshold[3])
}


#' @name BASiCS_VarThresholdSearchHVG
#' @aliases BASiCS_VarThresholdSearchHVG BASiCS_VarThresholdSearchHVG_LVG
#' 
#' @title Detection method for highly and lowly variable genes using a grid of variance contribution thresholds
#' 
#' @param Data an object of class \code{\link[BASiCS]{BASiCS_Data-class}}
#' @param object an object of class \code{\link[BASiCS]{BASiCS_Chain-class}}
#' @param VarThresholdsGrid Grid of values for the variance contribution threshold (they must be contained in (0,1))
#' @param PrintProgress If \code{PrintProgress = TRUE}, partial output is printed in the console.
#' 
#' @examples
#' 
#' # See
#' help(BASiCS_MCMC)
#' 
#' @details See vignette
#' 
#' @return 
#' \describe{
#' \item{\code{BASiCS_VarThresholdSearchHVG}}{A table displaying the results of highly variable genes detecting for different variance contribution thresholds.}
#' \item{\code{BASiCS_VarThresholdSearchLVG}}{A table displaying the results of lowly variable genes detecting for different variance contribution thresholds.oo}
#' }
#' 
#'    
#' @seealso \code{\link[BASiCS]{BASiCS_Chain-class}} 
#' 
#' @author Catalina A. Vallejos \email{catalina.vallejos@@mrc-bsu.cam.ac.uk}
#' 
#' @rdname BASiCS_VarThresholdSearchHVG_LVG
BASiCS_VarThresholdSearchHVG=function(
  Data,
  object, 
  VarThresholdsGrid, # 
  PrintProgress = FALSE)
{
  
  if(!is(object,"BASiCS_Chain")) stop("'object' is not a BASiCS_Chain class object.")    
  if(sum(VarThresholdsGrid<0)>0 | sum(VarThresholdsGrid>1)>0 | sum(!is.finite(VarThresholdsGrid))>0 )
    stop("Variance contribution thresholds for HVG and LVG detection must be contained in (0,1).")    
  
  Table=matrix(0,nrow=length(VarThresholdsGrid),ncol=5)
  colnames(Table)=c("Var. Threshold (%)","EFDR (%)", "EFNR (%)","Optimal evidence thres.","# Detected genes")
  
  for(i in 1:length(VarThresholdsGrid))
  {
    VarThreshold=VarThresholdsGrid[i]
    
    if(PrintProgress) {print(paste0("Evaluating variance contribution threshold = ",100*VarThreshold," % ...")); cat("\n")}
    
    DetectHVG <- BASiCS_DetectHVG(Data, object, VarThreshold = VarThreshold)
    
    Table[i,]=c(100*VarThreshold, round(100*DetectHVG$EFDR,2),
                round(100*DetectHVG$EFNR,2), 
                DetectHVG$EviThreshold,sum(as.numeric(DetectHVG$Table[,7])))
  }
  return(Table)
}

#' @name BASiCS_VarThresholdSearchLVG
#' @aliases BASiCS_VarThresholdSearchLVG BASiCS_VarThresholdSearchHVG_LVG
#' @rdname BASiCS_VarThresholdSearchHVG_LVG
BASiCS_VarThresholdSearchLVG=function(
  Data,
  object, 
  VarThresholdsGrid, # Range of values for the variance contribution threshold (they must be contained in (0,1))
  PrintProgress = FALSE)
{
  
  if(!is(object,"BASiCS_Chain")) stop("'object' is not a BASiCS_Chain class object.")    
  if(sum(VarThresholdsGrid<0)>0 | sum(VarThresholdsGrid>1)>0 | sum(!is.finite(VarThresholdsGrid))>0 )
    stop("Variance contribution thresholds for HVG and LVG detection must be contained in (0,1).")    
  
  Table=matrix(0,nrow=length(VarThresholdsGrid),ncol=5)
  colnames(Table)=c("Var. Threshold (%)","EFDR (%)", "EFNR (%)","Optimal evidence thres.","# Detected genes")
  
  for(i in 1:length(VarThresholdsGrid))
  {
    VarThreshold=VarThresholdsGrid[i]
    
    if(PrintProgress) {print(paste0("Evaluating variance contribution threshold = ",100*VarThreshold," % ...")); cat("\n")}
    
    DetectLVG <- BASiCS_DetectLVG(Data, object, VarThreshold = VarThreshold)
    
    Table[i,]=c(100*VarThreshold, round(100*DetectLVG$EFDR,2),round(100*DetectLVG$EFNR,2), 
                DetectLVG$EviThreshold,sum(as.numeric(DetectLVG$Table[,7])))
    
  }
  return(Table)
}

#' @name DenoisedRates
#' 
#' @title Calculates normalised and denoised expression rates 
#' 
#' @description Calculates normalised and denoised expression rates, by removing the effect of technical variation. 
#' 
#' @param Data an object of class \code{\link[BASiCS]{BASiCS_Data-class}}
#' @param Chain an object of class \code{\link[BASiCS]{BASiCS_Chain-class}}
#' @param PrintProgress If \code{T}, partial progress information is printed in the console. 
#' @param Propensities If \code{T}, returns underlying expression propensitites \eqn{\rho[ij]}. Otherwise, denoised rates \eqn{\mu[i] \rho[ij]} are returned.
#' 
#' @examples
#' 
#' # See
#' help(BASiCS_MCMC)
#' 
#' @details See vignette
#' 
#' @return A matrix of normalised and denoised expression counts. 
#'    
#' @seealso \code{\link[BASiCS]{BASiCS_Data-class}} \code{\link[BASiCS]{BASiCS_Chain-class}} 
#' 
#' @author Catalina A. Vallejos \email{catalina.vallejos@@mrc-bsu.cam.ac.uk}
#' 
#' @rdname DenoisedRates
BASiCS_DenoisedRates=function(
  Data, 
  Chain,
  PrintProgress = F,
  Propensities = T)
{
  if(!is(Data,"BASiCS_Data")) stop("'Data' is not a BASiCS_Data class object.") 
  if(!is(Chain,"BASiCS_Chain")) stop("'Chain' is not a BASiCS_Chain class object.") 
  
  N=dim(Chain@delta)[1]; q.bio=dim(Chain@delta)[2]; n=dim(Chain@phi)[2]
  
  print(paste("This calculation requires a loop across the",N, "MCMC iterations")); 
  print("You might need to be patient ... "); cat("\n")
  
  rho=matrix(0,ncol=n,nrow=q.bio)
  for(m in 1:N)
  {
    if(PrintProgress) {print(paste("Iteration",m,"has been completed."))}
    aux1=Data@Counts[1:q.bio,]+1/Chain@delta[m,]
    aux2=t(tcrossprod(Chain@phi[m,]*Chain@nu[m,],Chain@mu[m,1:q.bio]))+1/Chain@delta[m,]
    rho=rho+aux1/aux2
  }
  rho=rho/N
  
  if(Propensities) {out = rho}
  else {out = rho * apply(Chain@mu[,1:q.bio],2,median)}
  
  return(out)
  
}

#' @name DenoisedCounts
#' 
#' @title Calculates normalised and denoised expression expression counts 
#' 
#' @description Calculates normalised and denoised expression counts, by removing the effect of technical variation. 
#' 
#' @param Data an object of class \code{\link[BASiCS]{BASiCS_Data-class}}
#' @param Chain an object of class \code{\link[BASiCS]{BASiCS_Chain-class}}
#' 
#' @examples
#' 
#' # See
#' help(BASiCS_MCMC)
#' 
#' @details See vignette
#' 
#' @return A matrix of normalised and denoised expression counts. 
#'    
#' @seealso \code{\link[BASiCS]{BASiCS_Data-class}} \code{\link[BASiCS]{BASiCS_Chain-class}} 
#' 
#' @author Catalina A. Vallejos \email{catalina.vallejos@@mrc-bsu.cam.ac.uk}
#' 
#' @rdname DenoisedCounts
BASiCS_DenoisedCounts=function(
  Data, 
  Chain)
{
  if(!is(Data,"BASiCS_Data")) stop("'Data' is not a BASiCS_Data class object.") 
  if(!is(Chain,"BASiCS_Chain")) stop("'Chain' is not a BASiCS_Chain class object.") 
  
  N=dim(Chain@delta)[1]; q.bio=dim(Chain@delta)[2]; n=dim(Chain@phi)[2]
  
  Phi = apply(Chain@phi,2,median)
  Nu = apply(Chain@nu,2,median)
  
  out1 = t( t(Data@Counts[!Data@Tech,]) / (Phi * Nu) )
  out2 = t( t(Data@Counts[Data@Tech,]) / Nu )
  out = rbind(out1, out2)
 
  return(out)  
}