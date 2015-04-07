#' Simulates expression counts according to the model implemented in BASiCS 
#' 
#' \code{BASiCS_Sim} creates a simulated dataset from the model implemented in BASiCS. 
#' This function is used in order to illustrate the performance of the \code{BASiCS} library.
#' 
#' @param mu Gene-specific expression levels, defined as true input molecules in case of technical genes 
#' (vector of length \code{q}, technical genes located at the end of the vector, all elements must be positive numbers)
#' @param delta Gene-specific biological cell-to-cell heterogeneity hyper-parameters, biological genes only 
#' (vector of length \code{q.bio}, all elements must be positive numbers)
#' @param phi Cell-specific mRNA content normalising constants 
#' (vector of length \code{n}, all elements must be positive numbers and the sum of its elements must be equal to \code{n})
#' @param s Cell-specific capture efficiency (or amplification biases if not using UMI based counts) normalising constants
#' (vector of length \code{n}, all elements must be positive numbers)
#' @param theta Technical variability hyper-parameter (must be positive)
#' 
#' @return \code{BASiCS_Sim} returns a list with 3 elements:
#' 
#'  \code{Counts} Matrix of dimensions \code{q} times \code{n} whose elements corresponds to the simulated expression counts. 
#' First \code{q.bio} rows correspond to biological genes. Last \code{q-q.bio} rows correspond to technical spike-in genes. 
#' 
#'  \code{Tech} Vector of length \code{q}. If \code{F} the gene is biological; otherwise the gene is spike-in.
#'  
#'  \code{SpikeInput} Vector of length \code{q-q.bio} whose elements indicate the simulated input concentrations for the spike-in genes. 
#' 
#' @examples
#' # Loading estimated parameters for Mouse ESC dataset analysed in Vallejos et al (2015). 
#' data(ParamMouse)
#' attach(ParamMouse)
#' DataSimLong = BASiCS_Sim(mu = Mu, delta = Delta, phi = Phi, s = S, theta = Theta)
#' names(DataSimLong)
#' dim(DataSimLong$Counts)
#' table(DataSimLong$Tech)
#' length(DataSimLong$SpikeInput)
#' detach(ParamMouse)
#' 
#' # Simulating a smaller dataset that can be used for illustration purposes
#' # Include 100 biological genes and 46 spike-in genes in 41 cells
#' data(ParamMouseShort)
#' attach(ParamMouseShort)
#' DataSimShort = BASiCS_Sim(mu = MuShort, delta = DeltaShort, phi = PhiShort, s = SShort, theta = ThetaShort)
#' names(DataSimShort)
#' dim(DataSimShort$Counts)
#' table(DataSimShort$Tech)
#' length(DataSimShort$SpikeInput)
#' detach(ParamMouseShort)
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
  # Arguments checking
  testit::assert("Error in 'BASiCS_Sim': invalid argument type (mu)", is.vector(mu) & is.numeric(mu))
  testit::assert("Error in 'BASiCS_Sim': invalid argument type (delta)", is.vector(delta) & is.numeric(delta))
  testit::assert("Error in 'BASiCS_Sim': invalid argument type (phi)", is.vector(phi) & is.numeric(phi))
  testit::assert("Error in 'BASiCS_Sim': invalid argument type (s)", is.vector(s) & is.numeric(s))
  testit::assert("Error in 'BASiCS_Sim': invalid argument type (theta)", is.numeric(theta) & length(theta)==1)
  
  # Number of cells 
  n = length(phi)
  # Total number of genes, including biological and technical ones
  q = length(mu)
  # Number of biological genes
  q.bio = length(delta)
  
  testit::assert("Error in 'BASiCS_SIM': invalid number of cells and/or genes",all(c(n,q,q.bio,q-q.bio)>0))
  
  testit::assert("Error in 'BASiCS_SIM': invalid specification of parameter values", all(mu>0) & all(delta>=0) & 
          all(phi>0) & sum(phi)==n & all(s>0) & theta>=0)
  
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

  list("Counts" = Counts.sim, "Tech" = ifelse(1:q > q.bio, T, F), "SpikeInput" = mu[(q.bio+1):q])
}

#' Starting values for BASiCS_MCMC
#' 
#' Auxiliary function that creates starting values for BASiCS_MCMC. 
#' 
#' @details Typically, this function will not be directly called. It is internally used by \code{\link[BASiCS]{BASiCS_MCMC}}.
#' 
#' @param Data List containing 3 elements: 
#' \describe{
#'   \item{\code{Counts}}{Matrix of raw expression counts. Dimensions: \code{n} columns (cells) and \code{q} rows (genes). Spike-in genes located at the last \code{q-q.bio} rows.}
#'   \item{\code{Tech}}{Spike-in genes indicator. If \code{F} the gene is biological; otherwise the gene is spike-in. Length: \code{q} (genes).}
#'   \item{\code{SpikeInput}}{Input number of mRNA molecules of each spike-in gene that are added to each cell (in the same order as they appear in \code{Counts})}
#' }
#' @param ... Optional arguments (lower bounds for starting values of adaptive proposal variances, log-scale):
#' \describe{
#'  \item{\code{ls.mu0}}{Related to gene-specific expression levels}
#'  \item{\code{ls.delta0}}{Related to gene-specific biological cell-to-cell heterogeneity hyper-parameters}
#'  \item{\code{ls.kappa0}}{Related to cell-specific mRNA content normalising constants (logit-scale)}
#'  \item{\code{ls.nu0}}{Related to cell-specific random effects (technical noise)}
#'  \item{\code{ls.theta0}}{Related to technical variability hyper-parameter}
#' }
#' 
#' @return A list containing starting values for all model parameters and the corresponding variances of the adaptive proposals (log-scale).
#' @examples
#' # Simulating a smaller dataset that can be used for illustration purposes
#' # Include 100 biological genes and 46 spike-in genes in 41 cells
#' data(ParamMouseShort)
#' attach(ParamMouseShort)
#' DataSimShort = BASiCS_Sim(mu = MuShort, delta = DeltaShort, phi = PhiShort, s = SShort, theta = ThetaShort)
#' detach(ParamMouseShort)
#' 
#' # Creating starting values for BASiCS_MCMC
#' MCMC_Start <- BASiCS_MCMC_Start(Data = DataSimShort)
#' 
#' # Creating starting values for BASiCS_MCMC 
#' # With user-defined lower bounds for adaptive proposal variances
#' MCMC_Start2 <- BASiCS_MCMC_Start(Data = DataSimShort, ls.mu0 = -6)
#' 
#' @author Catalina A. Vallejos (catalina.vallejos 'at' mrc-bsu.cam.ac.uk)
#' @references Vallejos, Marioni and Richardson (2015). Bayesian Analysis of Single-Cell Sequencing data.
BASiCS_MCMC_Start<-function(
  Data,
  ...)
{
  
  testit::assert("Error in 'BASiCS_SIM': invalid argument type (Counts)", is.matrix(Data$Counts))
  testit::assert("Error in 'BASiCS_SIM': invalid argument type (Tech)", is.vector(Data$Tech))
  testit::assert("Error in 'BASiCS_SIM': invalid argument type (SpikeInput)", is.vector(Data$SpikeInput))

  testit::assert("Error in 'BASiCS_SIM': invalid dimension of the arguments", is.vector(Data$SpikeInput))
 
  # Number of instrinsic genes
  q.bio<-sum(Data$Tech==F)
  # Number of cells
  n <- dim(Data$Counts)[2]
  
  # Initialize phi as a vector of ones
  kappa0 = rep(0, times = n)
  # Initialize s as the empirical capture efficiency rates
  s0=colSums(Data$Counts[Data$Tech,])/sum(Data$SpikeInput); nu0=s0 
  # Initialize mu using average 'normalised counts' across cells 
  # (correcting by empirical capture efficiency rates)
  # and true input values for spike-in genes
  nCountsBio <- t( t(Data$Counts[!Data$Tech,]) / s0 )
  meansBio <- rowMeans( nCountsBio )
  mu0<-c(meansBio,Data$SpikeInput)
  
  # Random stating values for delta and theta
  delta0=rgamma(q.bio,1,1)
  theta0=rgamma(1,1,1)
  
  # If given, load default values for adaptive proposal variances
  args <- list(...)
  ls.mu0 = ifelse("ls.mu0" %in% names(args),args$ls.mu0,-4)
  ls.delta0 = ifelse("ls.delta0" %in% names(args),args$ls.delta0,-3)
  ls.kappa0 = ifelse("ls.kappa0" %in% names(args),args$ls.kappa0,-8)
  ls.nu0 = ifelse("ls.nu0" %in% names(args),args$ls.nu0,-10)
  ls.theta0 = ifelse("ls.theta0" %in% names(args),args$ls.theta0,-6)
  
  # Starting values for the proposal variances
  ls.mu0 =  pmax(2 * log (0.02 * abs(log(mu0))),ls.mu0)
  ls.delta0 =  pmax(2 * log (0.02 * abs(log(delta0))),ls.delta0)
  ls.kappa0 =  pmax(2 * log (0.02 * abs(kappa0)),ls.kappa0)
  ls.nu0 =  pmax(2 * log (0.02 * abs(log(nu0))),ls.nu0)
  ls.theta0 =  pmax(2 * log (0.02 * abs(log(theta0))),ls.theta0)
  
  return(list("mu0"=mu0,"delta0"=delta0,"kappa0"=kappa0,"s0"=s0,"nu0"=nu0,"theta0"=theta0,
              "ls.mu0"=ls.mu0,"ls.delta0"=ls.delta0,"ls.kappa0"=ls.kappa0,"ls.nu0"=ls.nu0,"ls.theta0"=ls.theta0))
}

#' BASiCS_MCMC
BASiCS_MCMC<- function( ){}
