#' @title Creates a SummarizedExperiment object from a matrix of expression counts and experimental information about spike-in genes
#'
#' @description \code{newBASiCS_Data} creates a \code{\linkS4class{SummarizedExperiment}} object from
#' a matrix of expression counts and experimental information about spike-in genes.
#'
#' @param Counts Matrix of dimensions \code{q} times \code{n} whose elements contain the expression counts to be analyses
#' (including biological and technical spike-in genes). Gene names must be stored as `rownames(Counts)`.
#' @param Tech Logical vector of length \code{q}. If \code{Tech = FALSE} the gene is biological; otherwise the gene is spike-in.
#' @param SpikeInfo \code{data.frame} whose first and second columns contain the gene names assigned to the spike-in genes
#' (they must match the ones in `rownames(Counts)`) and the associated input number of molecules, respectively.
#' @param BatchInfo Vector of length \code{n} whose elements indicate batch information.
#'
#' @return An object of class \code{\link[BASiCS]{BASiCS_Data-class}}.
#'
#' @examples
#'
#'
#' # Expression counts
#' set.seed(1)
#' Counts = Counts = matrix(rpois(50*10, 2), ncol = 10)
#' rownames(Counts) <- c(paste0("Gene", 1:40), paste0("Spike", 1:10))
#'
#' # Technical information
#' Tech = c(rep(FALSE,40),rep(TRUE,10))
#'
#' # Spikes input number of molecules
#' set.seed(2)
#' SpikeInfo <- data.frame(gene=rownames(Counts)[Tech],amount=rgamma(10,1,1))
#'
#' # Creating a BASiCS_Data object (no batch effect)
#' DataExample = newBASiCS_Data(Counts, Tech, SpikeInfo)
#'
#' # Creating a BASiCS_Data object (with batch effect)
#' BatchInfo = c(rep(1, 5), rep(2, 5))
#' DataExample = newBASiCS_Data(Counts, Tech, SpikeInfo, BatchInfo)
#'
#' # Thanks to Simon Andrews for reporting an issue in previous version of this documentation
#'
#' @seealso \code{\link[BASiCS]{BASiCS_Data-class}}
#'
#' @author Catalina A. Vallejos \email{cnvallej@@uc.cl}
#'
#' @references Vallejos, Marioni and Richardson (2015). Bayesian Analysis of Single-Cell Sequencing data.
newBASiCS_Data <- function(Counts, Tech, SpikeInfo, BatchInfo = NULL)
{
  
  if(is.null(BatchInfo)) {BatchInfo = rep(1, times = ncol(Counts))}
  if(is.factor(BatchInfo)) {
    BIunused <- length(levels(BatchInfo)) - length(unique(BatchInfo))
    if(BIunused > 0){
      message(sprintf(paste("'BatchInfo' was supplied as a 'factor'. All levels of 'BatchInfo' should be represented among cells. Otherwise, 'BASiCS_MCMC' will fail to store MCMC correctly. Therefore, ", BIunused, " unused batch levels have been removed from 'BatchInfo'. See 'help('droplevels')' for more information.")))
      BatchInfo <- droplevels(BatchInfo)
    }
  }
  
  # Re-ordering genes
  Counts = rbind(Counts[!Tech,], Counts[Tech,])
  Tech =c(Tech[!Tech], Tech[Tech])
  GeneNames <- rownames(Counts)
  
  if(!is.null(SpikeInfo))
  {
    # Extracting spike-in input molecules in the correct order
    if(sum(!(GeneNames[Tech] %in% SpikeInfo[,1])) > 0) stop("SpikeInfo is missing information for some of the spikes")
    if(sum(!(SpikeInfo[,1] %in% GeneNames[Tech])) > 0) stop("SpikeInfo includes spikes that are not in the Counts matrix")
    matching <- match(GeneNames[Tech], SpikeInfo[,1])
    SpikeInput <- SpikeInfo[matching,2]    
  }
  else {SpikeInput = 1}
  
  Data <- new("BASiCS_Data", Counts = Counts, Tech = Tech, SpikeInput = SpikeInput,
              GeneNames = GeneNames, BatchInfo = BatchInfo)
  show(Data)
  cat('\n')
  cat('NOTICE: BASiCS requires a pre-filtered dataset \n')
  cat('    - You must remove poor quality cells before creating the BASiCS data object \n')
  cat('    - We recommend to pre-filter very lowly expressed transcripts before creating the object. \n')
  cat('      Inclusion criteria may vary for each data. For example, remove transcripts \n')
  cat('          - with very low total counts across of all of the samples \n')
  cat('          - that are only expressed in a few cells \n')
  cat('            (by default genes expressed in only 1 cell are not accepted) \n')
  cat('          - with very low total counts across the samples where the transcript is expressed \n')
  cat('\n')
  cat(' BASiCS_Filter can be used for this purpose \n')
  return(Data)
}