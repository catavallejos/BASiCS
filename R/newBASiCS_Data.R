#' @title Creates a SingleCellExperiment object from a matrix of expression counts and experimental information about spike-in genes
#'
#' @description \code{newBASiCS_Data} creates a \code{\linkS4class{SingleCellExperiment}} object from
#' a matrix of expression counts and experimental information about spike-in genes.
#'
#' @param Counts Matrix of dimensions \code{q} times \code{n} whose elements contain the expression counts to be analyses
#' (including biological and technical spike-in genes). Gene names must be stored as \code{rownames(Counts)}.
#' @param Tech Logical vector of length \code{q}. If \code{Tech = FALSE} the gene is biological; otherwise the gene is spike-in.
#' @param SpikeInfo \code{data.frame} whose first and second columns contain the gene names assigned to the spike-in genes
#' (they must match the ones in \code{rownames(Counts)}) and the associated input number of molecules, respectively.
#' @param BatchInfo Vector of length \code{n} whose elements indicate batch information. Not required if a single batch is present on the data. 
#' Default value: \code{BatchInfo = NULL}. 
#'
#' @return An object of class \code{\linkS4class{SingleCellExperiment}}.
#'
#' @examples
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
#'
#' # Thanks to Simon Andrews for reporting an issue in previous version of this documentation
#'
#' @seealso \code{\linkS4class{SingleCellExperiment}}
#'
#' @author Catalina A. Vallejos \email{cnvallej@@uc.cl} and Nils Eling
#'
#' @references 
#' Vallejos, Marioni and Richardson (2015). Bayesian Analysis of Single-Cell Sequencing data. PLoS Computational Biology. 
#' 
newBASiCS_Data <- function(Counts, Tech, SpikeInfo, BatchInfo = NULL)
{
  # Validity checks for SpikeInfo
  if(!is.data.frame(SpikeInfo)) stop("'SpikeInfo' must be a 'data.frame'")
  if(data.table::is.data.table(SpikeInfo)) stop("'SpikeInfo' must be a 'data.frame'")
  
  if(is.null(BatchInfo)) {BatchInfo = rep(1, times = ncol(Counts))}
  
  # Re-ordering genes
  Counts = as.matrix(rbind(Counts[!Tech,], Counts[Tech,]))
  Tech = c(Tech[!Tech], Tech[Tech])
  GeneNames <- rownames(Counts)
  
  if(!is.null(SpikeInfo))
  {
    # Extracting spike-in input molecules in the correct order
    if(sum(!(GeneNames[Tech] %in% SpikeInfo[,1])) > 0) stop("'SpikeInfo' is missing information for some of the spikes")
    if(sum(!(SpikeInfo[,1] %in% GeneNames[Tech])) > 0) stop("'SpikeInfo' includes spikes that are not in the Counts matrix")
    matching <- match(GeneNames[Tech], SpikeInfo[,1])
    SpikeInput <- SpikeInfo[matching,2]    
  }
  else {SpikeInput = 1}
  
  # Checks for creating the SingleCellExperiment class
  errors <- character()
  
  if(!(is.numeric(Counts) & all(Counts>=0) & 
       sum(!is.finite(Counts))==0 )) errors <- c(errors, "Invalid value for 'Counts'")
  if(sum(Counts %% 1) > 0) errors <- c(errors, "Invalid value for 'Counts' (entries must be positive integers)")          
  if(!(is.logical(Tech))) errors <- c(errors, "Invalid value for 'Tech'")
  if(!(is.numeric(SpikeInput) & all(SpikeInput>0) & 
       sum(!is.finite(SpikeInput))==0 )) errors <- c(errors, "Invalid value for 'SpikeInput'.")
  
  q = nrow(Counts)
  q.bio = q - length(SpikeInput)
  n = ncol(Counts)
  
  # Checks valid for datasets with spikes only
  if(length(SpikeInput) > 1)
  {
    if(!( length(Tech) == q & sum(!Tech) == q.bio )) 
      errors <- c(errors, "Argument's dimensions are not compatible.")
    
    if(sum(apply(Counts[Tech,],2,sum) == 0) > 0) 
      errors <- c(errors, "Some cells have zero reads mapping back to the spike-in genes. Please remove them before creating the SingleCellExperiment object.")
    
    if(sum(apply(Counts[!Tech,],2,sum) == 0) > 0) 
      errors <- c(errors, "Some cells have zero reads mapping back to the intrinsic genes. Please remove them before creating the SingleCellExperiment object.")
    
    if(!( sum(Tech[1:q.bio]) == 0 & sum(Tech[(q.bio+1):q])==q-q.bio )) 
      errors <- c(errors, "Expression counts are not in the right format (spike-in genes must be at the bottom of the matrix).")
  }
  # Checks valid for datasets with no spikes only
  else
  {
    if(sum(apply(Counts,2,sum) == 0) > 0) 
      errors <- c(errors, "Some cells have zero reads mapping back to the intrinsic genes. Please remove them before creating the SingleCellExperiment object.")
    
    if(length(unique(BatchInfo)) == 1)
      errors <- c(errors, "If spike-in genes are not available, BASiCS requires the data to contain at least 2 batches of cells (for the same population)")
  }
  
  # Checks valid for any data
  if(length(Tech) != q)
    errors <- c(errors, "Argument's dimensions are not compatible.")
  
  if(length(GeneNames) != q)
    errors <- c(errors, "Incorrect length of the vector stored in the GeneNames slot.")
  
  if(sum(apply(Counts,1,sum) == 0) > 0) 
    warning("Some genes have zero counts across all cells. Unless running a differential expression analysis, please remove those genes. Otherwise, the BASiCS_Data object is still a valid object. However, due to the lack of counts, posterior estimates for mu[i] and delta[i] associated to those genes will be driven by the prior. In such case, you must specify `PriorDelta = 'log-normal' in BASiCS_MCMC function. ")
  
  if(length(BatchInfo) != n) 
    errors <- c(errors, "BatchInfo slot is not compatible with the number of cells contained in Counts slot.")
  
  if (length(errors) == 0) TRUE else stop(errors)
  
  # Create a SingleCellExperiment data object
  Data <- SingleCellExperiment::SingleCellExperiment(assays = list(Counts = as.matrix(Counts)),
                               metadata = list(SpikeInput = SpikeInput, BatchInfo = BatchInfo))
  isSpike(Data, "ERCC") <- Tech 
  
  show(Data)
  message('\n',
          'NOTICE: BASiCS requires a pre-filtered dataset \n',
          '    - You must remove poor quality cells before creating the BASiCS data object \n',
          '    - We recommend to pre-filter very lowly expressed transcripts before creating the object. \n',
          '      Inclusion criteria may vary for each data. For example, remove transcripts \n',
          '          - with very low total counts across of all of the samples \n',
          '          - that are only expressed in a few cells \n',
          '            (by default genes expressed in only 1 cell are not accepted) \n',
          '          - with very low total counts across the samples where the transcript is expressed \n',
          '\n',
          ' BASiCS_Filter can be used for this purpose \n')
  return(Data)
}