#' @title Creates a SingleCellExperiment object from a matrix of expression
#' counts and experimental information about spike-in genes
#'
#' @description \code{newBASiCS_Data} creates a
#' \code{\linkS4class{SingleCellExperiment}} object from a matrix of expression
#' counts and experimental information about spike-in genes.
#'
#' @param Counts Matrix of dimensions \code{q} times \code{n} whose elements
#' contain the expression counts to be analyses
#' (including biological and technical spike-in genes). Gene names must be
#' stored as \code{rownames(Counts)}.
#' @param Tech Logical vector of length \code{q}. If \code{Tech = FALSE} the
#' gene is biological; otherwise the gene is spike-in. Defaul value:
#' \code{Tech = rep(FALSE, nrow(Counts))}.
#' @param SpikeInfo \code{data.frame} whose first and second columns contain
#' the gene names assigned to the spike-in genes (they must match the ones in
#' \code{rownames(Counts)}) and the associated input number of molecules,
#' respectively. If \code{SpikeInfo = NULL}, only the horizontal integration
#' implementation (no spikes) can be run. Default value: \code{SpikeInfo = NULL}.
#' @param BatchInfo Vector of length \code{n} whose elements indicate batch
#' information. Not required if a single batch is present on the data.
#' Default value: \code{BatchInfo = NULL}.
#' @param SpikeType Character to indicate what type of spike-ins are in use.
#' Default value: \code{SpikeType = "ERCC"} (parameter is no longer used).
#'
#' @return An object of class \code{\linkS4class{SingleCellExperiment}}.
#'
#' @examples
#'
#' ## Data with spike-ins
#'
#' # Expression counts
#' set.seed(1)
#' Counts <- matrix(rpois(50*10, 2), ncol = 10)
#' rownames(Counts) <- c(paste0('Gene', 1:40), paste0('ERCC', 1:10))
#' # Technical information
#' Tech <- grepl("ERCC", rownames(Counts))
#' # Spikes input number of molecules
#' set.seed(2)
#' SpikeInfo <- data.frame(gene = rownames(Counts)[Tech],
#'                         amount = rgamma(10, 1, 1))
#'
#' # Creating a BASiCS_Data object (no batch effect)
#' DataExample <- newBASiCS_Data(Counts, Tech = Tech, SpikeInfo = SpikeInfo)
#'
#' # Creating a BASiCS_Data object (with batch effect)
#' BatchInfo <- c(rep(1, 5), rep(2, 5))
#' DataExample <- newBASiCS_Data(Counts, Tech = Tech,
#'                               SpikeInfo = SpikeInfo, BatchInfo = BatchInfo)
#'
#' ## Data without spike-ins (BatchInfo is required)
#'
#' # Expression counts
#' set.seed(1)
#' Counts <- matrix(rpois(50*10, 2), ncol = 10)
#' rownames(Counts) <- paste0('Gene', 1:50)
#' BatchInfo <- c(rep(1, 5), rep(2, 5))
#'
#' # Creating a BASiCS_Data object (with batch effect)
#' DataExample <- newBASiCS_Data(Counts, BatchInfo = BatchInfo)
#'
#' @seealso \code{\linkS4class{SingleCellExperiment}}
#'
#' @author Catalina A. Vallejos \email{cnvallej@@uc.cl}
#' @author Nils Eling \email{eling@@ebi.ac.uk}
#'
#' @export
newBASiCS_Data <- function(Counts,
                           Tech = rep(FALSE, nrow(Counts)),
                           SpikeInfo = NULL,
                           BatchInfo = NULL,
                           SpikeType = "ERCC") {

  if (!inherits(Counts, "matrix")) {
    stop("Counts must be a matrix.")
  }
  # Separating intrinsic from spike-in transcripts
  CountsBio <- Counts[!Tech, ]
  CountsTech <- Counts[Tech, ]
  # Extracting gene labels
  GeneName <- rownames(Counts)
  
  # Create a SingleCellExperiment data object
  Data <- SingleCellExperiment(assays = list(counts = as.matrix(CountsBio)))
  colnames(Data) <- colnames(CountsBio)
  rownames(Data) <- rownames(CountsBio)  
  
  # Adding metadata associated to batch information
  ## Setting a default value for BatchInfo when absent
  if (is.null(BatchInfo)) {
    BatchInfo <- rep(1, times = ncol(Counts))
  }  
  colData(Data) <- S4Vectors::DataFrame("BatchInfo" = BatchInfo)
  
  # Adding spike-ins information
  if (!is.null(SpikeInfo)) {
    # If SpikeInfo is provided, run validity checks
    if (!is.data.frame(SpikeInfo) | data.table::is.data.table(SpikeInfo)) {
      stop("'SpikeInfo' must be a 'data.frame'")
    }
    if (sum(Tech) == 0) {
      stop("'SpikeInfo' provided but no genes were marked as spike-ins \n",
           "Please revise the input value provided for 'Tech'")
    }
    # Extracting spike-in input molecules in the correct order
    if (sum(!(GeneName[Tech] %in% SpikeInfo[, 1])) > 0) {
      stop("'SpikeInfo' is missing information for some of the spikes")
    }
    if (sum(!(SpikeInfo[, 1] %in% GeneName[Tech])) > 0) {
      stop("'SpikeInfo' includes spikes that are not in 'Counts'")
    }
    if(any(GeneName[Tech] != SpikeInfo[, 1])) {
      # Re-order spike-ins in SpikeInfo if required
      matching <- match(GeneName[Tech], SpikeInfo[, 1])
      SpikeInfo <- SpikeInfo[matching, ]
    }
    metadata(Data)$SpikeInput <- SpikeInfo
    altExp(Data, "spike-ins") <- SummarizedExperiment(CountsTech)
    WithSpikes <- TRUE
  }  else {
    message("The data does not contain spike-in genes")
    WithSpikes <- FALSE
  }
  # Checks to assess if the data contains the required information
  errors <- HiddenChecksBASiCS_Data(Data, WithSpikes)
  if (length(errors) > 0) stop(errors) 
  
  message("\n", "NOTICE: BASiCS requires a pre-filtered dataset \n",
            "    - You must remove poor quality cells before hand \n",
            "    - We recommend to pre-filter lowly expressed transcripts. \n",
            "      Inclusion criteria may vary for each data. \n",
            "      For example, remove transcripts: \n",
            "          - with low total counts across of all of the cells \n",
            "          - that are only expressed in a few cells \n",
            "            (genes expressed in only 1 cell are not accepted) \n",
            "\n BASiCS_Filter can be used for this purpose. \n")
  Data
}
