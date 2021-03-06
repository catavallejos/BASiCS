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
#' implementation (no spikes) can be run. Default value:
#' \code{SpikeInfo = NULL}.
#' @param BatchInfo Vector of length \code{n} whose elements indicate batch
#' information. Not required if a single batch is present on the data.
#' Default value: \code{BatchInfo = NULL}.
#' @param SpikeType Character to indicate what type of spike-ins are in use.
#' Default value: \code{SpikeType = "ERCC"} (parameter is no longer used).
#'
#' @return An object of class \code{\linkS4class{SingleCellExperiment}}.
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
  .Deprecated()
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
  
  # Adding metadata associated to batch information
  ## Setting a default value for BatchInfo when absent
  if (is.null(BatchInfo)) {
    BatchInfo <- rep(1, times = ncol(Counts))
  }  
  colData(Data) <- S4Vectors::DataFrame(
    BatchInfo = BatchInfo,
    row.names = colnames(CountsBio)
  )
  colnames(Data) <- colnames(CountsBio)
  rownames(Data) <- rownames(CountsBio)
  
  # Adding spike-ins information
  if (!is.null(SpikeInfo)) {
    # If SpikeInfo is provided, run validity checks
    if (!is.data.frame(SpikeInfo)) {
      stop("'SpikeInfo' must be a 'data.frame'")
    }
    if (sum(Tech) == 0) {
      stop(
        "'SpikeInfo' provided but no genes were marked as spike-ins \n",
        "Please revise the input value provided for 'Tech'"
      )
    }
    # Extracting spike-in input molecules in the correct order
    if (any(!(GeneName[Tech] %in% SpikeInfo[, 1]))) {
      stop("'SpikeInfo' is missing information for some of the spikes")
    }
    if (any(!(SpikeInfo[, 1] %in% GeneName[Tech]))) {
      stop("'SpikeInfo' includes spikes that are not in 'Counts'")
    }
    if (any(GeneName[Tech] != SpikeInfo[, 1])) {
      # Re-order spike-ins in SpikeInfo if required
      matching <- match(GeneName[Tech], SpikeInfo[, 1])
      SpikeInfo <- SpikeInfo[matching, ]
    }
    altExp(Data, "spike-ins") <- SummarizedExperiment(CountsTech, 
      rowData = SpikeInfo)
    WithSpikes <- TRUE
  } else {
    message("The data does not contain spike-in genes")
    WithSpikes <- FALSE
  }
  # Checks to assess if the data contains the required information
  errors <- .ChecksBASiCS_Data(Data, WithSpikes)
  if (length(errors) > 0) stop(errors)
  
  message(
    "\n", "NOTICE: BASiCS requires a pre-filtered dataset \n",
    "    - You must remove poor quality cells before hand \n",
    "    - We recommend to pre-filter lowly expressed transcripts. \n",
    "      Inclusion criteria may vary for each data. \n",
    "      For example, remove transcripts: \n",
    "          - with low total counts across of all of the cells \n",
    "          - that are only expressed in a few cells \n",
    "            (genes expressed in only 1 cell are not accepted) \n",
    "\n BASiCS_Filter can be used for this purpose. \n"
  )
  
  Data
}
