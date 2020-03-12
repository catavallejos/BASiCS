#' @title Filter for input datasets
#'
#' @description \code{BASiCS_Filter} indicates which transcripts and
#' cells pass a pre-defined inclusion criteria. The output of this
#' function can be combined with \code{newBASiCS_Data} to generate a
#' the \code{\linkS4class{SingleCellExperiment}} object required to run BASiCS.
#' For more systematic tools for quality control, please refer to the
#' \code{scater} Bioconductor package.
#'
#' @param Counts Matrix of dimensions \code{q} times \code{n} whose elements
#' corresponds to the simulated expression counts.
#' First \code{q.bio} rows correspond to biological genes. Last \code{q-q.bio}
#' rows correspond to technical spike-in genes.
#' @param Tech Logical vector of length \code{q}. If \code{Tech = FALSE}
#' the gene is biological; otherwise the gene is spike-in.
#' @param SpikeInput Vector of length \code{q-q.bio} whose elements indicate
#' the simulated input concentrations for the spike-in genes.
#' @param BatchInfo Vector of length \code{n} whose elements indicate batch
#' information. Not required if a single batch is present on the data.
#' Default: \code{BatchInfo = NULL}.
#' @param MinTotalCountsPerCell Minimum value of total expression counts
#' required per cell (biological and technical).
#' Default: \code{MinTotalCountsPerCell = 2}.
#' @param MinTotalCountsPerGene Minimum value of total expression counts
#' required per transcript (biological and technical).
#' Default: \code{MinTotalCountsPerGene = 2}.
#' @param MinCellsWithExpression Minimum number of cells where expression
#' must be detected (positive count). Criteria applied to each transcript.
#' Default: \code{MinCellsWithExpression = 2}.
#' @param MinAvCountsPerCellsWithExpression Minimum average number of
#' counts per cells where expression is detected. Criteria applied to
#' each transcript. Default value: \code{MinAvCountsPerCellsWithExpression = 2}.
#'
#' @return A list of 2 elements
#' \describe{
#' \item{\code{Counts}}{Filtered matrix of expression counts}
#' \item{\code{Tech}}{Filtered vector of spike-in indicators}
#' \item{\code{SpikeInput}}{Filtered vector of spike-in genes input molecules}
#' \item{\code{BatchInfo}}{Filtered vector of the 'BatchInfo' argument}
#' \item{\code{IncludeGenes}}{Inclusion indicators for transcripts}
#' \item{\code{IncludeCells}}{Inclusion indicators for cells}
#' }
#'
#' @examples
#'
#' set.seed(1)
#' Counts <- matrix(rpois(50*10, 2), ncol = 10)
#' rownames(Counts) <- c(paste0('Gene', 1:40), paste0('Spike', 1:10))
#' Tech <- c(rep(FALSE,40),rep(TRUE,10))
#' set.seed(2)
#' SpikeInput <- rgamma(10,1,1)
#' SpikeInfo <- data.frame('SpikeID' = paste0('Spike', 1:10),
#'                         'SpikeInput' = SpikeInput)
#'
#' Filter <- BASiCS_Filter(Counts, Tech, SpikeInput,
#'                         MinTotalCountsPerCell = 2,
#'                         MinTotalCountsPerGene = 2,
#'                         MinCellsWithExpression = 2,
#'                         MinAvCountsPerCellsWithExpression = 2)
#' SpikeInfoFilter <- SpikeInfo[SpikeInfo$SpikeID %in% rownames(Filter$Counts),]
#' FilterData <- newBASiCS_Data(Filter$Counts, Filter$Tech, SpikeInfoFilter)
#'
#' @author Catalina A. Vallejos \email{cnvallej@@uc.cl}
#'
#' @export
BASiCS_Filter <- function(Counts,
                          Tech = rep(FALSE, nrow(Counts)),
                          SpikeInput = NULL,
                          BatchInfo = NULL,
                          MinTotalCountsPerCell = 2,
                          MinTotalCountsPerGene = 2,
                          MinCellsWithExpression = 2,
                          MinAvCountsPerCellsWithExpression = 2) {
  q <- length(Tech)
  n <- ncol(Counts)
  CellIndex <- seq_len(n)
  GeneIndex <- seq_len(q)
  colSumsAll <- Matrix::colSums(Counts)
  if((sum(Tech) > 0) & is.null(SpikeInput)) {
    stop("`SpikeInput` is required when the data contains spike-ins")
  }

  # Cell filter
  IncludeCells <- rep(TRUE, times = n)
  if(sum(Tech) > 0) {
    colSumsBio <- Matrix::colSums(Counts[!Tech, ])
    colSumsTech <- Matrix::colSums(Counts[Tech, ])
    # Remove cells with zero counts in either biological or technical genes
    IncludeCells[which((colSumsBio == 0) | (colSumsTech == 0))] <- FALSE
    if (sum(IncludeCells) == 0) {
      stop("All cells have zero biological or technical counts \n")
    }
  }
  IncludeCells[which(colSumsAll < MinTotalCountsPerCell)] <- FALSE
  Counts1 <- Counts[, IncludeCells]

  # Remove transcripts with low total counts across all cells
  IncludeGenes <- rep(TRUE, length = q)
  rowSumsAll <- matrixStats::rowSums2(Counts1)
  IncludeGenes[which(rowSumsAll < MinTotalCountsPerGene)] <- FALSE

  # Remove transcripts expressed in less than 'MinExpressedCells' cells
  rowSumsNonZero <- matrixStats::rowSums2(Counts1 > 0)
  IncludeGenes[which(rowSumsNonZero < MinCellsWithExpression)] <- FALSE

  # Remove transcripts with low counts in the cells where they are expressed
  if ((min(rowSumsNonZero) == 0) & (MinCellsWithExpression == 0)) {
    warning("Some genes have zero counts in all cells. \n",
            "These should be removed before running the analysis \n",
            "(use 'MinCellsWithExpression' > 0).")
  }
  ind <- which(rowSumsAll < MinAvCountsPerCellsWithExpression * rowSumsNonZero)
  IncludeGenes[ind] <- FALSE

  if (!is.null(BatchInfo)) { BatchInfo <- BatchInfo[IncludeCells] }
  if (!is.null(SpikeInput)) {
    IncludeTech <- IncludeGenes[Tech]
    SpikeInput <- SpikeInput[IncludeTech]
  }

  list(Counts = Counts1[IncludeGenes, ],
       Tech = Tech[IncludeGenes],
       SpikeInput = SpikeInput,
       BatchInfo = BatchInfo,
       IncludeGenes = IncludeGenes,
       IncludeCells = IncludeCells)
}
