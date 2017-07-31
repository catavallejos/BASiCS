#' @title Filter for input datasets
#'
#' @description \code{BASiCS_Filter} indicates which transcripts and cells pass a pre-defined inclusion criteria. The output of this function can be combined with \code{newBASiCS_Data} to generate a 
#' the \code{\linkS4class{SummarizedExperiment}} object required to run BASiCS. 
#'
#' @param Counts Matrix of dimensions \code{q} times \code{n} whose elements corresponds to the simulated expression counts.
#' First \code{q.bio} rows correspond to biological genes. Last \code{q-q.bio} rows correspond to technical spike-in genes.
#' @param Tech Logical vector of length \code{q}. If \code{Tech = F} the gene is biological; otherwise the gene is spike-in.
#' @param SpikeInput Vector of length \code{q-q.bio} whose elements indicate the simulated input concentrations for the spike-in genes.
#' @param BatchInfo Vector of length \code{n} whose elements indicate batch information. Not required if a single batch is present on the data. Default value: \code{BatchInfo = NULL}. 
#' @param MinTotalCountsPerCell Minimum value of total expression counts required per cell (biological and technical). Default value: \code{MinTotalCountsPerCell = 2}.
#' @param MinTotalCountsPerGene Minimum value of total expression counts required per transcript (biological and technical). Default value: \code{MinTotalCountsPerGene = 2}.
#' @param MinCellsWithExpression Minimum number of cells where expression must be detected (positive count). Criteria applied to each transcript. Default value: \code{MinCellsWithExpression = 2}.
#' @param MinAvCountsPerCellsWithExpression Minimum average number of counts per cells where expression is detected. Criteria applied to each transcript. Default value: \code{MinAvCountsPerCellsWithExpression = 2}.
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
#' Counts = Counts = matrix(rpois(50*10, 2), ncol = 10)
#' rownames(Counts) <- c(paste0("Gene", 1:40), paste0("Spike", 1:10))
#' Tech = c(rep(FALSE,40),rep(TRUE,10))
#' set.seed(2)
#' SpikeInput = rgamma(10,1,1)
#' SpikeInfo <- data.frame("SpikeID" = paste0("Spike", 1:10), "SpikeInput" = SpikeInput)
#'
#' Filter = BASiCS_Filter(Counts, Tech, SpikeInput,
#'                        MinTotalCountsPerCell = 2, MinTotalCountsPerGene = 2,
#'                        MinCellsWithExpression = 2, MinAvCountsPerCellsWithExpression = 2)
#' SpikeInfoFilter = SpikeInfo[SpikeInfo$SpikeID %in%
#'          names(Filter$IncludeGenes)[Filter$IncludeGenes == TRUE],]
#' FilterData = newBASiCS_Data(Filter$Counts, Filter$Tech, SpikeInfoFilter)
#'
#' @author Catalina A. Vallejos \email{cnvallej@@uc.cl}
#'
#' @references 
#' Vallejos, Marioni and Richardson (2015). Bayesian Analysis of Single-Cell Sequencing data. PLoS Computational Biology. 
#' 
#' Vallejos, Marioni and Richardson (2016). Beyond comparisons of means: understanding changes in gene expression at the single-cell level. Genome Biology.
BASiCS_Filter <- function(Counts, Tech, SpikeInput, BatchInfo = NULL,
                          MinTotalCountsPerCell = 2, MinTotalCountsPerGene = 2,
                          MinCellsWithExpression = 2, MinAvCountsPerCellsWithExpression = 2)
{
  q = length(Tech)
  q.bio = q - sum(SpikeInput)
  n = ncol(Counts)
  
  CellIndex = 1:n
  GeneIndex = 1:q
  
  # Remove cells with zero counts in either biological or technical genes
  IncludeCells = ifelse(apply(Counts[!Tech,],2,sum)>0 & apply(Counts[Tech,],2,sum)>0, TRUE, FALSE)
  if(sum(IncludeCells) == 0) stop('All cells have zero biological or technical counts (across all transcripts) \n')
  IncludeCells = ifelse(apply(Counts,2,sum) >= MinTotalCountsPerCell, IncludeCells, F)
  Counts1 = Counts[,IncludeCells]
  
  # Remove transcripts with zero counts across all cells
  IncludeBio = ifelse(apply(Counts1[!Tech,],1,sum) >= MinTotalCountsPerGene, TRUE, FALSE)
  IncludeTech = ifelse(apply(Counts1[Tech,],1,sum) >= MinTotalCountsPerGene, TRUE, FALSE)
  
  # Remove transcripts expressed in less than 'MinExpressedCells' cells
  NonZero <- I(Counts1>0)
  IncludeBio = ifelse(apply(NonZero[!Tech,], 1, sum) >= MinCellsWithExpression, IncludeBio, FALSE)
  IncludeTech = ifelse(apply(NonZero[Tech,], 1, sum) >= MinCellsWithExpression, IncludeTech, FALSE)
  
  # Remove transcripts with low counts in the cells where they are expressed
  if(min(apply(NonZero[!Tech,], 1, sum)) == 0 & MinCellsWithExpression == 0) warning("Some genes have zero counts in all cells. These should be removed before running the analysis (use 'MinCellsWithExpression' > 0).")
  IncludeBio = ifelse(apply(Counts1[!Tech,], 1, sum) >= MinAvCountsPerCellsWithExpression * apply(NonZero[!Tech,], 1, sum), IncludeBio, FALSE)
  IncludeTech = ifelse(apply(Counts1[Tech,], 1, sum) >= MinAvCountsPerCellsWithExpression * apply(NonZero[Tech,], 1, sum), IncludeTech, FALSE)

  if(!is.null(BatchInfo))
  {
    list("Counts" = Counts1[c(IncludeBio,IncludeTech),], 
         "Tech" = Tech[c(IncludeBio,IncludeTech)],
         "SpikeInput" = SpikeInput[IncludeTech], 
         "BatchInfo" = BatchInfo[IncludeCells],
         "IncludeGenes" = c(IncludeBio,IncludeTech), 
         "IncludeCells" = IncludeCells)
  }
  else
  {
    list("Counts" = Counts1[c(IncludeBio,IncludeTech),], 
         "Tech" = Tech[c(IncludeBio,IncludeTech)],
         "SpikeInput" = SpikeInput[IncludeTech], 
         "BatchInfo" = NULL,
         "IncludeGenes" = c(IncludeBio,IncludeTech), 
         "IncludeCells" = IncludeCells)    
  }
  
}