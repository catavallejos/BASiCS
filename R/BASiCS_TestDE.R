#' @name BASiCS_TestDE
#' 
#' @title Detection of genes with changes in expression
#' 
#' @description Function to assess changes in expression (mean and over-dispersion)
#' 
#' @param Chain1 an object of class \code{\link[BASiCS]{BASiCS_Chain-class}} containing parameter estimates for the first group of cells
#' @param Chain2 an object of class \code{\link[BASiCS]{BASiCS_Chain-class}} containing parameter estimates for the second group of cells
#' @param EpsilonM Minimum fold change tolerance threshold for detecting changes in overall expression (must be a positive real number)
#' @param EpsilonD Minimum fold change tolerance threshold for detecting changes in cell-to-cell biological over dispersion (must be a positive real number)
#' @param EviThresholdM Optional parameter. Evidence threshold for detecting changes in overall expression (must be a positive value, between 0 and 1)
#' @param EviThresholdD Optional parameter. Evidence threshold for detecting changes in cell-to-cell biological over dispersion (must be a positive value, between 0 and 1)
#' @param OrderVariable Ordering variable for output. Must take values in \code{c("GeneIndex", "GeneNames", "Prob")}.
#' @param GroupLabel1 Label assigned to reference group. Default: \code{GroupLabel1 = "Group1"}
#' @param GroupLabel2 Label assigned to reference group. Default: \code{GroupLabel2 = "Group2"}
#' @param Plot If \code{Plot = TRUE}, MA and volcano plots are generated. 
#' @param PlotOffset If \code{Plot = TRUE}, the offset effect is visualised.   
#' @param OffSet Optional argument to remove a fix offset effect (if not previously removed from the MCMC chains). This argument will be removed shorly, once offset removal is built as an internal step. 
#' @param EFDR_M Target for expected false discovery rate related to the comparison of means (default = 0.05)
#' @param EFDR_D Target for expected false discovery rate related to the comparison of dispersions (default = 0.05)
#' @param GenesSelect Optional argument to provide a user-defined list of genes to be considered for the comparison (default = NULL). When used, this argument must be a vector of 'TRUE' (include gene) / 'FALSE' (exclude gene) indicator, with the same length as the number of intrinsic genes and following the same order as how genes are displayed in the table of counts.  This argument is necessary in order to have a meaningful EFDR calibration when the user decides to exclude some genes from the comparison. 
#' @param ... Graphical parameters (see \code{\link[graphics]{par}}).
#' 
#' @return \code{BASiCS_TestDE} returns a list of 4 elements:
#' \describe{
#' \item{\code{TableMean}}{A \code{\link[base]{data.frame}} containing the results of the differential mean test}
#'    \describe{
#'    \item{\code{GeneNames}}{Gene name}
#'    \item{\code{MeanOverall}}{For each gene, the estimated mean expression parameter \eqn{\mu[i]} is averaged across both groups of cells (weighted by sample size).}
#'    \item{\code{Mean1}}{Estimated mean expression parameter \eqn{\mu[i]} for each biological gene in the first group of cells.}
#'    \item{\code{Mean2}}{Estimated mean expression parameter \eqn{\mu[i]} for each biological gene in the second group of cells.}
#'    \item{\code{MeanFC}}{Fold change in mean expression parameters between the first and second groups of cells.}
#'    \item{\code{MeanLog2FC}}{Log2-transformed fold change in mean expression between the first and second groups of cells.}
#'    \item{\code{ProbDiffMean}}{Posterior probability for mean expression difference between the first and second groups of cells.}
#'    \item{\code{ResultDiffExp}}{Indicator if a gene has a higher mean expression in the first or second groups of cells.}
#'    }
#' \item{\code{TableDisp}}{A \code{\link[base]{data.frame}} containing the results of the differential dispersion test (excludes genes for which the mean changes).}
#'    \describe{
#'    \item{\code{GeneNames}}{Gene name}
#'    \item{\code{MeanOverall}}{For each gene, the estimated mean expression parameter \eqn{\mu[i]} is averaged across both groups of cells (weighted by sample size).}
#'    \item{\code{DispOverall}}{For each gene, the estimated over-dispersion parameter \eqn{\delta[i]} is averaged across both groups of cells (weighted by sample size).}
#'    \item{\code{Disp1}}{Estimated over-dispersion parameter \eqn{\delta[i]} for each biological gene in the first group of cells.}
#'    \item{\code{Disp2}}{Estimated over-dispersion parameter \eqn{\delta[i]} for each biological gene in the second group of cells.}
#'    \item{\code{DispFC}}{Fold change in over-dispersion parameters between the between the first and second groups of cells.}
#'    \item{\code{DispLog2FC}}{Log-transformed fold change in over-dispersion between the first and second groups of cells.}
#'    \item{\code{ProbDiffDisp}}{Posterior probability for over-dispersion difference between the first and second groups of cells.}
#'    \item{\code{ResultDiffDisp}}{Indicator if a gene has a higher over-dispersion in the first or second groups of cells.}
#'    }
#'  \item{\code{DiffExpSummary}}{A list containing the following information for the differential mean expression test:}
#'    \describe{
#'   \item{\code{EviThreshold}}{Evidence thresholds.}
#'   \item{\code{EFDR}}{Expected false discovery rate for the given thresholds.}
#'   \item{\code{EFNR}}{Expected false negative rate for the given thresholds.}
#'   }
#'  \item{\code{DiffOverDispSummary}}{A list containing the following information for the differential over-dispersion test:}
#'    \describe{
#'   \item{\code{EviThreshold}}{Evidence thresholds.}
#'   \item{\code{EFDR}}{Expected false discovery rate for the given thresholds.}
#'   \item{\code{EFNR}}{Expected false negative rate for the given thresholds.}
#'   }
#'  \item{\code{Chain1_offset}}{an \code{\link[BASiCS]{BASiCS_Chain-class}} object: \code{Chain1} after offset removal.}
#'  \item{\code{Chain2_offset}}{an \code{\link[BASiCS]{BASiCS_Chain-class}} object: \code{Chain2} after offset removal
#'  (this is only provided for completeness; \code{Chain2} is not affected by the offset).}
#'  \item{\code{OffsetChain}}{MCMC chain calculated for the offset effect.}
#'  \item{\code{Offset}}{Estimated offset (posterior median of \code{OffsetChain}). Default value set equal to 1 when offset correction is not performed.}
#' }
#'  
#' @examples
#' 
#' # Loading two 'BASiCS_Chain' objects (obtained using the 'BASiCS_MCMC' function)
#' data(ChainSC)
#' data(ChainRNA)
#' 
#' Test <- BASiCS_TestDE(Chain1 = ChainSC, Chain2 = ChainRNA,
#'                       GroupLabel1 = "SC", GroupLabel2 = "P&S",
#'                       EpsilonM = 0.4*log(2), EpsilonD = 0.4*log(2), OffSet = TRUE)
#'                       
#' # Results for the differential mean test
#' head(Test$TableMean)
#' 
#' # Results for the differential over-dispersion test
#' # This only includes genes marked as 'NoDiff' in Test$TableMean 
#' head(Test$TableDisp)                    
#' 
#' @author Catalina A. Vallejos \email{cnvallej@@uc.cl} and Nils Eling
#' 
#' @rdname BASiCS_TestDE
BASiCS_TestDE <- function(Chain1,
                          Chain2,
                          EpsilonM = 0.40 * log(2),
                          EpsilonD = 0.40 * log(2),
                          EviThresholdM = NULL,
                          EviThresholdD = NULL,
                          OrderVariable = "Prob",
                          GroupLabel1 = "Group1",
                          GroupLabel2 = "Group2",
                          Plot = TRUE, 
                          PlotOffset = TRUE, 
                          OffSet = TRUE, 
                          EFDR_M = 0.05,
                          EFDR_D = 0.05,
                          GenesSelect = NULL, 
                            ...)
{
  # Checking validity of input arguments
  if(!is(Chain1,"BASiCS_Chain")) stop("'Chain1' is not a 'BASiCS_Chain' class object.") 
  if(!is(Chain2,"BASiCS_Chain")) stop("'Chain2' is not a 'BASiCS_Chain' class object.") 
  
  # Test compatibility of both BASiCS_Chain objects
  if(ncol(Chain1@mu) != ncol(Chain2@mu)) stop("The two BASiCS_Chain objects do not contain the same number of genes.")
  if(!identical(colnames(Chain1@mu), colnames(Chain2@mu))) stop("The gene names of both BASiCS_Chain objects are not in the same order.")
  
  GeneNames = colnames(Chain1@mu)
  
  if(EpsilonM < 0 | !is.finite(EpsilonM)) stop("Minimum tolerance of fold-change 'EpsilonM' must be a positive real value")
  if(EpsilonD < 0 | !is.finite(EpsilonD)) stop("Minimum tolerance of fold-change 'EpsilonD' must be a positive real value")
  if(!is.logical(Plot) | length(Plot)!=1) stop("Please insert TRUE or FALSE for Plot parameter")
  if(!is.null(EviThresholdM) | !is.null(EviThresholdD))
  {
    if(EviThresholdM < 0 | EviThresholdM > 1 | !is.finite(EviThresholdM)) 
      stop("Evidence threshold 'EviThresholdM' must be contained in (0,1) \n For automatic threshold search use EviThresholdM = NULL.")    
    if(EviThresholdD < 0 | EviThresholdD > 1 | !is.finite(EviThresholdD)) 
      stop("Evidence threshold 'EviThresholdD' must be contained in (0,1) \n For automatic threshold search use EviThresholdD = NULL.")    
  }
  if(!(OrderVariable %in% c("GeneIndex", "GeneNames", "Prob"))) stop("Invalid 'OrderVariable' value")
  if(!is.character(GroupLabel1) | length(GroupLabel1) > 1) stop("Invalid value for 'GroupLabel1'")
  if(!is.character(GroupLabel2) | length(GroupLabel2) > 1) stop("Invalid value for 'GroupLabel2'")
  if(!is.null(GenesSelect) & (length(GenesSelect) != length(GeneNames))) stop("Invalid value for 'GenesSelect'")
  if(!is.null(GenesSelect) & !is.logical(GenesSelect)) stop("Invalid value for 'GenesSelect'")
  # if(!is.null(GenesSelect) & (length(GenesSelect) == sum(!GenesSelect))) stop("Invalid value for 'GenesSelect'")

  message("Log-fold change thresholds are now set in a log2 scale. \n Previous versions of BASiCS used a natural logarithm scale.")
  
  n1 = ncol(Chain1@nu)
  n2 = ncol(Chain2@nu)
  n = n1 + n2
  # With offset correction
  if(OffSet)
  {
    # Calculating iteration-specific offset
    OffsetChain <- rowSums(Chain1@mu)/rowSums(Chain2@mu)
    # Offset point estimate
    OffsetEst <- median(OffsetChain)
    
    # Offset correction
    Chain1_offset <- Chain1
    Chain1_offset@mu <- Chain1@mu / OffsetEst
    Chain1_offset@phi <- Chain1@phi * OffsetEst
    Chain2_offset <- Chain2 # Chain2 requires no change
    Summary1 <- Summary(Chain1_offset)
    Summary2 <- Summary(Chain2_offset)
    
    # Pre-offset correction LFC estimates
    MuBase_old = (matrixStats::colMedians(Chain1@mu) * n1 + matrixStats::colMedians(Chain2@mu) * n2)/n
    ChainTau_old = log2(Chain1@mu  / Chain2@mu  )
    MedianTau_old = matrixStats::colMedians(ChainTau_old)

    # Offset corrected LFC estimates
    MuBase = (Summary1@mu[,1] * n1 + Summary2@mu[,1] * n2)/n
    ChainTau = log2(Chain1_offset@mu  / Chain2_offset@mu  )
    MedianTau = matrixStats::colMedians(ChainTau)
    
    if(!PlotOffset)
    {
      message(paste("--------------------------------------------------------------------- \n",
                    "Offset estimate:", round(OffsetEst,4), "(ratio", GroupLabel1, "vs", GroupLabel2, ").\n",
                    "To visualise its effect, please use 'PlotOffset = TRUE'.\n",
                    "--------------------------------------------------------------------- \n"))      
    }
    else
    {    
      message(paste("--------------------------------------------------------------------- \n",
                    "Offset estimate:", round(OffsetEst,4), "(ratio", GroupLabel1, "vs", GroupLabel2, ").\n",
                    "--------------------------------------------------------------------- \n"))
    }

    
    if(PlotOffset == TRUE)
    {
      par(ask=TRUE)
      # Offset uncertainty
      message("Posterior uncertainty for offset estimate ... ")
      graphics::boxplot(OffsetChain, frame = FALSE, 
              main = "Offset MCMC chain", ylab = "Offset estimate") 
      # Mean expression parameters before/after offset correction
      message("Mean expression parameters before/after offset correction ...")
      par(mfrow = c(1,2))
      graphics::boxplot(cbind(matrixStats::colMedians(Chain1@mu), Summary2@mu[,1]), 
              frame = FALSE, main = "Before correction", 
              names = c(GroupLabel1, GroupLabel2),
              ylab = "Mean expression", log = "y")
      graphics::boxplot(cbind(Summary1@mu[,1], Summary2@mu[,1]), 
              frame = FALSE, main = "After correction", 
              names = c(GroupLabel1, GroupLabel2),
              ylab = "Mean expression", log = "y")
      # MA plot pre/after offset
      message("MA plot before/after offset correction ...")
      par(mfrow = c(1,2))
      graphics::smoothScatter(MuBase_old, MedianTau_old, 
           bty = "n", xlab = "Mean expresssion",  
           ylab = paste("Log2 fold change", GroupLabel1, "vs", GroupLabel2), 
           main = "Before correction (log2offset = red line)", log = "x")
      abline(h = 0, lty = 2)
      abline(h = log(OffsetEst), lty = 1, col = "red")
      
      graphics::smoothScatter(MuBase, MedianTau, 
           bty = "n", xlab = "Mean expresssion", 
           ylab = paste("Log2 fold change", GroupLabel1, "vs", GroupLabel2),
           main = "After correction", log = "x")
      abline(h = 0, lty = 2)
      par(ask=FALSE)
    }
  }
  else
  {
    message(paste("--------------------------------------------------------------------- \n",
                  "It is recomended to perform a global offset correction between the two groups of cells \n",
                  "Default offset value set equal to 1.\n",
                  "To perform offset correction, please set 'Offset = TRUE'. \n",
                  "--------------------------------------------------------------------- \n"))
    
    Summary1 <- Summary(Chain1)
    Summary2 <- Summary(Chain2)
    MuBase = (Summary1@mu[,1] * n1 + Summary2@mu[,1] * n2)/n
    ChainTau = log2(Chain1@mu / Chain2@mu)
    MedianTau = matrixStats::colMedians(ChainTau) 
    # Default values when no offset correction is applied
    OffsetEst = 1; OffsetChain = NULL
    Chain1_offset = NULL; Chain2_offset = NULL

  }
  
  Search = is.null(EviThresholdM)
  
  # Changes in mean expression
  AuxMean <- HiddenThresholdSearchTestDE(ChainTau, EpsilonM, EviThresholdM,
                                         GenesSelect, EFDR_M,
                                         Task = "Differential mean")  
  ProbM = AuxMean$Prob; OptThresholdM = AuxMean$OptThreshold
  # Test results
  MeanPlus1 = which(ProbM > OptThresholdM[1] & MedianTau > 0)
  MeanPlus2 = which(ProbM > OptThresholdM[1] & MedianTau < 0)  
  ResultDiffMean = rep("NoDiff", length(MedianTau))
  ResultDiffMean[MeanPlus1] = paste0(GroupLabel1,"+")
  ResultDiffMean[MeanPlus2] = paste0(GroupLabel2,"+") 
  if(!is.null(GenesSelect)) {ResultDiffMean[!GenesSelect] = "ExcludedByUser"}
  # Output table
  TableMean = cbind.data.frame("GeneNames" = GeneNames,
                               "MeanOverall"= round(as.numeric(MuBase),3),
                               "Mean1"= round(as.numeric(Summary1@mu[,1]),3),
                               "Mean2"= round(as.numeric(Summary2@mu[,1]),3),
                               "MeanFC"= round(as.numeric(exp(MedianTau)),3),
                               "MeanLog2FC"= round(as.numeric(MedianTau),3),
                               "ProbDiffMean"= round(as.numeric(ProbM),3),
                               "ResultDiffMean" = ResultDiffMean,
                               stringsAsFactors = FALSE)
  
  # Genes for which mean doesn't change
  NotDE <- which(ResultDiffMean == "NoDiff")

  # Changes in over dispersion
  ChainOmega = log2(Chain1@delta[,NotDE] / Chain2@delta[,NotDE])
  MedianOmega = matrixStats::colMedians(ChainOmega)
  DeltaBase=(Summary1@delta[NotDE,1] * n1 + Summary2@delta[NotDE,1] * n2)/n
  AuxDisp <- HiddenThresholdSearchTestDE(ChainOmega, EpsilonD, EviThresholdD,
                                         NULL, EFDR_D,
                                         Task = "Differential dispersion")  
  ProbD = AuxDisp$Prob; OptThresholdD = AuxDisp$OptThreshold
  # Test results
  DispPlus1 = which(ProbD > OptThresholdD[1] & MedianOmega > 0)
  DispPlus2 = which(ProbD > OptThresholdD[1] & MedianOmega < 0)
  ResultDiffDisp = rep("NoDiff", length(MedianOmega))
  ResultDiffDisp[DispPlus1] = paste0(GroupLabel1,"+")
  ResultDiffDisp[DispPlus2] = paste0(GroupLabel2,"+")
  # Output table
  TableDisp = cbind.data.frame("GeneNames" = GeneNames[NotDE],
                               "MeanOverall" = round(as.numeric(MuBase[NotDE]),3),
                               "DispOverall"= round(as.numeric(DeltaBase),3),
                               "Disp1"= round(as.numeric(Summary1@delta[NotDE,1]),3),
                               "Disp2"= round(as.numeric(Summary2@delta[NotDE,1]),3),  
                               "DispFC"= round(as.numeric(exp(MedianOmega)),3), 
                               "DispLog2FC"= round(as.numeric(MedianOmega),3),
                               "ProbDiffDisp"= round(as.numeric(ProbD),3),
                               "ResultDiffDisp" = ResultDiffDisp,
                               stringsAsFactors = FALSE)
  
  # Update after removing DE genes from Disp table!
  # Reordering the tables
  GeneIndex = 1:length(MuBase)
  
  if(OrderVariable == "GeneIndex") orderVar = GeneIndex
  if(OrderVariable == "GeneNames") orderVar = GeneNames
  if(OrderVariable == "Prob") orderVar = ProbM
  TableMean = TableMean[order(orderVar, decreasing = TRUE),]
  
  if(OrderVariable == "GeneIndex") orderVar = GeneIndex[NotDE]
  if(OrderVariable == "GeneNames") orderVar = GeneNames[NotDE]
  if(OrderVariable == "Prob") orderVar = ProbD
  TableDisp = TableDisp[order(orderVar, decreasing = TRUE),]
  
  if(!is.null(GenesSelect))
  {
    message("--------------------------------------------------------------------- \n",
            paste("The user excluded ", sum(!GenesSelect), " genes from the comparison. \n"),
            "These genes are marked as 'ExcludedByUser' in the results table and excluded from EFDR calibration. \n",
            "--------------------------------------------------------------------- \n")    
  }
  
  nMeanPlus1 = sum(ResultDiffMean == paste0(GroupLabel1,"+"))
  nMeanPlus2 = sum(ResultDiffMean == paste0(GroupLabel2,"+"))
  nDispPlus1 = sum(ResultDiffDisp == paste0(GroupLabel1,"+"))
  nDispPlus2 = sum(ResultDiffDisp == paste0(GroupLabel2,"+"))
  
  message("--------------------------------------------------------------------- \n",
          paste(nMeanPlus1 + nMeanPlus2, " genes with a change on the overall expression:  \n"),
          paste("- Higher expression in ",GroupLabel1,"samples:", nMeanPlus1, "\n"),
          paste("- Higher expression in ",GroupLabel2,"samples:", nMeanPlus2, "\n"),
          paste("- Fold change tolerance = ", round(100*EpsilonM,2), "% \n"),
          paste("- Evidence threshold = ", OptThresholdM[1], "\n"),
          paste("- EFDR = ", round(100*OptThresholdM[2],2), "% \n"),
          paste("- EFNR = ", round(100*OptThresholdM[3],2), "% \n"),
          "--------------------------------------------------------------------- \n",
          "\n",
          "--------------------------------------------------------------------- \n",
          paste(nDispPlus1 + nDispPlus2, " genes with a change on the cell-to-cell biological over dispersion:  \n"),
          paste("- Higher over dispersion in ",GroupLabel1,"samples:", nDispPlus1, "\n"),
          paste("- Higher over dispersion in ",GroupLabel2,"samples:", nDispPlus2, "\n"),
          paste("- Fold change tolerance = ", round(100*EpsilonD,2), "% \n"),
          paste("- Evidence threshold = ", OptThresholdD[1], "\n"),
          paste("- EFDR = ", round(100*OptThresholdD[2],2), "% \n"),
          paste("- EFNR = ", round(100*OptThresholdD[3],2), "% \n"),
          paste("NOTE: differential dispersion results only include the", length(MedianOmega),
                "for which the mean did not change. \n"),
          "--------------------------------------------------------------------- \n")
  
  if(Plot)
  {    
    args <- list(...)
    
    par(ask=TRUE)
        
    if(Search)
    {    
        message("EFDR/EFNR control ...")
        par(mfrow = c(1,2))
        EviThresholds <- seq(0.5,0.9995,by=0.00025)
        plot(EviThresholds, AuxMean$EFDRgrid, type = "l", lty = 1, 
             bty = "n", ylab = "Error rate", xlab = "Evidence threshold", 
             ylim = c(0,1), main = "Differential mean test")
        lines(EviThresholds, AuxMean$EFNRgrid, lty = 2)  
        abline(h = EFDR_M, col = "blue", lwd = 2, lty = 1)
        abline(v = OptThresholdM[1], col = "red", lwd = 2, lty = 1)
        legend('top', c("EFDR", "EFNR", "Target EFDR"), lty = c(1:2,1), 
               col = c("black","black", "blue"), bty = "n")
        plot(EviThresholds, AuxDisp$EFDRgrid, type = "l", lty = 1, 
             bty = "n", ylab = "Error rate", xlab = "Evidence threshold", 
             ylim = c(0,1), main = "Differential dispersion test")
        lines(EviThresholds, AuxDisp$EFNRgrid, lty = 2)  
        abline(h = EFDR_D, col = "blue", lwd = 2, lty = 1)
        abline(v = OptThresholdD[1], col = "red", lwd = 2, lty = 1)
        legend('top', c("EFDR", "EFNR", "Target EFDR"), lty = c(1:2,1), 
               col = c("black","black", "blue"), bty = "n")
    }
    
    message("MA plot ...")
    
    par(mfrow = c(1,2))
    graphics::smoothScatter(MuBase, MedianTau, 
                            bty = "n", xlab = "Mean expresssion",  
                            ylab = paste("Log2 fold change", GroupLabel1, "vs", GroupLabel2), 
                            main = "Differential mean test", log = "x")
    with(TableMean[!(TableMean$ResultDiffMean %in% c("ExcludedByUser", "NoDiff")), ],
         points(MeanOverall, MeanLog2FC, pch = 16, col = "red"))
    abline(h = 0, lty = 2)
    graphics::smoothScatter(MuBase[NotDE], MedianOmega, 
                            bty = "n", xlab = "Mean expresssion", 
                            ylab = paste("Log2 fold change", GroupLabel1, "vs", GroupLabel2),
                            main = "Differential dispersion test", log = "x")
    with(TableDisp[!(TableDisp$ResultDiffDisp %in% c("ExcludedByUser", "NoDiff")), ],
         points(MeanOverall, DispLog2FC, pch = 16, col = "red"))
    abline(h = 0, lty = 2)
    par(ask=FALSE)
  }
  
  list("TableMean" = TableMean, "TableDisp" = TableDisp, 
       "DiffMeanSummary" = list("EviThreshold" = OptThresholdM[1],
                                "EFDR" = OptThresholdM[2], 
                                "EFNR" = OptThresholdM[3]),
       "DiffDispSummary" = list("EviThreshold" = OptThresholdD[1],
                                "EFDR" = OptThresholdD[2], 
                                "EFNR" = OptThresholdD[3]),
       "Chain1_offset" = Chain1_offset, "Chain2_offset" = Chain2_offset,
       "OffsetChain" = OffsetChain, "Offset" = OffsetEst)
}