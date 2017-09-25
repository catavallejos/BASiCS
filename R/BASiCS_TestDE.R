#' @name BASiCS_TestDE
#' 
#' @title Detection of genes with changes in expression
#' 
#' @description Function to assess changes in expression between two groups
#' of cells (mean and over-dispersion)
#' 
#' @param Chain1 an object of class \code{\linkS4class{BASiCS_Chain}} 
#' containing parameter estimates for the first group of cells
#' @param Chain2 an object of class \code{\linkS4class{BASiCS_Chain}} 
#' containing parameter estimates for the second group of cells
#' @param EpsilonM Minimum fold change tolerance threshold for detecting 
#' changes in overall expression (must be a positive real number). 
#' Default value: \code{EpsilonM = log2(1.5)} (i.e. 50\% increase).
#' @param EpsilonD Minimum fold change tolerance threshold for detecting 
#' changes in biological over-dispersion (must be a positive real number). 
#' Default value: \code{EpsilonM = log2(1.5)} (i.e. 50\% increase).
#' @param ProbThresholdM Optional parameter. Probdence threshold for detecting 
#' changes in overall expression (must be a positive value, between 0 and 1)
#' @param ProbThresholdD Optional parameter. Probdence threshold for detecting 
#' changes in cell-to-cell biological over-dispersion 
#' (must be a positive value, between 0 and 1)
#' @param OrderVariable Ordering variable for output. 
#' Possible values: \code{'GeneIndex'}, \code{'GeneName'} and \code{'Prob'}.
#' @param GroupLabel1 Label assigned to reference group. 
#' Default: \code{GroupLabel1 = 'Group1'}
#' @param GroupLabel2 Label assigned to reference group. 
#' Default: \code{GroupLabel2 = 'Group2'}
#' @param Plot If \code{Plot = TRUE}, MA and volcano plots are generated. 
#' @param PlotOffset If \code{Plot = TRUE}, the offset effect is visualised.   
#' @param OffSet Optional argument to remove a fix offset effect (if not 
#' previously removed from the MCMC chains). This argument will be removed 
#' shorly, once offset removal is built as an internal step. 
#' @param EFDR_M Target for expected false discovery rate related to 
#' the comparison of means. Default \code{EFDR_M = 0.10}.
#' @param EFDR_D Target for expected false discovery rate related to 
#' the comparison of dispersions. Default \code{EFDR_D = 0.10}.
#' @param GenesSelect Optional argument to provide a user-defined list 
#' of genes to be considered for the comparison.
#' Default: \code{GenesSelect = NULL}. When used, this argument must be a vector 
#' of \code{TRUE} (include gene) / \code{FALSE} (exclude gene) indicator, 
#' with the same length as the number of intrinsic genes and following the same 
#' order as how genes are displayed in the table of counts.  
#' This argument is necessary in order to have a meaningful EFDR calibration
#' when the user decides to exclude some genes from the comparison. 
#' @param ... Graphical parameters (see \code{\link[graphics]{par}}).
#' 
#' @return \code{BASiCS_TestDE} returns a list of 4 elements:
#' \describe{
#' \item{\code{TableMean}}{A \code{\link[base]{data.frame}} containing the 
#'       results of the differential mean test
#'    \describe{
#'    \item{\code{GeneName}}{Gene name}
#'    \item{\code{MeanOverall}}{For each gene, the estimated mean expression 
#'          parameter \eqn{\mu_i} is averaged across both groups of cells 
#'          (weighted by sample size).}
#'    \item{\code{Mean1}}{Estimated mean expression parameter \eqn{\mu_i} for 
#'          each biological gene in the first group of cells.}
#'    \item{\code{Mean2}}{Estimated mean expression parameter \eqn{\mu_i} for 
#'          each biological gene in the second group of cells.}
#'    \item{\code{MeanFC}}{Fold change in mean expression parameters between 
#'          the first and second groups of cells.}
#'    \item{\code{MeanLog2FC}}{Log2-transformed fold change in mean expression
#'          between the first and second groups of cells.}
#'    \item{\code{ProbDiffMean}}{Posterior probability for mean expression 
#'          difference between the first and second groups of cells.}
#'    \item{\code{ResultDiffExp}}{Indicator if a gene has a higher mean 
#'          expression in the first or second groups of cells.}
#'    }}
#' \item{\code{TableDisp}}{A \code{\link[base]{data.frame}} containing the 
#'       results of the differential dispersion test (excludes genes for which 
#'       the mean does not changes).
#'    \describe{
#'    \item{\code{GeneName}}{Gene name}
#'    \item{\code{MeanOverall}}{For each gene, the estimated mean expression 
#'          parameter \eqn{\mu_i} is averaged across both groups of cells 
#'          (weighted by sample size).}
#'    \item{\code{DispOverall}}{For each gene, the estimated over-dispersion 
#'          parameter \eqn{\delta_i} is averaged across both groups of cells 
#'          (weighted by sample size).}
#'    \item{\code{Disp1}}{Estimated over-dispersion parameter \eqn{\delta_i} 
#'          for each biological gene in the first group of cells.}
#'    \item{\code{Disp2}}{Estimated over-dispersion parameter \eqn{\delta_i} 
#'          for each biological gene in the second group of cells.}
#'    \item{\code{DispFC}}{Fold change in over-dispersion parameters between 
#'          the between the first and second groups of cells.}
#'    \item{\code{DispLog2FC}}{Log-transformed fold change in over-dispersion 
#'          between the first and second groups of cells.}
#'    \item{\code{ProbDiffDisp}}{Posterior probability for over-dispersion 
#'          difference between the first and second groups of cells.}
#'    \item{\code{ResultDiffDisp}}{Indicator if a gene has a higher 
#'          over-dispersion in the first or second groups of cells.}
#'    }}
#'  \item{\code{DiffExpSummary}}{A list containing the following information 
#'        for the differential mean expression test:
#'    \describe{
#'   \item{\code{ProbThreshold}}{Posterior probability threshold.}
#'   \item{\code{EFDR}}{Expected false discovery rate for the given thresholds.}
#'   \item{\code{EFNR}}{Expected false negative rate for the given thresholds.}
#'   }}
#'  \item{\code{DiffOverDispSummary}}{A list containing the following 
#'        information for the differential over-dispersion test:
#'    \describe{
#'   \item{\code{ProbThreshold}}{Posterior probability threshold.}
#'   \item{\code{EFDR}}{Expected false discovery rate for the given thresholds.}
#'   \item{\code{EFNR}}{Expected false negative rate for the given thresholds.}
#'   }}
#'  \item{\code{Chain1_offset}}{an \code{\link[BASiCS]{BASiCS_Chain}} object: 
#'        \code{Chain1} after offset removal.}
#'  \item{\code{Chain2_offset}}{an \code{\link[BASiCS]{BASiCS_Chain}} object: 
#'        \code{Chain2} after offset removal (this is only provided for 
#'        completeness; \code{Chain2} is not affected by the offset).}
#'  \item{\code{OffsetChain}}{MCMC chain calculated for the offset effect.}
#'  \item{\code{Offset}}{Estimated offset (posterior median of 
#'        \code{OffsetChain}). Default value set equal to 1 when offset 
#'        correction is not performed.}
#' }
#'  
#' @examples
#' 
#' # Loading two 'BASiCS_Chain' objects (obtained using 'BASiCS_MCMC')
#' data(ChainSC)
#' data(ChainRNA)
#' 
#' Test <- BASiCS_TestDE(Chain1 = ChainSC, Chain2 = ChainRNA,
#'                       GroupLabel1 = 'SC', GroupLabel2 = 'P&S',
#'                       EpsilonM = log2(1.5), EpsilonD = log2(1.5), 
#'                       OffSet = TRUE)
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
#' @references 
#' Vallejos, Richardson and Marioni (2016). Genome Biology.
#'
#' @rdname BASiCS_TestDE
BASiCS_TestDE <- function(Chain1, 
                          Chain2, 
                          EpsilonM = log2(1.5), 
                          EpsilonD = log2(1.5), 
                          ProbThresholdM = NULL, 
                          ProbThresholdD = NULL, 
                          OrderVariable = "Prob", 
                          GroupLabel1 = "Group1", 
                          GroupLabel2 = "Group2", 
                          Plot = TRUE, 
                          PlotOffset = TRUE, 
                          OffSet = TRUE, 
                          EFDR_M = 0.1, 
                          EFDR_D = 0.1, 
                          GenesSelect = NULL, ...) 
{
  # Checking validity of input arguments
  if (!is(Chain1, "BASiCS_Chain")) 
    stop("'Chain1' is not a 'BASiCS_Chain' class object.")
  if (!is(Chain2, "BASiCS_Chain")) 
    stop("'Chain2' is not a 'BASiCS_Chain' class object.")
    
  # Test compatibility of both BASiCS_Chain objects
  if (ncol(Chain1@mu) != ncol(Chain2@mu)) 
    stop("The 'BASiCS_Chain' objects contain different number of genes.")
  if (!identical(colnames(Chain1@mu), colnames(Chain2@mu))) 
    stop("The  'BASiCS_Chain' objects contain genes in different order.")
    
  GeneName <- colnames(Chain1@mu)
    
  if (EpsilonM < 0 | !is.finite(EpsilonM)) 
    stop("'EpsilonM' must be a positive real value")
  if (EpsilonD < 0 | !is.finite(EpsilonD)) 
    stop("'EpsilonD' must be a positive real value")
  if (!is.logical(Plot) | length(Plot) != 1) 
    stop("Please insert TRUE or FALSE for 'Plot' parameter")
  if (!is.null(ProbThresholdM) | !is.null(ProbThresholdD)) 
  {
    if (ProbThresholdM < 0 | ProbThresholdM > 1 | !is.finite(ProbThresholdM)) 
      stop("'ProbThresholdM' must be contained in (0,1) \n 
            For automatic threshold search use ProbThresholdM = NULL.")
    if (ProbThresholdD < 0 | ProbThresholdD > 1 | !is.finite(ProbThresholdD)) 
      stop("'ProbThresholdD' must be contained in (0,1) \n 
            For automatic threshold search use ProbThresholdD = NULL.")
    }
    if (!(OrderVariable %in% c("GeneIndex", "GeneName", "Prob"))) 
      stop("Invalid 'OrderVariable' value")
    if (!is.character(GroupLabel1) | length(GroupLabel1) > 1) 
      stop("Invalid value for 'GroupLabel1'")
    if (!is.character(GroupLabel2) | length(GroupLabel2) > 1) 
      stop("Invalid value for 'GroupLabel2'")
    if (!is.null(GenesSelect) & (length(GenesSelect) != length(GeneName))) 
      stop("Invalid value for 'GenesSelect'")
    if (!is.null(GenesSelect) & !is.logical(GenesSelect)) 
      stop("Invalid value for 'GenesSelect'")

    message("Log-fold change thresholds are now set in a log2 scale. \n", 
            "Original BASiCS release used a natural logarithm scale.")
    
    n1 <- ncol(Chain1@nu)
    n2 <- ncol(Chain2@nu)
    n <- n1 + n2
    # With offset correction
    if (OffSet) 
    {
      # Calculating iteration-specific offset
      OffsetChain <- matrixStats::rowSums2(Chain1@mu) / 
                      matrixStats::rowSums2(Chain2@mu)
      # Offset point estimate
      OffsetEst <- median(OffsetChain)
        
      # Offset correction
      Chain1_offset <- Chain1
      Chain1_offset@mu <- Chain1@mu / OffsetEst
      Chain1_offset@phi <- Chain1@phi * OffsetEst
      Chain2_offset <- Chain2  # Chain2 requires no change
      Summary1 <- Summary(Chain1_offset)
      Summary2 <- Summary(Chain2_offset)
        
      # Pre-offset correction LFC estimates
      Summary1_old <- Summary(Chain1)
      Summary2_old <- Summary(Chain2)
      MuBase_old <- (Summary1_old@mu[, 1] * n1 + Summary2_old@mu[, 1] * n2)/n
      ChainTau_old <- log2(Chain1@mu / Chain2@mu)
      MedianTau_old <- matrixStats::colMedians(ChainTau_old)
        
      # Offset corrected LFC estimates
      MuBase <- (Summary1@mu[, 1] * n1 + Summary2@mu[, 1] * n2)/n
      ChainTau <- log2(Chain1_offset@mu / Chain2_offset@mu)
      MedianTau <- matrixStats::colMedians(ChainTau)
        
      if (!PlotOffset) 
      {
        message("-------------------------------------------------------------\n", 
                "Offset estimate:", round(OffsetEst, 4), "\n",  
                "(ratio", GroupLabel1, "vs", GroupLabel2, ").\n", 
                "To visualise its effect, please use 'PlotOffset = TRUE'.\n", 
                "-------------------------------------------------------------\n")
      } 
      else 
      {
        message("-------------------------------------------------------------\n", 
                "Offset estimate:", round(OffsetEst, 4), "\n", 
                "(ratio", GroupLabel1, "vs", GroupLabel2, ").\n", 
                "-------------------------------------------------------------\n")
      }

      if (PlotOffset == TRUE) 
      {
        message("Plots to follow: \n",
                "1. Posterior uncertainty for offset estimate \n",
                "2. Mean expression estimates before/after offset correction \n",
                "3. MA plot before/after offset correction \n")
        
        par(ask = TRUE)
        # Offset uncertainty
        graphics::boxplot(OffsetChain, frame = FALSE, 
                          main = "Offset MCMC chain", ylab = "Offset estimate")
        # Mean expression parameters before/after offset correction
        par(mfrow = c(1, 2))
        graphics::boxplot(cbind(Summary1_old@mu[,1], Summary2_old@mu[,1]), 
                          frame = FALSE, main = "Before correction", 
                          names = c(GroupLabel1, GroupLabel2), 
                          ylab = "Mean expression", log = "y")
        graphics::boxplot(cbind(Summary1@mu[, 1], Summary2@mu[, 1]), 
                          frame = FALSE, main = "After correction", 
                          names = c(GroupLabel1, GroupLabel2), 
                          ylab = "Mean expression", log = "y")
        # MA plot pre/after offset
        par(mfrow = c(1, 2))
        graphics::smoothScatter(log2(MuBase_old), MedianTau_old, bty = "n", 
                                xlab = "Mean expresssion (log2)", 
                                ylab = paste("Log2 fold change", GroupLabel1, 
                                             "vs", GroupLabel2), 
                                main = "Before correction")
        abline(h = 0, lty = 2)
        abline(h = log2(OffsetEst), lty = 1, col = "red")
        legend('topright', "log2offset", lty = 1, col = "red")
        graphics::smoothScatter(log2(MuBase), MedianTau, bty = "n", 
                                xlab = "Mean expresssion (log2)", 
                                ylab = paste("Log2 fold change", GroupLabel1, 
                                             "vs", GroupLabel2), 
                                main = "After correction")
        abline(h = 0, lty = 2)
        par(ask = FALSE)
      }
    } 
    else 
    {
      message("-------------------------------------------------------------\n", 
              "It is recomended to perform a global offset correction \n",
              "to remove global changes between the two groups of cells \n", 
              "Default offset value set equal to 1.\n", 
              "To perform offset correction, please set 'Offset = TRUE'. \n", 
              "-------------------------------------------------------------\n")
      
      Summary1 <- Summary(Chain1)
      Summary2 <- Summary(Chain2)
      MuBase <- (Summary1@mu[, 1] * n1 + Summary2@mu[, 1] * n2)/n
      ChainTau <- log2(Chain1@mu / Chain2@mu)
      MedianTau <- matrixStats::colMedians(ChainTau)
      
      # Default values when no offset correction is applied
      OffsetEst <- 1
      OffsetChain <- NULL
      Chain1_offset <- NULL
      Chain2_offset <- NULL
        
    }
    
    Search <- is.null(ProbThresholdM)
    
    # Changes in mean expression
    AuxMean <- HiddenThresholdSearchTestDE(ChainTau, EpsilonM, 
                                           ProbThresholdM, 
                                           GenesSelect, EFDR_M, 
                                           Task = "differential mean")
    ProbM <- AuxMean$Prob
    OptThresholdM <- AuxMean$OptThreshold
    
    # Test results
    MeanPlus1 <- which(ProbM > OptThresholdM[1] & MedianTau > 0)
    MeanPlus2 <- which(ProbM > OptThresholdM[1] & MedianTau < 0)
    ResultDiffMean <- rep("NoDiff", length(MedianTau))
    ResultDiffMean[MeanPlus1] <- paste0(GroupLabel1, "+")
    ResultDiffMean[MeanPlus2] <- paste0(GroupLabel2, "+")
    if (!is.null(GenesSelect)) {
        ResultDiffMean[!GenesSelect] <- "ExcludedByUser"
    }
    # Output table
    TableMean <- cbind.data.frame(GeneName = GeneName, 
                                  MeanOverall = as.numeric(MuBase), 
                                  Mean1 = as.numeric(Summary1@mu[,1]), 
                                  Mean2 = as.numeric(Summary2@mu[,1]), 
                                  MeanFC = as.numeric(2^(MedianTau)), 
                                  MeanLog2FC = as.numeric(MedianTau), 
                                  ProbDiffMean = as.numeric(ProbM), 
                                  ResultDiffMean = ResultDiffMean, 
                                  stringsAsFactors = FALSE)
    # Rounding to 3 decimal points
    TableMean[, seq(2,7)] <- round(TableMean[, seq(2,7)], 3)
    
    # Genes for which mean doesn't change
    NotDE <- which(ResultDiffMean == "NoDiff")
    
    # Changes in over dispersion
    ChainOmega <- log2(Chain1@delta[, NotDE] / Chain2@delta[, NotDE])
    MedianOmega <- matrixStats::colMedians(ChainOmega)
    DeltaBase <- (Summary1@delta[NotDE,1] * n1 + Summary2@delta[NotDE,1] * n2)/n
    
    AuxDisp <- HiddenThresholdSearchTestDE(ChainOmega, EpsilonD, 
                                           ProbThresholdD, NULL, 
                                           EFDR_D, 
                                           Task = "Differential dispersion")
    ProbD <- AuxDisp$Prob
    OptThresholdD <- AuxDisp$OptThreshold
    
    # Test results
    DispPlus1 <- which(ProbD > OptThresholdD[1] & MedianOmega > 0)
    DispPlus2 <- which(ProbD > OptThresholdD[1] & MedianOmega < 0)
    ResultDiffDisp <- rep("NoDiff", length(MedianOmega))
    ResultDiffDisp[DispPlus1] <- paste0(GroupLabel1, "+")
    ResultDiffDisp[DispPlus2] <- paste0(GroupLabel2, "+")
    # Output table
    TableDisp <- cbind.data.frame(GeneName = GeneName[NotDE], 
                                  MeanOverall = as.numeric(MuBase[NotDE]), 
                                  DispOverall = as.numeric(DeltaBase), 
                                  Disp1 = as.numeric(Summary1@delta[NotDE, 1]), 
                                  Disp2 = as.numeric(Summary2@delta[NotDE, 1]), 
                                  DispFC = as.numeric(2^(MedianOmega)), 
                                  DispLog2FC = as.numeric(MedianOmega), 
                                  ProbDiffDisp = as.numeric(ProbD), 
                                  ResultDiffDisp = ResultDiffDisp, 
                                  stringsAsFactors = FALSE)
    # Rounding to 3 decimal points
    TableDisp[, seq(2,8)] <- round(TableDisp[, seq(2,8)], 3)
    
    # Update after removing DE genes from Disp table!  Reordering the tables
    GeneIndex <- seq_len(length(MuBase))
    
    if (OrderVariable == "GeneIndex") 
        orderVar = GeneIndex
    if (OrderVariable == "GeneName") 
        orderVar = GeneName
    if (OrderVariable == "Prob") 
        orderVar = ProbM
    TableMean <- TableMean[order(orderVar, decreasing = TRUE), ]
    
    if (OrderVariable == "GeneIndex") 
        orderVar = GeneIndex[NotDE]
    if (OrderVariable == "GeneName") 
        orderVar = GeneName[NotDE]
    if (OrderVariable == "Prob") 
        orderVar = ProbD
    TableDisp <- TableDisp[order(orderVar, decreasing = TRUE), ]
    
    if (!is.null(GenesSelect)) 
    {
      message("-------------------------------------------------------------\n", 
              "The user excluded ", sum(!GenesSelect), " genes. \n", 
              "These genes are marked as 'ExcludedByUser' \n", 
              "and excluded from EFDR calibration. \n", 
              "-------------------------------------------------------------\n")
    }
    
    
    if (Plot) 
    {
      args <- list(...)
        
      if (Search) 
      {
        message("Plots to follow: \n",
                "1. EFDR/EFNR control plots \n",
                "2. MA plots \n",
                "3. Volcano plots \n")
      } 
      else 
      {
        message("Plots to follow: \n",
                "1. MA plots \n",
                "2. Volcano plots \n")
      }
        
        par(ask = TRUE)
        
        if (Search) 
        {
          par(mfrow = c(1, 2))
          ProbThresholds <- seq(0.5, 0.9995, by = 0.00025)
          plot(ProbThresholds, AuxMean$EFDRgrid, type = "l", lty = 1, bty = "n", 
               ylab = "Error rate", xlab = "Probdence threshold", 
               ylim = c(0, 1), main = "Differential mean")
          lines(ProbThresholds, AuxMean$EFNRgrid, lty = 2)
          abline(h = EFDR_M, col = "blue", lwd = 2, lty = 1)
          abline(v = OptThresholdM[1], col = "red", lwd = 2, lty = 1)
          legend("top", c("EFDR", "EFNR", "Target EFDR"), lty = c(1:2, 1), 
                 col = c("black", "black", "blue"), bty = "n")
          plot(ProbThresholds, AuxDisp$EFDRgrid, type = "l", lty = 1, bty = "n", 
               ylab = "Error rate", xlab = "Probdence threshold", 
               ylim = c(0, 1), main = "Differential dispersion")
          lines(ProbThresholds, AuxDisp$EFNRgrid, lty = 2)
          abline(h = EFDR_D, col = "blue", lwd = 2, lty = 1)
          abline(v = OptThresholdD[1], col = "red", lwd = 2, lty = 1)
          legend("top", c("EFDR", "EFNR", "Target EFDR"), lty = c(1:2, 1), 
                 col = c("black", "black", "blue"), bty = "n")
        }
        
        par(mfrow = c(1, 2))
        with(TableMean, 
             graphics::smoothScatter(log2(MeanOverall), MeanLog2FC, 
                                     bty = "n", 
                                     xlab = "Mean expresssion (log2)", 
                                     ylab = paste("Log2 fold change", 
                                                  GroupLabel1, "vs", 
                                                  GroupLabel2), 
                                     main = "Differential mean"))
        with(TableMean[!(TableMean$ResultDiffMean %in% 
                           c("ExcludedByUser", "NoDiff")), ], 
             points(log2(MeanOverall), MeanLog2FC, pch = 16, col = "red"))
        abline(h = c(-EpsilonM, EpsilonM), lty = 2)
        with(TableDisp, 
             graphics::smoothScatter(log2(MeanOverall), DispLog2FC, 
                                     bty = "n", 
                                     xlab = "Mean expresssion (log2)", 
                                     ylab = paste("Log2 fold change", 
                                                  GroupLabel1, "vs", 
                                                  GroupLabel2), 
                                     main = "Differential dispersion"))
        with(TableDisp[TableDisp$ResultDiffDisp != "NoDiff", ],
             points(log2(MeanOverall), DispLog2FC, pch = 16, col = "red"))
        abline(h = c(-EpsilonD, EpsilonD), lty = 2)
        
        par(mfrow = c(1, 2))
        with(TableMean, 
             graphics::smoothScatter(MeanLog2FC, ProbDiffMean, 
                                     bty = "n", ylim = c(0, 1), 
                                     ylab = "Posterior probability", 
                                     xlab = paste("Log2 fold change", 
                                                  GroupLabel1, "vs", 
                                                  GroupLabel2), 
                                     main = "Differential mean test"))
        with(TableMean[!(TableMean$ResultDiffMean %in% 
                           c("ExcludedByUser", "NoDiff")), ], 
             points(MeanLog2FC, ProbDiffMean, pch = 16, col = "red"))
        abline(v = c(-EpsilonM, EpsilonM), lty = 2)
        with(TableDisp, 
             graphics::smoothScatter(DispLog2FC, ProbDiffDisp, 
                                     bty = "n", ylim = c(0, 1), 
                                     ylab = "Posterior probability", 
                                     xlab = paste("Log2 fold change", 
                                                  GroupLabel1, "vs", 
                                                  GroupLabel2), 
                                     main = "Differential dispersion test"))
        with(TableDisp[TableDisp$ResultDiffDisp != "NoDiff", ], 
             points(DispLog2FC, ProbDiffDisp, pch = 16, col = "red"))
        abline(v = c(-EpsilonD, EpsilonD), lty = 2)
        
        par(ask = FALSE)
    }
    
    nMeanPlus1 <- sum(ResultDiffMean == paste0(GroupLabel1, "+"))
    nMeanPlus2 <- sum(ResultDiffMean == paste0(GroupLabel2, "+"))
    nDispPlus1 <- sum(ResultDiffDisp == paste0(GroupLabel1, "+"))
    nDispPlus2 <- sum(ResultDiffDisp == paste0(GroupLabel2, "+"))
    
    message("-------------------------------------------------------------\n", 
            nMeanPlus1 + nMeanPlus2," genes with a change in mean expression:\n", 
            "- Higher expression in ", GroupLabel1, "samples:", nMeanPlus1, "\n", 
            "- Higher expression in ", GroupLabel2, "samples:", nMeanPlus2, "\n", 
            "- Fold change tolerance = ", round(2^(EpsilonM), 2), "% \n",  
            "- Probability threshold = ", OptThresholdM[1], "\n", 
            "- EFDR = ", round(100 * OptThresholdM[2], 2), "% \n", 
            "- EFNR = ", round(100 * OptThresholdM[3], 2), "% \n", 
            "-------------------------------------------------------------\n\n", 
            "-------------------------------------------------------------\n", 
            nDispPlus1 + nDispPlus2," genes with a change in over dispersion:\n", 
            "- Higher dispersion in ", GroupLabel1, "samples:", nDispPlus1,"\n", 
            "- Higher dispersion in ", GroupLabel2, "samples:", nDispPlus2,"\n", 
            "- Fold change tolerance = ", round(2^(EpsilonD), 2), "% \n", 
            "- Probability threshold = ", OptThresholdD[1], "\n", 
            "- EFDR = ", round(100 * OptThresholdD[2], 2), "% \n", 
            "- EFNR = ", round(100 * OptThresholdD[3], 2), "% \n", 
            "NOTE: differential dispersion assessment only applied to the \n", 
            length(MedianOmega), " genes for which the mean did not change. \n", 
            "--------------------------------------------------------------\n")
    
    list(TableMean = TableMean, 
         TableDisp = TableDisp, 
         DiffMeanSummary = list(ProbThreshold = OptThresholdM[1], 
                                EFDR = OptThresholdM[2], 
                                EFNR = OptThresholdM[3]), 
         DiffDispSummary = list(ProbThreshold = OptThresholdD[1], 
                                EFDR = OptThresholdD[2], 
                                EFNR = OptThresholdD[3]), 
         Chain1_offset = Chain1_offset, 
         Chain2_offset = Chain2_offset,
         OffsetChain = OffsetChain, 
         Offset = OffsetEst)
}
