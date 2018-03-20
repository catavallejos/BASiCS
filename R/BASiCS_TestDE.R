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
#' @param EpsilonR Minimum distance threshold for detecting 
#' changes in residual over-dispersion (must be a positive real number). 
#' Default value: \code{EpsilonR= log2(1.5)/log2(exp(1))} (i.e. 50\% increase).
#' @param ProbThresholdM Optional parameter. Probdence threshold for detecting 
#' changes in overall expression (must be a positive value, between 0 and 1)
#' @param ProbThresholdD Optional parameter. Probdence threshold for detecting 
#' changes in cell-to-cell biological over-dispersion 
#' (must be a positive value, between 0 and 1)
#' @param ProbThresholdR Optional parameter. Probdence threshold for detecting 
#' changes in residual over-dispersion 
#' (must be a positive value, between 0 and 1)
#' @param OrderVariable Ordering variable for output. 
#' Possible values: \code{'GeneIndex'}, \code{'GeneName'} and \code{'Prob'}.
#' @param GroupLabel1 Label assigned to reference group. 
#' Default: \code{GroupLabel1 = 'Group1'}
#' @param GroupLabel2 Label assigned to reference group. 
#' Default: \code{GroupLabel2 = 'Group2'}
#' @param Plot If \code{Plot = TRUE}, MA and volcano plots are generated. 
#' @param PlotOffset If \code{Plot = TRUE}, the offset effect is visualised.   
#' @param Offset Optional argument to remove a fix offset effect (if not 
#' previously removed from the MCMC chains). Default: \code{Offset = TRUE}. 
#' @param EFDR_M Target for expected false discovery rate related to 
#' the comparison of means. Default \code{EFDR_M = 0.10}.
#' @param EFDR_D Target for expected false discovery rate related to 
#' the comparison of dispersions. Default \code{EFDR_D = 0.10}.
#' @param EFDR_R Target for expected false discovery rate related to 
#' the comparison of residual over-dispersions. Default \code{EFDR_D = 0.10}.
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
#' \item{\code{TableResDisp}}{A \code{\link[base]{data.frame}} containing the 
#'       results of the differential residual over-dispersion test 
#'       (excludes genes that are not expressed in at least 2 cells).
#'    \describe{
#'    \item{\code{GeneName}}{Gene name}
#'    \item{\code{MeanOverall}}{For each gene, the estimated mean expression 
#'          parameter \eqn{\mu_i} is averaged across both groups of cells 
#'          (weighted by sample size).}
#'    \item{\code{ResDispOverall}}{For each gene, the estimated residual 
#'          over-dispersion parameter \eqn{\delta_i} is averaged across both 
#'          groups of cells (weighted by sample size).}
#'    \item{\code{ResDisp1}}{Estimated residual over-dispersion parameter 
#'          \eqn{\epsilon_i} for each biological gene in the first group of 
#'          cells.}
#'    \item{\code{ResDisp2}}{Estimated residual over-dispersion parameter 
#'          \eqn{\epsilon_i} for each biological gene in the second group of 
#'          cells.}
#'    \item{\code{ResDispDistance}}{Difference in residual over-dispersion 
#'          between the first and second groups of cells.}
#'    \item{\code{ProbDiffResDisp}}{Posterior probability for residual 
#'          over-dispersion difference between the first and second groups of 
#'          cells.}
#'    \item{\code{ResultDiffResDisp}}{Indicator if a gene has a higher 
#'          residual over-dispersion in the first or second groups of cells.}
#'    }}
#'  \item{\code{DiffMeanSummary}}{A list containing the following information 
#'        for the differential mean expression test:
#'    \describe{
#'   \item{\code{ProbThreshold}}{Posterior probability threshold.}
#'   \item{\code{EFDR}}{Expected false discovery rate for the given thresholds.}
#'   \item{\code{EFNR}}{Expected false negative rate for the given thresholds.}
#'   }}
#'  \item{\code{DiffDispSummary}}{A list containing the following 
#'        information for the differential over-dispersion test:
#'    \describe{
#'   \item{\code{ProbThreshold}}{Posterior probability threshold.}
#'   \item{\code{EFDR}}{Expected false discovery rate for the given thresholds.}
#'   \item{\code{EFNR}}{Expected false negative rate for the given thresholds.}
#'   }}
#'  \item{\code{DiffResDispSummary}}{A list containing the following 
#'        information for the differential residual over-dispersion test:
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
#' # For testing differences in residual over-dispersion, two chains obtained 
#' # via 'BASiCS_MCMC(Data, N, Thin, Burn, Regression=TRUE)' need to be provided
#' data(ChainSCReg)
#' data(ChainRNAReg)
#' 
#' Test <- BASiCS_TestDE(Chain1 = ChainSCReg, Chain2 = ChainRNAReg,
#'                       GroupLabel1 = 'SC', GroupLabel2 = 'P&S',
#'                       EpsilonM = log2(1.5), EpsilonD = log2(1.5), 
#'                       EpsilonR = log2(1.5)/log2(exp(1)),
#'                       OffSet = TRUE)        
#' 
#' @author Catalina A. Vallejos \email{cnvallej@@uc.cl} 
#' @author Nils Eling \email{eling@@ebi.ac.uk}
#'
#' @rdname BASiCS_TestDE
BASiCS_TestDE <- function(Chain1, 
                          Chain2, 
                          EpsilonM = log2(1.5), 
                          EpsilonD = log2(1.5), 
                          EpsilonR = log2(1.5)/log2(exp(1)),
                          ProbThresholdM = NULL, 
                          ProbThresholdD = NULL, 
                          ProbThresholdR = NULL,
                          OrderVariable = "GeneIndex", 
                          GroupLabel1 = "Group1", 
                          GroupLabel2 = "Group2", 
                          Plot = TRUE, 
                          PlotOffset = TRUE, 
                          Offset = TRUE, 
                          EFDR_M = 0.1, 
                          EFDR_D = 0.1, 
                          EFDR_R = 0.1, 
                          GenesSelect = NULL, ...) 
{
  
  HiddenHeaderTest_DE(Chain1, 
                      Chain2,
                      EpsilonM, 
                      EpsilonD,
                      EFDR_M,
                      EFDR_D, 
                      ProbThresholdM, 
                      ProbThresholdD,
                      OrderVariable,
                      GroupLabel1,
                      GroupLabel2,
                      GenesSelect,
                      Plot, 
                      PlotOffset,
                      Offset)
  
  GeneName <- colnames(Chain1@parameters$mu)
  GeneIndex <- seq_len(length(GeneName))
  
  message("-------------------------------------------------------------\n",
          "Log-fold change thresholds are now set in a log2 scale. \n", 
          "Original BASiCS release used a natural logarithm scale.")
  
  n1 <- ncol(Chain1@parameters$nu)
  n2 <- ncol(Chain2@parameters$nu)
  n <- n1 + n2
  # With offset correction
  if (Offset) 
  {
    # Calculating iteration-specific offset
    OffsetChain <- matrixStats::rowSums2(Chain1@parameters$mu) / 
                    matrixStats::rowSums2(Chain2@parameters$mu)
    # Offset point estimate
    OffsetEst <- median(OffsetChain)
        
    # Offset correction
    Chain1_offset <- Chain1
    Chain1_offset@parameters$mu <- Chain1@parameters$mu / OffsetEst
    Chain1_offset@parameters$phi <- Chain1@parameters$phi * OffsetEst
    Chain2_offset <- Chain2  # Chain2 requires no change
    Summary1 <- Summary(Chain1_offset)
    Summary2 <- Summary(Chain2_offset)
        
    # Pre-offset correction LFC estimates
    Summary1_old <- Summary(Chain1)
    Summary2_old <- Summary(Chain2)
    MuBase_old <- (Summary1_old@parameters$mu[, 1] * n1 + Summary2_old@parameters$mu[, 1] * n2)/n
    ChainTau_old <- log2(Chain1@parameters$mu / Chain2@parameters$mu)
    MedianTau_old <- matrixStats::colMedians(ChainTau_old)

    # Offset corrected LFC estimates
    MuBase <- (Summary1@parameters$mu[, 1] * n1 + Summary2@parameters$mu[, 1] * n2)/n
    ChainTau <- log2(Chain1_offset@parameters$mu / Chain2_offset@parameters$mu)
    MedianTau <- matrixStats::colMedians(ChainTau)
        
    if (!PlotOffset) 
    {
      message("-------------------------------------------------------------\n", 
              "Offset estimate: ", round(OffsetEst, 4), "\n",  
              "(ratio ", GroupLabel1, " vs ", GroupLabel2, ").\n", 
              "To visualise its effect, please use 'PlotOffset = TRUE'.\n", 
              "-------------------------------------------------------------\n")
    } 
    else 
    {
      message("-------------------------------------------------------------\n", 
              "Offset estimate: ", round(OffsetEst, 4), "\n", 
              "(ratio ", GroupLabel1, " vs ", GroupLabel2, ").\n", 
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
      graphics::boxplot(cbind(Summary1_old@parameters$mu[,1], Summary2_old@parameters$mu[,1]), 
                        frame = FALSE, main = "Before correction", 
                        names = c(GroupLabel1, GroupLabel2), 
                        ylab = "Mean expression", log = "y")
      graphics::boxplot(cbind(Summary1@parameters$mu[, 1], Summary2@parameters$mu[, 1]), 
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
    MuBase <- (Summary1@parameters$mu[, 1] * n1 + Summary2@parameters$mu[, 1] * n2)/n
    ChainTau <- log2(Chain1@parameters$mu / Chain2@parameters$mu)
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
                                         Task = "Differential mean")
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
                                Mean1 = as.numeric(Summary1@parameters$mu[,1]), 
                                Mean2 = as.numeric(Summary2@parameters$mu[,1]), 
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
  ChainOmega <- log2(Chain1@parameters$delta[, NotDE] / Chain2@parameters$delta[, NotDE])
  MedianOmega <- matrixStats::colMedians(ChainOmega)
  DeltaBase <- (Summary1@parameters$delta[NotDE,1] * n1 + Summary2@parameters$delta[NotDE,1] * n2)/n

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
                                Disp1 = as.numeric(Summary1@parameters$delta[NotDE, 1]), 
                                Disp2 = as.numeric(Summary2@parameters$delta[NotDE, 1]), 
                                DispFC = as.numeric(2^(MedianOmega)), 
                                DispLog2FC = as.numeric(MedianOmega), 
                                ProbDiffDisp = as.numeric(ProbD), 
                                ResultDiffDisp = ResultDiffDisp, 
                                stringsAsFactors = FALSE)
  # Rounding to 3 decimal points
  TableDisp[, seq(2,8)] <- round(TableDisp[, seq(2,8)], 3)
  
  # Changes in residual over-dispersion - if regression approach was used
  if(!is.null(Chain1@parameters$epsilon)){
    NotExcluded <- !(is.na(Chain1@parameters$epsilon[1,]) | 
                               is.na(Chain2@parameters$epsilon[1,]))
      
    ChainPsi <- Chain1@parameters$epsilon[,NotExcluded] - 
                    Chain2@parameters$epsilon[,NotExcluded]
    MedianPsi <- matrixStats::colMedians(ChainPsi)
    EpsilonBase <- (Summary1@parameters$epsilon[NotExcluded,1] * n1 + 
                      Summary2@parameters$epsilon[NotExcluded,1] * n2)/n
    
    AuxResDisp <- HiddenThresholdSearchTestDE(ChainPsi, EpsilonR, 
                                              ProbThresholdR, 
                                              GenesSelect[NotExcluded], 
                                              EFDR_R, 
                                              Task = "Differential residual dispersion")
    ProbE <- AuxResDisp$Prob
    OptThresholdE <- AuxResDisp$OptThreshold
    
    # Test results
    ResDispPlus1 <- which(ProbE > OptThresholdE[1] & MedianPsi > 0)
    ResDispPlus2 <- which(ProbE > OptThresholdE[1] & MedianPsi < 0)
    ResultDiffResDisp <- rep("NoDiff", length(MedianPsi))
    ResultDiffResDisp[ResDispPlus1] <- paste0(GroupLabel1, "+")
    ResultDiffResDisp[ResDispPlus2] <- paste0(GroupLabel2, "+")
    
    if (!is.null(GenesSelect)) {
      cur_GenesSelect <- GenesSelect[NotExcluded]
      ResultDiffResDisp[!cur_GenesSelect] <- "ExcludedByUser"
    }
    else{
      cur_GenesSelect <- NotExcluded
    }
    
    # Output table
    TableResDisp <- cbind.data.frame(GeneName = GeneName[NotExcluded], 
      MeanOverall = as.numeric(MuBase[NotExcluded]), 
      ResDispOverall = as.numeric(EpsilonBase), 
      ResDisp1 = as.numeric(Summary1@parameters$epsilon[NotExcluded, 1]), 
      ResDisp2 = as.numeric(Summary2@parameters$epsilon[NotExcluded, 1]),  
      ResDispDistance = as.numeric(MedianPsi), 
      ProbDiffResDisp = as.numeric(ProbE), 
      ResultDiffResDisp = ResultDiffResDisp, 
      stringsAsFactors = FALSE)
  
    # Rounding to 3 decimal points
    TableResDisp[, seq(2,7)] <- round(TableResDisp[, seq(2,7)], 3)
      
    if (OrderVariable == "GeneIndex") 
      orderVar <- order(GeneIndex[NotExcluded], decreasing = FALSE)
    if (OrderVariable == "GeneName") 
      orderVar <- order(GeneName[NotExcluded], decreasing = TRUE)
    if (OrderVariable == "Prob") 
      orderVar <- order(ProbE[NotExcluded], decreasing = TRUE)
    TableResDisp <- TableResDisp[orderVar, ]
  }
  
  if (OrderVariable == "GeneIndex") 
      orderVar <- order(GeneIndex, decreasing = FALSE)
  if (OrderVariable == "GeneName") 
    orderVar <- order(GeneName, decreasing = TRUE)
  if (OrderVariable == "Prob") 
    orderVar <- order(ProbM, decreasing = TRUE)
  TableMean <- TableMean[orderVar, ]
    
  if (OrderVariable == "GeneIndex") 
    orderVar <- order(GeneIndex[NotDE], decreasing = FALSE)
  if (OrderVariable == "GeneName") 
    orderVar <- order(GeneName[NotDE], decreasing = TRUE)
  if (OrderVariable == "Prob") 
    orderVar <- order(ProbD[NotDE], decreasing = TRUE)
  TableDisp <- TableDisp[orderVar, ]
    
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
      legend("top", c("EFDR", "EFNR", "Target EFDR"), lty = c(1, 2, 1), 
             col = c("black", "black", "blue"), bty = "n")
      plot(ProbThresholds, AuxDisp$EFDRgrid, type = "l", lty = 1, bty = "n", 
           ylab = "Error rate", xlab = "Probdence threshold", 
           ylim = c(0, 1), main = "Differential dispersion")
      lines(ProbThresholds, AuxDisp$EFNRgrid, lty = 2)
      abline(h = EFDR_D, col = "blue", lwd = 2, lty = 1)
      abline(v = OptThresholdD[1], col = "red", lwd = 2, lty = 1)
      legend("top", c("EFDR", "EFNR", "Target EFDR"), lty = c(1, 2, 1), 
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
    if(!is.null(Chain1@parameters$epsilon)){
      nResDispPlus1 <- sum(ResultDiffResDisp == paste0(GroupLabel1, "+"))
      nResDispPlus2 <- sum(ResultDiffResDisp == paste0(GroupLabel2, "+"))
    }
    
    if(!is.null(Chain1@parameters$epsilon)){
      message("-------------------------------------------------------------\n", 
            nMeanPlus1 + nMeanPlus2," genes with a change in mean expression:\n", 
            "- Higher expression in ", GroupLabel1, " samples: ", nMeanPlus1, "\n", 
            "- Higher expression in ", GroupLabel2, " samples: ", nMeanPlus2, "\n", 
            "- Fold change tolerance = ", round(2^(EpsilonM)*100, 2), "% \n",  
            "- Probability threshold = ", OptThresholdM[1], "\n", 
            "- EFDR = ", round(100 * OptThresholdM[2], 2), "% \n", 
            "- EFNR = ", round(100 * OptThresholdM[3], 2), "% \n", 
            "-------------------------------------------------------------\n\n", 
            "-------------------------------------------------------------\n", 
            nDispPlus1 + nDispPlus2," genes with a change in over dispersion:\n", 
            "- Higher dispersion in ", GroupLabel1, " samples: ", nDispPlus1,"\n", 
            "- Higher dispersion in ", GroupLabel2, " samples: ", nDispPlus2,"\n", 
            "- Fold change tolerance = ", round(2^(EpsilonD)*100, 2), "% \n", 
            "- Probability threshold = ", OptThresholdD[1], "\n", 
            "- EFDR = ", round(100 * OptThresholdD[2], 2), "% \n", 
            "- EFNR = ", round(100 * OptThresholdD[3], 2), "% \n", 
            "NOTE: differential dispersion assessment only applied to the \n", 
            length(MedianOmega), " genes for which the mean did not change. \n", 
            "--------------------------------------------------------------\n",
            "-------------------------------------------------------------\n", 
            nResDispPlus1 + nResDispPlus2," genes with a change in residual over dispersion:\n", 
            "- Higher residual dispersion in ", GroupLabel1, " samples: ", nResDispPlus1,"\n", 
            "- Higher residual dispersion in ", GroupLabel2, " samples: ", nResDispPlus2,"\n", 
            "- Distance tolerance = ", round(EpsilonR, 2), "\n", 
            "- Probability threshold = ", OptThresholdE[1], "\n", 
            "- EFDR = ", round(100 * OptThresholdE[2], 2), "% \n", 
            "- EFNR = ", round(100 * OptThresholdE[3], 2), "% \n", 
            "NOTE: differential residual dispersion assessment applied to \n", 
            sum(cur_GenesSelect), " genes expressed in at least 2 cells per condition \n",
            "and included for testing. \n",
            "--------------------------------------------------------------\n")
    
      list(TableMean = TableMean, 
         TableDisp = TableDisp, 
         TableResDisp = TableResDisp,
         DiffMeanSummary = list(ProbThreshold = OptThresholdM[1], 
                                EFDR = OptThresholdM[2], 
                                EFNR = OptThresholdM[3]), 
         DiffDispSummary = list(ProbThreshold = OptThresholdD[1], 
                                EFDR = OptThresholdD[2], 
                                EFNR = OptThresholdD[3]), 
         DiffResDispSummary = list(ProbThreshold = OptThresholdE[1], 
                                EFDR = OptThresholdE[2], 
                                EFNR = OptThresholdE[3]), 
         Chain1_offset = Chain1_offset, 
         Chain2_offset = Chain2_offset,
         OffsetChain = OffsetChain, 
         Offset = OffsetEst)
    }
    else{
      message("-------------------------------------------------------------\n", 
              nMeanPlus1 + nMeanPlus2," genes with a change in mean expression:\n", 
              "- Higher expression in ", GroupLabel1, " samples: ", nMeanPlus1, "\n", 
              "- Higher expression in ", GroupLabel2, " samples: ", nMeanPlus2, "\n", 
              "- Fold change tolerance = ", round(2^(EpsilonM)*100, 2), "% \n",  
              "- Probability threshold = ", OptThresholdM[1], "\n", 
              "- EFDR = ", round(100 * OptThresholdM[2], 2), "% \n", 
              "- EFNR = ", round(100 * OptThresholdM[3], 2), "% \n", 
              "-------------------------------------------------------------\n\n", 
              "-------------------------------------------------------------\n", 
              nDispPlus1 + nDispPlus2," genes with a change in over dispersion:\n", 
              "- Higher dispersion in ", GroupLabel1, " samples: ", nDispPlus1,"\n", 
              "- Higher dispersion in ", GroupLabel2, " samples: ", nDispPlus2,"\n", 
              "- Fold change tolerance = ", round(2^(EpsilonD)*100, 2), "% \n", 
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

}
