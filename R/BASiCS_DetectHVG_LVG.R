#' @name BASiCS_DetectHVG
#' @aliases BASiCS_DetectHVG BASiCS_DetectHVG_LVG
#'
#' @title Detection method for highly and lowly variable genes
#'
#' @description Functions to detect highly and lowly variable genes
#'
#' @param Chain an object of class \code{\link[BASiCS]{BASiCS_Chain-class}}
#' @param VarThreshold Variance contribution threshold (must be a positive value, between 0 and 1)
#' @param ProbThreshold Optional parameter. Posterior probability threshold (must be a positive value, between 0 and 1)
#' @param EFDR Target for expected false discovery rate related to HVG/LVG detection (default = 0.10)
#' @param OrderVariable Ordering variable for output. Possible values: \code{"GeneIndex"}, \code{"Mu"},
#'  \code{"Delta"}, \code{"Sigma"} and \code{"Prob"}.
#' @param Plot If \code{Plot = TRUE} error control and expression versus HVG/LVG probability plots are generated
#' @param ... Graphical parameters (see \code{\link[graphics]{par}}).
#'
#' @return \code{BASiCS_DetectHVG} returns a list of 4 elements:
#' \describe{
#' \item{\code{Table}}{Matrix whose columns contain}
#'    \describe{
#'    \item{\code{GeneIndex}}{Vector of length \code{q.bio}. Gene index as in the order present in the analysed \code{\link[SingleCellExperiment]{SingleCellExperiment}}}
#'    \item{\code{GeneNames}}{Vector of length \code{q.bio}. Gene name as in the order present in the analysed \code{\link[SingleCellExperiment]{SingleCellExperiment}}}
#'    \item{\code{Mu}}{Vector of length \code{q.bio}. For each biological gene, posterior median of gene-specific expression levels \eqn{\mu_i}}
#'    \item{\code{Delta}}{Vector of length \code{q.bio}. For each biological gene, posterior median of gene-specific biological over-dispersion parameter \eqn{\delta_i}}
#'    \item{\code{Sigma}}{Vector of length \code{q.bio}. For each biological gene, proportion of the total variability that is due to a cell-to-cell biological heterogeneity component. }
#'    \item{\code{Prob}}{Vector of length \code{q.bio}. For each biological gene, probability of being highly variable according to the given thresholds.}
#'    \item{\code{HVG}}{Vector of length \code{q.bio}. For each biological gene, indicator of being detected as highly variable according to the given thresholds. }
#'    }
#' \item{\code{ProbThreshold}}{Posterior probability threshold.}
#' \item{\code{EFDR}}{Expected false discovery rate for the given thresholds.}
#' \item{\code{EFNR}}{Expected false negative rate for the given thresholds.}
#' }
#' \code{BASiCS_DetectLVG} produces a similar output, replacing the column \code{HVG} by \code{LVG}, an indicator of a gene being detected as lowly variable according to the given thresholds.
#'
#' @examples
#'
#' # See
#' help(BASiCS_MCMC)
#'
#' @details See vignette
#'
#'
#' @seealso \code{\link[BASiCS]{BASiCS_Chain-class}}
#'
#' @author Catalina A. Vallejos \email{cnvallej@@uc.cl}
#'
#' @references 
#' Vallejos, Marioni and Richardson (2015). Bayesian Analysis of Single-Cell Sequencing data. PLoS Computational Biology. 
#'
#' @rdname BASiCS_DetectHVG_LVG
BASiCS_DetectHVG <- function(Chain,
                             VarThreshold,
                             ProbThreshold = NULL,
                             EFDR = 0.10, 
                             OrderVariable = "Prob",
                             Plot = FALSE,
                             ...)
{
  # Safety checks
  HiddenHeaderDetectHVG_LVG(Chain, VarThreshold,
                            ProbThreshold, EFDR, OrderVariable, Plot)
    
  Search = FALSE
  if(is.null(ProbThreshold)) Search = TRUE
  
  # Variance decomposition
  VarDecomp <- HiddenVarDecomp(Chain)
  
  # HVG probability for a given variance threshold
  Prob <- HiddenProbHVG(VarThreshold = VarThreshold, VarDecomp = VarDecomp)
  
  # Threshold search
  Aux <- HiddenThresholdSearchDetectHVG_LVG(ProbThreshold, VarThreshold, Prob, EFDR)
  EFDRgrid <- Aux$EFDRgrid 
  EFNRgrid <- Aux$EFNRgrid
  ProbThresholds <- Aux$ProbThresholds
  OptThreshold <- Aux$OptThreshold
  
  # Output preparation
  Sigma <- apply(VarDecomp$BioVarGlobal, 2, median)
  Mu <- apply(Chain@mu, 2, median)
  Delta <- apply(Chain@delta, 2, median)
  HVG <- ifelse(Prob > OptThreshold[1], TRUE, FALSE)

  # Output preparation
  qbio = length(Sigma)
  Genes = 1:qbio
  GeneNames = colnames(Chain@mu)
  
  GeneIndex = 1:length(Mu)
  Table = cbind.data.frame("GeneIndex" = Genes,
                           "GeneNames" = GeneNames,
                           "Mu" = Mu,
                           "Delta" = Delta,
                           "Sigma" = Sigma,
                           "Prob" = Prob,
                           "HVG" = HVG, stringsAsFactors = FALSE)
  rownames(Table) = Genes
  if(OrderVariable == "GeneNames") orderVar = GeneNames
  if(OrderVariable == "Mu") orderVar = Mu
  if(OrderVariable == "Delta") orderVar = Delta
  if(OrderVariable == "Sigma") orderVar = Sigma
  if(OrderVariable == "Prob") orderVar = Prob
  Table = Table[order(orderVar, decreasing = TRUE),]
  
  if(Plot)
  {
    args <- list(...)
    
    if(Search)
    {
      par(ask=TRUE)
      # EFDR / EFNR plot
      HiddenPlot1DetectHVG_LVG(ProbThresholds, EFDRgrid, EFNRgrid, EFDR)
    }
    
    # Output plot : mean vs prob
    HiddenPlot2DetectHVG_LVG(args, Task = "HVG", 
                             Mu, Prob, OptThreshold, Hits = HVG)
    
    par(ask=FALSE)
  }
  
  message(paste(sum(HVG), " genes classified as highly variable using: \n"),
          paste("- Variance contribution threshold = ", round(100*VarThreshold,2), "% \n"),
          paste("- Evidence threshold = ", OptThreshold[1], "\n"),
          paste("- EFDR = ", round(100*OptThreshold[2],2), "% \n"),
          paste("- EFNR = ", round(100*OptThreshold[3],2), "% \n"))
  
  
  
  list("Table" = Table,
       "EviThreshold" = OptThreshold[1], "EFDR" = OptThreshold[2], "EFNR" = OptThreshold[3])
}

#' @name BASiCS_DetectLVG
#' @aliases BASiCS_DetectLVG BASiCS_DetectHVG_LVG
#' @rdname BASiCS_DetectHVG_LVG
BASiCS_DetectLVG <- function(Chain,
                             VarThreshold,
                             ProbThreshold = NULL,
                             EFDR = 0.05,
                             OrderVariable = "Prob",
                             Plot = FALSE,
                             ...)
{
  # Safety checks
  HiddenHeaderDetectHVG_LVG(Chain, VarThreshold,
                            ProbThreshold, EFDR, OrderVariable, Plot)
  
  Search = FALSE
  if(is.null(ProbThreshold)) Search = TRUE
  
  # Variance decomposition
  VarDecomp <- HiddenVarDecomp(Chain)
  
  # LVG probability for a given variance threshold
  Prob <- HiddenProbLVG(VarThreshold = VarThreshold, VarDecomp = VarDecomp)
  
  # Threshold search
  Aux <- HiddenThresholdSearchDetectHVG_LVG(ProbThreshold, VarThreshold, Prob, EFDR)
  EFDRgrid <- Aux$EFDRgrid 
  EFNRgrid <- Aux$EFNRgrid
  ProbThresholds <- Aux$ProbThresholds
  OptThreshold <- Aux$OptThreshold
  
  # Output preparation
  Sigma <- apply(VarDecomp$BioVarGlobal, 2, median)
  Mu <- apply(Chain@mu, 2, median)
  Delta <- apply(Chain@delta, 2, median)
  LVG <- ifelse(Prob > OptThreshold[1], TRUE, FALSE)

  qbio = length(Sigma)
  Genes = 1:qbio
  GeneNames = colnames(Chain@mu)
  GeneIndex = 1:length(Mu)
  Table = cbind.data.frame("GeneIndex" = Genes,
                           "GeneNames" = GeneNames,
                           "Mu" = Mu,
                           "Delta" = Delta,
                           "Sigma" = Sigma,
                           "Prob" = Prob,
                           "LVG" = LVG, stringsAsFactors = FALSE)
  rownames(Table) = Genes
  
  if(OrderVariable == "GeneNames") orderVar = GeneNames
  if(OrderVariable == "Mu") orderVar = Mu
  if(OrderVariable == "Delta") orderVar = Delta
  if(OrderVariable == "Sigma") orderVar = Sigma
  if(OrderVariable == "Prob") orderVar = Prob
  Table = Table[order(orderVar, decreasing = TRUE),]
  
  if(Plot)
  {
    args <- list(...)
    
    if(Search)
    {
      par(ask=TRUE)
      # EFDR / EFNR plot
      HiddenPlot1DetectHVG_LVG(ProbThresholds, EFDRgrid, EFNRgrid, EFDR)
    }
    
    # Output plot : mean vs prob
    HiddenPlot2DetectHVG_LVG(args, Task = "LVG", 
                             Mu, Prob, OptThreshold, Hits = LVG)
    
    par(ask=FALSE)
    
  }
  
  message(paste(sum(LVG), " genes classified as lowly variable using: \n"),
          paste("- Variance contribution threshold = ", round(100*VarThreshold,2), "% \n"),
          paste("- Evidence threshold = ", OptThreshold[1], "\n"),
          paste("- EFDR = ", round(100*OptThreshold[2],2), "% \n"),
          paste("- EFNR = ", round(100*OptThreshold[3],2), "% \n"))
  
  list("Table" = Table,
       "EviThreshold" = OptThreshold[1], "EFDR" = OptThreshold[2], "EFNR" = OptThreshold[3])
}
