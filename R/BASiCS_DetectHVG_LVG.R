#' @name BASiCS_DetectHVG
#' @aliases BASiCS_DetectHVG BASiCS_DetectHVG_LVG
#'
#' @title Detection method for highly (HVG) and lowly (LVG) variable genes
#'
#' @description Functions to detect highly and lowly variable genes
#'
#' @param Chain an object of class \code{\linkS4class{BASiCS_Chain}}
#' @param PercentileThreshold Threshold to detect a percentile of variable genes
#' (must be a positive value, between 0 and 1). Defaults: 0.9 for HVG, 0.1 for LVG
#' @param VarThreshold Variance contribution threshold
#' (must be a positive value, between 0 and 1). This is only used when the 
#' BASiCS non-regression model was used to generate the Chain object.
#' @param ProbThreshold Optional parameter. Posterior probability threshold
#' (must be a positive value, between 0 and 1)
#' @param EFDR Target for expected false discovery rate related
#' to HVG/LVG detection (default = 0.10)
#' @param OrderVariable Ordering variable for output.
#' Possible values: \code{'GeneIndex'}, \code{'GeneName'} and \code{'Prob'}.
#' @param Plot If \code{Plot = TRUE} error control and
#' expression versus HVG/LVG probability plots are generated
#' @param ... Graphical parameters (see \code{\link[graphics]{par}}).
#'
#' @return \code{BASiCS_DetectHVG} returns a list of 4 elements:
#' \describe{
#' \item{\code{Table}}{Matrix whose columns can contain}
#'    \describe{
#'    \item{\code{GeneIndex}}{Vector of length \code{q.bio}.
#'         Gene index as in the order present in the analysed
#'         \code{\link[SingleCellExperiment]{SingleCellExperiment}}}
#'    \item{\code{GeneName}}{Vector of length \code{q.bio}.
#'                Gene name as in the order present in the analysed
#'          \code{\link[SingleCellExperiment]{SingleCellExperiment}}}
#'    \item{\code{Mu}}{Vector of length \code{q.bio}. For each biological gene,
#'          posterior median of gene-specific mean expression
#'          parameters \eqn{\mu_i}}
#'    \item{\code{Delta}}{Vector of length \code{q.bio}. For each biological
#'          gene, posterior median of gene-specific biological
#'          over-dispersion parameter \eqn{\delta_i}}
#'    \item{\code{Sigma}}{Vector of length \code{q.bio}.
#'          For each biological gene, proportion of the total variability
#'          that is due to a biological heterogeneity component. }
#'    \item{\code{Epsilon}}{Vector of length \code{q.bio}.
#'          For each biological gene, posterior median of gene-specific residual
#'          over-dispersion parameter \eqn{\epsilon_i}. }
#'    \item{\code{Prob}}{Vector of length \code{q.bio}.
#'          For each biological gene, probability of being highly variable
#'          according to the given thresholds.}
#'    \item{\code{HVG}}{Vector of length \code{q.bio}.
#'          For each biological gene, indicator of being detected as highly
#'          variable according to the given thresholds. }
#'    \item{\code{LVG}}{Vector of length \code{q.bio}.
#'          For each biological gene, indicator of being detected as lowly
#'          variable according to the given thresholds. }
#'    }
#' \item{\code{ProbThreshold}}{Posterior probability threshold.}
#' \item{\code{EFDR}}{Expected false discovery rate for the given thresholds.}
#' \item{\code{EFNR}}{Expected false negative rate for the given thresholds.}
#' }
#'
#' @examples
#'
#' # See
#' help(BASiCS_MCMC)
#'
#' @details See vignette
#'
#'
#' @seealso \code{\linkS4class{BASiCS_Chain}}
#'
#' @author Catalina A. Vallejos \email{cnvallej@@uc.cl}
#' @author Nils Eling \email{eling@@ebi.ac.uk}
#'
#' @references
#'
#' Vallejos, Marioni and Richardson (2015). PLoS Computational Biology.
#'
#' @rdname BASiCS_DetectHVG_LVG
#' @export
BASiCS_DetectHVG <- function(Chain,
                             PercentileThreshold =  0.9,
                             VarThreshold = NULL,
                             ProbThreshold = NULL,
                             EFDR = 0.1,
                             OrderVariable = "Prob",
                             Plot = FALSE, ...)
{
  # Safety checks
  HiddenHeaderDetectHVG_LVG(Chain, PercentileThreshold, VarThreshold,
                            ProbThreshold, EFDR, OrderVariable, Plot)

  Search <- ifelse(is.null(ProbThreshold), TRUE, FALSE)
  
  # If the BASiCS_Chain object was generated using the regression approach,
  # BASiCS finds the top highly variable genes based on the posteriors of the 
  # epsilon parameters. Otherwise, the old approach is used, which initially 
  # performs the variance decomposition.
  
  if(!is.null(Chain@parameters$beta) & !is.null(PercentileThreshold)){
    # Find the epsilon threshold that correspond to the top 'PercentileThreshold'
    # genes
    
    nGenes <- ncol(Chain@parameters$epsilon)
    
    Epsilon <- matrixStats::colMedians(Chain@parameters$epsilon)
    EpsilonThreshold <- quantile(Epsilon, PercentileThreshold, na.rm = TRUE)
    
    # HVG probability for a given epsilon threshold
    Prob <- matrixStats::colMeans2(ifelse(Chain@parameters$epsilon >
                                     EpsilonThreshold, 1, 0))
    
    # Threshold search
    Aux <- HiddenThresholdSearchDetectHVG_LVG(ProbThreshold, Prob[!is.na(Prob)], EFDR)
    
    if(Search)
    {
      EFDRgrid <- Aux$EFDRgrid
      EFNRgrid <- Aux$EFNRgrid
      ProbThresholds <- Aux$ProbThresholds
    }
    OptThreshold <- Aux$OptThreshold
    
    # Output preparation
    Mu <- matrixStats::colMedians(Chain@parameters$mu)
    HVG <- ifelse(Prob > OptThreshold[1], TRUE, FALSE)
    
    GeneIndex <- seq_along(Mu)
    GeneName <- colnames(Chain@parameters$mu)
    
    Table <- cbind.data.frame(GeneIndex = GeneIndex, GeneName = GeneName,
                              Mu = Mu, Epsilon = Epsilon,
                              Prob = Prob,
                              HVG = HVG, stringsAsFactors = FALSE)
    
    # Re-order the table of results
    if (OrderVariable == "GeneIndex") { orderVar <- GeneIndex }
    if (OrderVariable == "GeneName") { orderVar <- GeneName }
    if (OrderVariable == "Prob") { orderVar <- Prob }
    Table <- Table[order(orderVar, decreasing = TRUE, na.last = TRUE), ]
    
    if (Plot)
    {
      if (Search)
      {
        # EFDR / EFNR plot
        par(ask = TRUE)
        HiddenPlot1DetectHVG_LVG(ProbThresholds, EFDRgrid, EFNRgrid, EFDR)
      }
      
      # Output plot : mean vs prob
      HiddenPlot2DetectHVG_LVG(Task = "HVG", Mu, Prob,
                               OptThreshold, Hits = HVG, ...)
      
      par(ask = FALSE)
    }
    
    message(sum(HVG, na.rm = TRUE), " genes classified as highly variable using: \n",
            "- The ",
            round(100 * PercentileThreshold, 2), "percentile of variable genes \n",
            "- Evidence threshold = ", OptThreshold[1], "\n",
            "- EFDR = ", round(100 * OptThreshold[2], 2), "% \n",
            "- EFNR = ", round(100 * OptThreshold[3], 2), "% \n")
    
    list(Table = Table, EviThreshold = OptThreshold[1],
         EFDR = OptThreshold[2], EFNR = OptThreshold[3])
  }
  else{
    # Variance decomposition
    VarDecomp <- HiddenVarDecomp(Chain)

    # HVG probability for a given variance threshold
    Prob <- matrixStats::colMeans2(ifelse(VarDecomp$BioVarGlobal >
                                            VarThreshold, 1, 0))

    # Threshold search
    Aux <- HiddenThresholdSearchDetectHVG_LVG(ProbThreshold, Prob, EFDR)
    if(Search)
    {
      EFDRgrid <- Aux$EFDRgrid
      EFNRgrid <- Aux$EFNRgrid
      ProbThresholds <- Aux$ProbThresholds
    }
    OptThreshold <- Aux$OptThreshold

    # Output preparation
    Sigma <- matrixStats::colMedians(VarDecomp$BioVarGlobal)
    Mu <- matrixStats::colMedians(Chain@parameters$mu)
    Delta <- matrixStats::colMedians(Chain@parameters$delta)
    HVG <- ifelse(Prob > OptThreshold[1], TRUE, FALSE)

    GeneIndex <- seq_along(Mu)
    GeneName <- colnames(Chain@parameters$mu)

    Table <- cbind.data.frame(GeneIndex = GeneIndex, GeneName = GeneName,
                              Mu = Mu, Delta = Delta,
                              Sigma = Sigma, Prob = Prob,
                              HVG = HVG, stringsAsFactors = FALSE)

    # Re-order the table of results
    if (OrderVariable == "GeneName") { orderVar <- GeneName }
    if (OrderVariable == "Mu") { orderVar <- Mu }
    if (OrderVariable == "Delta") { orderVar <- Delta }
    if (OrderVariable == "Sigma") { orderVar <- Sigma }
    if (OrderVariable == "Prob") { orderVar <- Prob }
    Table <- Table[order(orderVar, decreasing = TRUE), ]

    if (Plot)
    {
      if (Search)
      {
        # EFDR / EFNR plot
        par(ask = TRUE)
        HiddenPlot1DetectHVG_LVG(ProbThresholds, EFDRgrid, EFNRgrid, EFDR)
      }

      # Output plot : mean vs prob
      HiddenPlot2DetectHVG_LVG(Task = "HVG", Mu, Prob,
                               OptThreshold, Hits = HVG, ...)

      par(ask = FALSE)
    }

    message(sum(HVG), " genes classified as highly variable using: \n",
            "- Variance contribution threshold = ",
            round(100 * VarThreshold, 2), "% \n",
            "- Evidence threshold = ", OptThreshold[1], "\n",
            "- EFDR = ", round(100 * OptThreshold[2], 2), "% \n",
            "- EFNR = ", round(100 * OptThreshold[3], 2), "% \n")

    list(Table = Table, EviThreshold = OptThreshold[1],
         EFDR = OptThreshold[2], EFNR = OptThreshold[3])
  }
}

#' @name BASiCS_DetectLVG
#' @aliases BASiCS_DetectLVG BASiCS_DetectHVG_LVG
#' @rdname BASiCS_DetectHVG_LVG
#' @export
BASiCS_DetectLVG <- function(Chain,
                             PercentileThreshold = 0.1,
                             VarThreshold = NULL,
                             ProbThreshold = NULL,
                             EFDR = 0.1,
                             OrderVariable = "Prob",
                             Plot = FALSE, ...)
{
  HiddenHeaderDetectHVG_LVG(Chain, PercentileThreshold, VarThreshold,
                            ProbThreshold, EFDR, OrderVariable, Plot)
  
  Search <- ifelse(is.null(ProbThreshold), TRUE, FALSE)

  if(!is.null(Chain@parameters$beta) & !is.null(PercentileThreshold)){
    # Find the epsilon threshold that correspond to the 'PercentileThreshold'
    
    nGenes <- ncol(Chain@parameters$epsilon)
    
    Epsilon <- matrixStats::colMedians(Chain@parameters$epsilon)
    EpsilonThreshold <- quantile(Epsilon, PercentileThreshold, na.rm = TRUE)
    
    # HVG probability for a given epsilon threshold
    Prob <- matrixStats::colMeans2(ifelse(Chain@parameters$epsilon <
                                            EpsilonThreshold, 1, 0))
    
    # Threshold search
    Aux <- HiddenThresholdSearchDetectHVG_LVG(ProbThreshold, Prob[!is.na(Prob)], EFDR)
    
    if(Search)
    {
      EFDRgrid <- Aux$EFDRgrid
      EFNRgrid <- Aux$EFNRgrid
      ProbThresholds <- Aux$ProbThresholds
    }
    OptThreshold <- Aux$OptThreshold
    
    # Output preparation
    Mu <- matrixStats::colMedians(Chain@parameters$mu)
    LVG <- ifelse(Prob > OptThreshold[1], TRUE, FALSE)
    
    GeneIndex <- seq_along(Mu)
    GeneName <- colnames(Chain@parameters$mu)
    
    Table <- cbind.data.frame(GeneIndex = GeneIndex, GeneName = GeneName,
                              Mu = Mu, Epsilon = Epsilon,
                              Prob = Prob,
                              LVG = LVG, stringsAsFactors = FALSE)
    
    # Re-order the table of results
    if (OrderVariable == "GeneIndex") { orderVar <- GeneIndex }
    if (OrderVariable == "GeneName") { orderVar <- GeneName }
    if (OrderVariable == "Prob") { orderVar <- Prob }
    Table <- Table[order(orderVar, decreasing = TRUE, na.last = TRUE), ]
    
    if (Plot)
    {
      if (Search)
      {
        # EFDR / EFNR plot
        par(ask = TRUE)
        HiddenPlot1DetectHVG_LVG(ProbThresholds, EFDRgrid, EFNRgrid, EFDR)
      }
      
      # Output plot : mean vs prob
      HiddenPlot2DetectHVG_LVG(Task = "LVG", Mu, Prob,
                               OptThreshold, Hits = LVG, ...)
      
      par(ask = FALSE)
    }
    
    message(sum(LVG, na.rm = TRUE), " genes classified as lowly variable using: \n",
            "- The ",
            round(100 * PercentileThreshold, 2), "percentile of variable genes \n",
            "- Evidence threshold = ", OptThreshold[1], "\n",
            "- EFDR = ", round(100 * OptThreshold[2], 2), "% \n",
            "- EFNR = ", round(100 * OptThreshold[3], 2), "% \n")
    
    list(Table = Table, EviThreshold = OptThreshold[1],
         EFDR = OptThreshold[2], EFNR = OptThreshold[3])
  }
  else{
    # Variance decomposition
    VarDecomp <- HiddenVarDecomp(Chain)

    # LVG probability for a given variance threshold
    Prob <- matrixStats::colMeans2(ifelse(VarDecomp$BioVarGlobal <
                                            VarThreshold, 1, 0))

    # Threshold search
    Aux <- HiddenThresholdSearchDetectHVG_LVG(ProbThreshold, Prob, EFDR)
    if(Search)
    {
      EFDRgrid <- Aux$EFDRgrid
      EFNRgrid <- Aux$EFNRgrid
      ProbThresholds <- Aux$ProbThresholds
    }
    OptThreshold <- Aux$OptThreshold

    # Output preparation
    Sigma <- matrixStats::colMedians(VarDecomp$BioVarGlobal)
    Mu <- matrixStats::colMedians(Chain@parameters$mu)
    Delta <- matrixStats::colMedians(Chain@parameters$delta)
    LVG <- ifelse(Prob > OptThreshold[1], TRUE, FALSE)

    GeneIndex <- seq_along(Mu)
    GeneName <- colnames(Chain@parameters$mu)

    Table <- cbind.data.frame(GeneIndex = GeneIndex, GeneName = GeneName,
                              Mu = Mu, Delta = Delta,
                              Sigma = Sigma, Prob = Prob,
                              LVG = LVG, stringsAsFactors = FALSE)

    # Re-order the table of results
    if (OrderVariable == "GeneName") { orderVar <- GeneName }
    if (OrderVariable == "Mu") { orderVar <- Mu }
    if (OrderVariable == "Delta") { orderVar <- Delta }
    if (OrderVariable == "Sigma") { orderVar <- Sigma }
    if (OrderVariable == "Prob") { orderVar <- Prob }
    Table <- Table[order(orderVar, decreasing = TRUE), ]

    if (Plot)
    {
      if (Search)
      {
        # EFDR / EFNR plot
        par(ask = TRUE)
        HiddenPlot1DetectHVG_LVG(ProbThresholds, EFDRgrid, EFNRgrid, EFDR)
      }

      # Output plot : mean vs prob
      HiddenPlot2DetectHVG_LVG(Task = "LVG", Mu, Prob,
                               OptThreshold, Hits = LVG, ...)

      par(ask = FALSE)
    }

    message(sum(LVG), " genes classified as lowly variable using: \n",
            "- Variance contribution threshold = ",
            round(100 * VarThreshold, 2), "% \n",
            "- Evidence threshold = ", OptThreshold[1], "\n",
            "- EFDR = ", round(100 * OptThreshold[2], 2), "% \n",
            "- EFNR = ", round(100 * OptThreshold[3], 2), "% \n")

      list(Table = Table, EviThreshold = OptThreshold[1],
           EFDR = OptThreshold[2], EFNR = OptThreshold[3])
  }
}
