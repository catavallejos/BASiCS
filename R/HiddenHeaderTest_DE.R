HiddenHeaderTest_DE <- function(Chain1,
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
                                Offset,
                                ...)
{
  # Checking validity of input arguments
  if (!is(Chain1, "BASiCS_Chain"))
    stop("'Chain1' is not a 'BASiCS_Chain' class object.")
  if (!is(Chain2, "BASiCS_Chain"))
    stop("'Chain2' is not a 'BASiCS_Chain' class object.")

  # Test compatibility of both BASiCS_Chain objects
  if (ncol(Chain1@parameters$mu) != ncol(Chain2@parameters$mu))
    stop("The 'BASiCS_Chain' objects contain different number of genes.")
  if (!identical(colnames(Chain1@parameters$mu), colnames(Chain2@parameters$mu)))
    stop("The  'BASiCS_Chain' objects contain genes in different order.")

  if (EpsilonM < 0 | !is.finite(EpsilonM))
    stop("'EpsilonM' must be a positive real value")
  if (EpsilonD < 0 | !is.finite(EpsilonD))
    stop("'EpsilonD' must be a positive real value")

  if (!is.logical(Plot) | length(Plot) != 1)
    stop("Please insert TRUE or FALSE for 'Plot' parameter")
  if (!is.logical(PlotOffset) | length(PlotOffset) != 1)
    stop("Please insert TRUE or FALSE for 'PlotOffset' parameter")
  if (!is.logical(Offset) | length(Offset) != 1)
    stop("Please insert TRUE or FALSE for 'Offset' parameter")

  if("OffSet" %in% names(list(...))) {
    stop("'OffSet' is no longer a valid argument. Use 'Offset' instead.\n")
  }

  if (!is.null(ProbThresholdM))
  {
    if (ProbThresholdM < 0 | ProbThresholdM > 1 | !is.finite(ProbThresholdM))
      stop("'ProbThresholdM' must be contained in (0,1) \n
           For automatic threshold search use ProbThresholdM = NULL.")
  }
  if (!is.null(ProbThresholdD))
  {
    if (ProbThresholdD < 0 | ProbThresholdD > 1 | !is.finite(ProbThresholdD))
      stop("'ProbThresholdD' must be contained in (0,1) \n
           For automatic threshold search use ProbThresholdD = NULL.")
  }

  if (!is.null(EFDR_M) | !is.null(EFDR_D))
  {
    if (EFDR_M < 0 | EFDR_M > 1 | !is.finite(EFDR_M))
      stop("'EFDR_M' must be contained in (0,1). \n")
    if (EFDR_D < 0 | EFDR_D > 1 | !is.finite(EFDR_D))
      stop("'EFDR_D' must be contained in (0,1). \n")
  }

  if (!(OrderVariable %in% c("GeneIndex", "GeneName", "Prob")))
    stop("Invalid 'OrderVariable' value")
  if (!is.character(GroupLabel1) | length(GroupLabel1) > 1)
    stop("Invalid value for 'GroupLabel1'")
  if (!is.character(GroupLabel2) | length(GroupLabel2) > 1)
    stop("Invalid value for 'GroupLabel2'")

  GeneName <- colnames(Chain1@parameters$mu)
  if (!is.null(GenesSelect) & (length(GenesSelect) != length(GeneName)))
    stop("Invalid value for 'GenesSelect'")
  if (!is.null(GenesSelect) & !is.logical(GenesSelect))
    stop("Invalid value for 'GenesSelect'")

  if(!is.null(ProbThresholdM))
    message("A value has been provided for `ProbThresholdM` \n",
            "EFDR will not be calibrated in differential mean test. \n")

  if(!is.null(ProbThresholdD))
    message("A value has been provided for `ProbThresholdD` \n",
            "EFDR will not be calibrated in differential over-dispersion test. \n")
}
