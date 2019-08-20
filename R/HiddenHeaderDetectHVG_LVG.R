HiddenHeaderDetectHVG_LVG <- function(Chain,
                                      PercentileThreshold,
                                      VarThreshold,
                                      ProbThreshold,
                                      EFDR,
                                      OrderVariable,
                                      Plot)
{
  if (!is(Chain, "BASiCS_Chain")) {
    stop("'Chain' is not a BASiCS_Chain class object.")
  }
  
  # Test if the chain contains beta parameters
  if(is.null(Chain@parameters$beta) & is.null(VarThreshold)){
    stop("'Chain' was not generated using the BASiCS regression model.
         Please supply a values between 0 and 1 for 'VarThreshold'.")
  }
  
  # Add a warning that by default, the variance decomposition threshold is
  # used if the user does not supply the PercentileThreshold parameter
  if(is.null(Chain@parameters$beta) & !is.null(PercentileThreshold)){
    warning("'Chain' was not generated using the BASiCS regression model.\n", 
            "By default, variable genes are detected by testing against ",
            "a variance threshold of ", 100*VarThreshold, "%")
  }
  
  # Add a warning that it's better to use the regression trend when estimating 
  # highly variable genes 
  if(!is.null(Chain@parameters$beta) & !is.null(VarThreshold) & 
     !is.null(PercentileThreshold)){
    warning("'Chain' was generated using the BASiCS regression model.\n",
            "By default, the ", 100*PercentileThreshold, 
            " percentile of variable genes will be returned.")
  }
  
  if (!is.null(VarThreshold)){ 
    if(VarThreshold < 0 | VarThreshold > 1 | !is.finite(VarThreshold)) {
      stop("Variance contribution threshold must be in (0,1)")
    }
  }
  
  if (!is.null(PercentileThreshold)){
    if(PercentileThreshold < 0 | PercentileThreshold > 1 | 
       !is.finite(PercentileThreshold)) {
    stop("Percentile threshold must be in (0,1)")
    }
  }
  
  if (!is.logical(Plot) | length(Plot) != 1) {
    stop("Please insert `TRUE` or `FALSE` for `Plot` parameter")
  }
  if (!is.null(ProbThreshold)) {
    if (ProbThreshold < 0 | ProbThreshold > 1 | !is.finite(ProbThreshold)) {
      stop("Posterior probability threshold must be contained in (0,1) \n
           For automatic threshold search use `ProbThreshold = NULL`.")
    }
    else {
      message("Posterior probability threshold has been supplied \n",
              "EFDR will not be calibrated. \n")
    }
  }
  if (!(OrderVariable %in% c("GeneNames", "GeneIndex", "Prob")))
    stop("Invalid 'OrderVariable' value")

  if (is.null(ProbThreshold)) {
    message("Posterior probability threshold set by EFDR = ",
            100 * EFDR, "% (+-2.5% tolerance) ...")
  }
}
