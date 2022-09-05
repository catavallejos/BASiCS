.VG <- function(x) {
  x@Table[[x@Name]]
}

.HeaderDetectHVG_LVG <- function(Chain,
                                 PercentileThreshold,
                                 VarThreshold,
                                 EpsilonThreshold,
                                 Plot) {
  
  if (!is(Chain, "BASiCS_Chain")) 
    stop("'Chain' is not a BASiCS_Chain class object.")
  
  if(is.null(PercentileThreshold) 
    & is.null(VarThreshold) 
    & is.null(EpsilonThreshold)) {
    stop("A value must be provided for 'PercentileThreshold', 'VarThreshold' or 'EpsilonThreshold'")
  }
  
  # Test if the chain does not contain epsilon parameters
  if (is.null(Chain@parameters$epsilon)) {
    
    if(!is.null(PercentileThreshold)) {
      stop("'Chain' does not include residual over-dispersion parameters.
           'PercentileThreshold' will be ignored.
           'VarThreshold' must be provided instead.")
    } 
    if(!is.null(EpsilonThreshold)) {
      stop("'Chain' does not include residual over-dispersion parameters.
           'EpsilonThreshold' will be ignored.
           'VarThreshold' must be provided instead.")
    }

    if (!is.null(VarThreshold)) {
      if (VarThreshold < 0 | VarThreshold > 1 | !is.finite(VarThreshold))
        stop("Variance contribution threshold must be in (0,1)")
    }
  }
  
  # Test if the chain contains beta parameters
  if (!is.null(Chain@parameters$epsilon)) {
    
    if(!is.null(VarThreshold)) 
      stop("'Chain' includes residual over-dispersion parameters.\n
           'VarThreshold' will be ignored. \n
           'PercentileThreshold'must be provided instead. \n")
    
    if(!is.null(PercentileThreshold)) {
      if(PercentileThreshold < 0 | PercentileThreshold > 1 | 
         !is.finite(PercentileThreshold)) 
        stop("Percentile threshold must be in (0,1)")
    }
  }
  
  if(!is.logical(Plot)) 
    stop("`Plot` must be TRUE or FALSE")
}
