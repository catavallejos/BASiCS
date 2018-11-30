HiddenHeaderDetectHVG_LVG <- function(object,
                                      VarThreshold,
                                      ProbThreshold = NULL,
                                      EFDR = 0.05,
                                      OrderVariable = "Prob",
                                      Plot = FALSE)
{
  if (!is(object, "BASiCS_Chain")) {
    stop("'object' is not a BASiCS_Chain class object.")
  }
  if (VarThreshold < 0 | VarThreshold > 1 | !is.finite(VarThreshold)) {
    stop("Variance contribution threshold must be in (0,1)")
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
  if (!(OrderVariable %in% c("GeneNames", "Mu", "Delta", "Sigma", "Prob")))
    stop("Invalid 'OrderVariable' value")

  if (is.null(ProbThreshold)) {
    message("Posterior probability threshold set by EFDR = ",
            100 * EFDR, "% (+-2.5% tolerance) ...")
  }
}
