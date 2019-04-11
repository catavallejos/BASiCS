#' Calculate effective sample size for BASiCS_Chain parameters
#' 
#' @param object A BASiCS_Chain object.
#' @param Param The parameter to use to calculate effectiveSize.
#' @param na.rm Remove NA values before calculating effectiveSize. Default is TRUE.
#' 
#' @examples
#' data(ChainSCReg)
#' # Will often fail for real chains, when epsilon is NA
#' effectiveSize(ChainSCReg@parameters[["epsilon"]])
#' # This has an na.rm argument which removes these NA values before
#' # calculating effectiveSize
#' BASiCS_effectiveSize(ChainSC, Param = "mu")
#' 
#' @export
BASiCS_effectiveSize <- function(object, Param, na.rm = TRUE) {
  HiddenGetMeasure(object, Param, Measure = "effectiveSize", na.rm)
}
