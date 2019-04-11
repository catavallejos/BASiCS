#' Calculate effective sample size for BASiCS_Chain parameters
#' 
#' @param object an object of class \code{\linkS4class{BASiCS_Chain}}.
#' @param Param The parameter to use to calculate effectiveSize. Possible
#' values: \code{'mu'}, \code{'delta'}, \code{'phi'}, \code{'s'}, 
#' \code{'nu'}, \code{'theta'}, \code{'beta'}, \code{'sigma2'} and 
#' \code{'epsilon'}. 
#' @param na.rm Remove \code{NA} values before calculating effectiveSize. 
#' Only relevant when \code{Param = "epsilon"} (genes with very low 
#' expression are excluding when infering the mean/over-dispersion trend. 
#' Default: \code{na.rm = TRUE}.
#' 
#' @examples
#' data(ChainSCReg)
#' # Will often fail for real chains, when epsilon is NA
#' # Only relevant when \code{Param = "epsilon"} (genes with very low 
#' # expression are excluding when infering the mean/over-dispersion trend.
#' effectiveSize(ChainSCReg@parameters[["epsilon"]])
#'
#' # This has an na.rm argument which removes these NA values before
#' # calculating effectiveSize
#' BASiCS_effectiveSize(ChainSC, Param = "mu")
#' 
#' @export
BASiCS_effectiveSize <- function(object, Param, na.rm = TRUE) {
  HiddenGetMeasure(object, Param, Measure = "effectiveSize", na.rm)
}
