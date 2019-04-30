#' Calculate effective sample size for BASiCS_Chain parameters
#' 
#' @description A wrapper of \code{coda::effectiveSize} to be used with 
#' \code{\linkS4class{BASiCS_Chain}} objects.
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
#' @return A vector with effective sample sizes for all the elements 
#' of \code{Param}
#' 
#' @examples
#' 
#' data(ChainSC)
#' BASiCS_EffectiveSize(ChainSC, Param = "mu")
#' 
#' @export
BASiCS_EffectiveSize <- function(object, Param, na.rm = TRUE) {
  HiddenGetMeasure(object, Param, Measure = "effectiveSize", na.rm)
}
#' @rdname BASiCS_EffectiveSize
#' @export
BASiCS_effectiveSize <- BASiCS_EffectiveSize
