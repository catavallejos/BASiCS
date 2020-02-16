#' Calculate effective sample size for BASiCS_Chain parameters
#' 
#' @description A function to calculate effective sample size
#' \code{\linkS4class{BASiCS_Chain}} objects.
#' 
#' @param object an object of class \code{\linkS4class{BASiCS_Chain}}.
#' @param Param The parameter to use to calculate effective sample size.
#' Possible
#' values: \code{'mu'}, \code{'delta'}, \code{'phi'}, \code{'s'}, 
#' \code{'nu'}, \code{'theta'}, \code{'beta'}, \code{'sigma2'} and 
#' \code{'epsilon'}. 
#' @param na.rm Remove \code{NA} values before calculating effective sample 
#' size. Only relevant when \code{Param = "epsilon"} (genes with very low 
#' expression are excluding when infering the mean/over-dispersion trend. 
#' Default: \code{na.rm = TRUE}.
#' @param ... Unused.
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
  HiddenGetMeasure(object, Param, Measure = "ess", na.rm)
}

#' @rdname BASiCS_EffectiveSize
#' @export
BASiCS_effectiveSize <- function(...) {
  .Deprecated("BASiCS_EffectiveSize")
  BASiCS_EffectiveSize(...)
}
