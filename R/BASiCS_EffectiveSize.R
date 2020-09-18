#' Calculate effective sample size for BASiCS_Chain parameters
#' 
#' @description A function to calculate effective sample size
#' \code{\linkS4class{BASiCS_Chain}} objects.
#' 
#' @param object an object of class \code{\linkS4class{BASiCS_Chain}}.
#' @param Parameter The parameter to use to calculate effective sample size.
#' Possible
#' values: \code{'mu'}, \code{'delta'}, \code{'phi'}, \code{'s'}, 
#' \code{'nu'}, \code{'theta'}, \code{'beta'}, \code{'sigma2'} and 
#' \code{'epsilon'}. 
#' @param na.rm Remove \code{NA} values before calculating effective sample 
#' size. Only relevant when \code{Parameter = "epsilon"} (genes with very low 
#' expression are excluding when infering the mean/over-dispersion trend. 
#' Default: \code{na.rm = TRUE}.
#' @param ... Unused.
#' 
#' @return A vector with effective sample sizes for all the elements 
#' of \code{Parameter}
#' 
#' @examples
#' 
#' data(ChainSC)
#' BASiCS_EffectiveSize(ChainSC, Parameter = "mu")
#' 
#' @export
BASiCS_EffectiveSize <- function(object, Parameter, na.rm = TRUE) {
  .GetMeasure(object, Parameter, Measure = "ess", na.rm)
}

#' @rdname BASiCS_EffectiveSize
#' @export
BASiCS_effectiveSize <- function(...) {
  .Deprecated("BASiCS_EffectiveSize")
  BASiCS_EffectiveSize(...)
}
