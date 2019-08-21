#' @title Create diagnostic plots of MCMC parameters
#'
#' @description Plot a histogram of effective sample size or Geweke's diagnostic
#' z-statistic. See \link[coda]{effectiveSize} and \link[coda]{geweke.diag} for
#' more details.
#'
#' @param object an object of class \code{\linkS4class{BASiCS_Summary}}
#' @param Param Optional name of a chain parameter to restrict the histogram;
#' if not supplied, all parameters will be assessed.
#' Possible values: \code{'mu'}, \code{'delta'}, \code{'phi'},
#' \code{'s'}, \code{'nu'}, \code{'theta'}, \code{'beta'},
#' \code{'sigma2'} and \code{'epsilon'}. Default \code{Param = NULL}.
#' @param na.rm Logical value indicating whether NA values should be removed
#' before calculating effective sample size.
#' @param ... Unused
#' @return A ggplot object.
#'
#' @examples
#'
#' # Built-in example chain
#' data(ChainSC)
#' 
#' # See effective sample size distribution across all parameters
#' BASiCS_DiagHist(ChainSC)
#' # For mu only
#' BASiCS_DiagHist(ChainSC, Param = "mu")
#' 
#' @seealso \code{\linkS4class{BASiCS_Chain}}
#'
#' @author Alan O'Callaghan \email{a.b.ocallaghan@sms.ed.ac.uk}
#'
#' @export
BASiCS_DiagHist <- function(object, Param = NULL, na.rm = TRUE) {
  if (!inherits(object, "BASiCS_Chain")) {
    stop(paste0("Incorrect class for object:", class(object)))
  }
  Measure <- "effectiveSize"

  if (is.null(Param)) {
    metric <- lapply(names(object@parameters), function(param) {
      try(HiddenGetMeasure(object, param, Measure, na.rm), silent = TRUE)
    })
    ind_error <- vapply(
      metric,
      function(x) inherits(x, "try-error"),
      logical(1)
    )
    metric <- metric[!ind_error]
    if (all(ind_error)) {
      stop("coda::effectiveSize failed for all parameters.")
    } else if (any(ind_error)) {
      warning(
        paste("coda::effectiveSize failed for some parameters:",
          paste(names(object@parameters)[ind_error], collapse = ", ")
        )
      )
    }
    metric <- Reduce(c, metric)              
  } else {
    metric <- HiddenGetMeasure(object, Param, Measure, na.rm)
  }
  if (length(metric) == 1) {
    stop(paste0("Cannot produce histogram of a single value (", metric, ")"))
  }
  ggplot2::ggplot(mapping = ggplot2::aes(x = metric)) + 
    ggplot2::geom_histogram(bins = grDevices::nclass.FD(metric)) +
    ggplot2::labs(x = HiddenScaleName(Measure, Param),
                  y = "Count")
}

#' @export
#' @rdname BASiCS_DiagHist
BASiCS_diagHist <- function(...) {
  .Deprecated("BASiCS_DiagHist")
  BASiCS_DiagHist(...)
}

