#' @title Create diagnostic plots of MCMC parameters
#'
#' @description Plot a histogram of effective sample size or Geweke's diagnostic
#' z-statistic. See \link[coda]{effectiveSize} and \link[coda]{geweke.diag} for
#' more details.
#'
#' @param object an object of class \code{\linkS4class{BASiCS_Summary}}
#' @param Parameter Optional name of a chain parameter to restrict the
#' histogram; if not supplied, all parameters will be assessed.
#' Default \code{Parameter = NULL}.
#' @param Measure Character scalar specifying the diagnostic measure to plot.
#' Current options are effective sample size, the Geweke diagnostic
#' criterion, and the \code{\link[posterior]{rhat}} diagnostic.
#' @param VLine Numeric scalar indicating a threshold value to be displayed as
#' a dashed line on the plot.
#' Alternatively, can be set to \code{FALSE} to disable line drawing,
#' or \code{TRUE} to use the default thresholds.
#' @param na.rm Logical scalar indicating whether NA values should be removed
#' before calculating effective sample size.
#' @param ... Unused.
#'
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
#' BASiCS_DiagHist(ChainSC, Parameter = "mu")
#'
#' @seealso \code{\linkS4class{BASiCS_Chain}}
#' @references
#'  Geweke, J. Evaluating the accuracy of sampling-based approaches to
#'  calculating posterior moments. In _Bayesian Statistics 4_ (ed JM
#'  Bernado, JO Berger, AP Dawid and AFM Smith). Clarendon Press,
#'  Oxford, UK.
#'
#' @author Alan O'Callaghan
#'
#' @export
BASiCS_DiagHist <- function(
    object,
    Parameter = NULL,
    Measure = c("ess", "geweke.diag", "rhat"),
    VLine = TRUE,
    na.rm = TRUE) {

  if (!inherits(object, "BASiCS_Chain")) {
    stop(paste0("Incorrect class for object:", class(object)))
  }
  Measure <- match.arg(Measure)
  ## if it's logical, we only care if it's FALSE
  DrawVLine <- TRUE
  if (is.logical(VLine)) {
    DrawVLine <- VLine
  }
  ## not numeric, assume it's NULL or logical, set to default
  if (!is.numeric(VLine)) {
    VLine <- .LineAt(Measure)
    DrawVLine <- TRUE
  }

  if (is.null(Parameter)) {
    metric <- lapply(names(object@parameters), function(param) {
      try(.GetMeasure(object, param, Measure, na.rm), silent = TRUE)
    })
    ind_error <- vapply(
      metric,
      function(x) inherits(x, "try-error"),
      logical(1)
    )
    metric <- metric[!ind_error]
    if (all(ind_error)) {
      stop("ess failed for all parameters.")
    } else if (any(ind_error)) {
      warning(
        paste("ess failed for some parameters:",
          paste(names(object@parameters)[ind_error], collapse = ", ")
        )
      )
    }
    metric <- Reduce(c, metric)
  } else {
    metric <- .GetMeasure(object, Parameter, Measure, na.rm)
  }
  if (length(metric) == 1) {
    stop(paste0("Cannot produce histogram of a single value (", metric, ")"))
  }
  ggplot(mapping = aes(x = metric)) +
    geom_histogram(bins = grDevices::nclass.FD(metric)) +
    (if (DrawVLine) geom_vline(xintercept = VLine, linetype = "dashed") else list()) +
    labs(x = .ScaleName(Measure, Parameter),
                  y = "Count")
}

#' @export
#' @rdname BASiCS_DiagHist
BASiCS_diagHist <- function(...) {
  .Deprecated("BASiCS_DiagHist")
  BASiCS_DiagHist(...)
}

.LineAt <- function(Measure) {
  if (Measure == "ess") {
    100
  } else if (Measure == "geweke.diag") {
    c(-3, 3)
  } else if (Measure == "rhat") {
    1.05
  }
}
