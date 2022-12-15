#' @title Create diagnostic plots of MCMC parameters
#'
#' @description Plot parameter values and effective sample size.
#' See \link[coda]{effectiveSize}
#' for more details on this diagnostic measure.
#'
#' @param object an object of class \code{\linkS4class{BASiCS_Summary}}
#' @param Parameter Name of the parameter to be plotted.
#' Default \code{Parameter = 'mu'}
#' @param Measure Character scalar specifying the diagnostic measure to plot.
#' Current options are effective sample size, the Geweke diagnostic
#' criterion, and the \code{\link[posterior]{rhat}} diagnostic.
#' @param x,y Optional MCMC parameter values to be plotted on the x or y axis,
#' respectively. If neither is supplied, Parameter will be plotted on the x axis
#' and effective sample size will be plotted on the y axis as
#' a density plot.
#' @param LogX,LogY A logical value indicating whether to use a log10
#' transformation for the x or y axis, respectively.
#' @param Smooth A logical value indicating whether to use smoothing
#' (specifically hexagonal binning using \code{\link[ggplot2]{geom_hex}}).
#' @param HLine Numeric scalar or vector indicating threshold value(s) to be
#' displayed as a dashed line on the plot when \code{DrawHLine = TRUE}.
#' Alternatively, can be set to \code{FALSE} to disable line drawing,
#' or \code{TRUE} to use the default thresholds.
#' @param na.rm Logical value indicating whether NA values should be removed
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
#' # Point estimates versus effective sample size
#' BASiCS_DiagPlot(ChainSC, Parameter = "mu")
#' # Effective sample size as colour, mu as x, delta as y.
#' BASiCS_DiagPlot(ChainSC, x = "mu", y = "delta")
#' 
#' # Point estimates versus Geweke diagnostic
#' BASiCS_DiagPlot(ChainSC, Parameter = "mu", Measure = "geweke.diag")
#'
#' @seealso \code{\linkS4class{BASiCS_Chain}}
#'
#' @author Alan O'Callaghan
#'
#' @export
BASiCS_DiagPlot <- function(object, 
                            Parameter = "mu",
                            Measure = c("ess", "geweke.diag", "rhat"),
                            x = NULL,
                            y = NULL,
                            LogX = isTRUE(x %in% c("mu", "delta")),
                            LogY = isTRUE(y %in% c("mu", "delta")),
                            Smooth = TRUE,
                            HLine = TRUE,
                            na.rm = TRUE) {

  Measure <- match.arg(Measure)
  if (!inherits(object, "BASiCS_Chain")) {
    stop(paste0("Incorrect class for object:", class(object)))
  }
  if (!is.null(x) || !is.null(y)) {
    if (!is.null(x) && is.null(y)) {
      stop("Must specify both x and y or neither!")
    }
  } else {
    LogX <- Parameter %in% c("mu", "delta")
  }
  .CheckValidCombination(x, y, Parameter)
  metric <- .GetMeasure(object, Parameter, Measure, na.rm)
  sX <- if (LogX) scale_x_log10() else scale_x_continuous()
  sY <- if (LogY) scale_y_log10() else scale_y_continuous()

  DrawHLine <- TRUE
  ## if it's logical, we only care if it's FALSE
  if (is.logical(HLine)) {
    DrawHLine <- HLine
  }
  ## not numeric, assume it's NULL or logical, set to default
  if (!is.numeric(HLine)) {
    HLine <- .LineAt(Measure)
  }

  ## x is one param, y is another, colour is measure
  if (!is.null(x)) {
    xMat <- .GetParam(object, x)
    yMat <- .GetParam(object, y)
    xMat <- xMat[, names(metric), drop = FALSE]
    yMat <- yMat[, names(metric), drop = FALSE]

    df <- data.frame(
      x = matrixStats::colMedians(xMat),
      y = matrixStats::colMedians(yMat),
      metric = metric
    )
    df <- df[order(df$metric), ]
    g <- ggplot(df, aes(x = .data$x, y = .data$y, colour = .data$metric)) + 
      geom_point(alpha = 0.5, shape = 16) +
      viridis::scale_colour_viridis(name = .ScaleName(Measure, Parameter)
        #, trans="log10"
        ) +
      sX + sY +
      labs(x = x, y = y)
  ## x is param, y is measure
  } else {
    xMat <- .GetParam(object, Parameter)
    xMat <- xMat[, names(metric), drop = FALSE]

    g <- ggplot() +
      aes(x = matrixStats::colMedians(xMat), y = metric) +
      sX + sY +
      labs(x = Parameter, y = .ScaleName(Measure))
    if (Smooth) {
      g <- g +
        geom_hex() +
        viridis::scale_fill_viridis(name = "Count", trans = "log10")
    } else {
      g <- g +
        geom_point()
    }
    if (DrawHLine) {
      g <- g +
        geom_hline(yintercept = HLine, linetype = "dashed")
    }
    g
  }
}
#' @export
#' @rdname BASiCS_DiagPlot
BASiCS_diagPlot <- function(...) {
  .Deprecated("BASiCS_DiagPlot")
  BASiCS_DiagPlot(...)
}
