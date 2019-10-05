#' @title Create diagnostic plots of MCMC parameters
#'
#' @description Plot parameter values and effective sample size.
#' See \link[coda]{effectiveSize}
#' for more details on this diagnostic measure.
#'
#' @param object an object of class \code{\linkS4class{BASiCS_Summary}}
#' @param Param Optional name of a chain parameter to restrict the histogram;
#' if not supplied, all parameters will be assessed.
#' Possible values: \code{'mu'}, \code{'delta'}, \code{'phi'},
#' \code{'s'}, \code{'nu'}, \code{'theta'}, \code{'beta'},
#' \code{'sigma2'} and \code{'epsilon'}. Default \code{Param = 'mu'}
#' @param x,y Optional MCMC parameter values to be plotted on the x or y axis, 
#' respectively. If neither is supplied, Param will be plotted on the x axis
#' and \code{coda::effectiveSize(Param)} will be plotted on the y axis as
#' a density plot.
#' @param LogX,LogY A logical value indicating whether to use a log10
#' transformation for the x or y axis, respectively.
#' @param Smooth A logical value indicating whether to use smoothing 
#' (specifically hexagonal binning using \code{\link[ggplot2]{geom_hex}}).
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
#' BASiCS_DiagPlot(ChainSC, Param = "mu")
#' # Effective sample size as colour, mu as x, delta as y.
#' BASiCS_DiagPlot(ChainSC, x = "mu", y = "delta")
#' 
#' @seealso \code{\linkS4class{BASiCS_Chain}}
#'
#' @author Alan O'Callaghan \email{a.b.ocallaghan@sms.ed.ac.uk}
#' 
#' @export
BASiCS_DiagPlot <- function(object, 
                            Param = "mu", 
                            x = NULL, 
                            y = NULL,
                            LogX = isTRUE(x %in% c("mu", "delta")),
                            LogY = isTRUE(y %in% c("mu", "delta")),
                            Smooth = TRUE,
                            na.rm = TRUE) {


  if (!inherits(object, "BASiCS_Chain")) {
    stop(paste0("Incorrect class for object:", class(object)))
  }
  if (!is.null(x) || !is.null(y)) {
    if (!is.null(x) && is.null(y)) {
      stop("Must specify both x and y or neither!")
    }
  } else {
    LogX <- Param %in% c("mu", "delta")
  }
  Measure <- "effectiveSize"
  HiddenCheckValidCombination(x, y, Param)
  metric <- HiddenGetMeasure(object, Param, Measure, na.rm)
  sX <- if (LogX) ggplot2::scale_x_log10() else ggplot2::scale_x_continuous()
  sY <- if (LogY) ggplot2::scale_y_log10() else ggplot2::scale_y_continuous()

  if (!is.null(x)) {
    xMat <- HiddenGetParam(object, x)
    yMat <- HiddenGetParam(object, y)
    xMat <- xMat[, names(metric), drop = FALSE]
    yMat <- yMat[, names(metric), drop = FALSE]

    df <- data.frame(
      x = matrixStats::colMedians(xMat),
      y = matrixStats::colMedians(yMat),
      metric = metric
    )
    df <- df[order(df$metric), ]
    g <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, color = metric)) + 
      ggplot2::geom_point(alpha = 0.5, shape = 16) +
      viridis::scale_color_viridis(name = HiddenScaleName(Measure, Param)
        #, trans="log10"
        ) +
      sX + sY +
      ggplot2::labs(x = x, y = y)              
  } else {
    xMat <- HiddenGetParam(object, Param)
    xMat <- xMat[, names(metric), drop = FALSE]
    df <- data.frame(
      x = matrixStats::colMedians(xMat),
      y = metric
    )
    g <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = metric)) + 
      sX + sY +
      ggplot2::labs(x = Param, y = HiddenScaleName(Measure))
    if (Smooth) {
      g <- g +
        ggplot2::geom_hex() +
        viridis::scale_fill_viridis(name = "Count", trans = "log10")
    } else {
      g <- g +
        ggplot2::geom_point()
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
