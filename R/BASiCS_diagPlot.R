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
#' \code{'sigma2'} and \code{'epsilon'}.
#' @param x,y Optional MCMC parameter values to be plotted on the x or y axis, 
#' respectively. If neither is supplied, Param will be plotted on the x axis
#' and \code{coda::effectiveSize(Param)} will be plotted on the y axis as
#' a density plot.
#' @param LogX,LogY A boolean value indicating whether to use a log10
#' transformation for the x or y axis, respectively.
#' @param na.rm Logical value indicating whether NA values should be removed
#' before calculating effective sample size.
#' 
#' @return A ggplot object.
#'
#' @examples
#'
#' # Built-in simulated dataset
#' Data <- makeExampleBASiCS_Data()
#' # To analyse real data, please refer to the instructions in:
#' # https://github.com/catavallejos/BASiCS/wiki/2.-Input-preparation
#'
#' # Only a short run of the MCMC algorithm for illustration purposes
#' # Longer runs migth be required to reach convergence
#' Chain <- BASiCS_MCMC(Data, N = 10, Thin = 2, Burn = 4, Regression = FALSE,
#'                      PrintProgress = FALSE, WithSpikes = TRUE)
#' # Plot effective sample size for mu (y) against mu, coloured by density
#' BASiCS_diagPlot(Chain)
#' # Effective sample size as colour, mu as x, delta as y.
#' BASiCS_diagPlot(Chain, x = "mu", y = "delta")
#' 
#' @seealso \code{\linkS4class{BASiCS_Chain}}
#'
#' @author Alan O'Callaghan \email{a.b.ocallaghan@sms.ed.ac.uk}
#' 
#' @export
BASiCS_diagPlot <- function(object, 
                            Param = "mu", 
                            x = NULL, 
                            y = NULL,
                            LogX = isTRUE(x %in% c("mu", "delta")),
                            LogY = isTRUE(y %in% c("mu", "delta")),
                            na.rm = FALSE) {


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
    df <- data.frame(
      x = matrixStats::colMedians(xMat),
      y = matrixStats::colMedians(yMat),
      metric = metric
    )

    ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, color = metric)) + 
      ggplot2::geom_point(alpha = 0.5, shape = 16) +
      viridis::scale_color_viridis(name = HiddenScaleName(Measure, Param)) +
      sX + sY +
      ggplot2::labs(x = x, y = y)              
  } else {
    xMat <- HiddenGetParam(object, Param)
    df <- data.frame(
      x = matrixStats::colMedians(xMat),
      y = metric
    )
    ggplot2::ggplot(df, ggplot2::aes(x = x, y = metric)) + 
      ggplot2::geom_hex( #ggplot2::aes_string(fill = "..density..")
        ) +
      viridis::scale_fill_viridis(name = "Density", trans = "log10") +
      sX + sY +
      ggplot2::labs(x = Param, y = HiddenScaleName(Measure))
  }
}
