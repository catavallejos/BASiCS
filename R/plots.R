
#' @rdname BASiCS_showFit
#'
#' @title Plotting the trend after Bayesian regression
#'
#' @description Plotting the trend after Bayesian regression using a
#' \code{\linkS4class{BASiCS_Chain}} object
#'
#' @param object an object of class \code{\linkS4class{BASiCS_Chain}}
#' @param xlab As in \code{\link[graphics]{par}}.
#' @param ylab As in \code{\link[graphics]{par}}.
#' @param pch As in \code{\link[graphics]{par}}. Default value \code{pch = 16}.
#' @param smooth Logical to indicate wether the smoothScatter function is used
#' to plot the scatter plot. Default value \code{smooth = TRUE}. 
#' @param variance Variance used to build GRBFs for regression. Default value
#' \code{variance = 1.2} 
#' @param colour colour used to denote genes within the scatterplot. Only used 
#' when \code{smooth = TRUE}. Default value 
#' \code{colour = "dark blue"}. 
#' @param markExcludedGenes Whether or not lowly expressed genes that were
#' excluded from the regression fit are included in the scatterplot.
#' Default value \code{markExcludedGenes = TRUE}.
#' @param GenesSel Vector of gene names to be highlighted in the scatterplot.
#' Only used when \code{smooth = TRUE}. Default value \code{GenesSel = NULL}.
#' @param colourGenesSel colour used to denote the genes listed in 
#' \code{GenesSel} within the scatterplot. Default value 
#' \code{colourGenesSel = "dark red"}.
#' @param ... Additional parameters for plotting.
#'
#' @examples
#' data(ChainRNAReg)
#' BASiCS_showFit(ChainRNAReg)
#'
#' @return A ggplot2 object
#'
#' @author Nils Eling \email{eling@@ebi.ac.uk}
#' @author Catalina Vallejos \email{cnvallej@@uc.cl}
#'
#' @references
#' Eling et al (2018). Cell Systems
#' https://doi.org/10.1016/j.cels.2018.06.011
#' @export
BASiCS_showFit <- function(object,
                           xlab = "log(mu)",
                           ylab = "log(delta)",
                           pch = 16,
                           smooth = TRUE,
                           variance = 1.2,
                           colour = "dark blue",
                           markExcludedGenes = TRUE,
                           GenesSel = NULL,
                           colourGenesSel = "dark red",
                           ...) {

  if (!inherits(object, "BASiCS_Chain")) {
    stop(paste0("Incorrect class for object:", class(object)))
  }

  if (!("beta" %in% names(object@parameters))) {
    stop("'beta' is missing. Regression was not performed.")
  }


  m <- log(object@parameters$mu[1,])
  grid.mu <- seq(round(min(m), digits = 2),
                 round(max(m), digits = 2), length.out = 1000)

  # Create design matrix across the grid
  n <- ncol(object@parameters$beta)
  range <- diff(range(grid.mu))
  myu <- seq(min(grid.mu), by = range/(n-3), length.out = n-2)
  h <- diff(myu)*variance

  B <- matrix(1,length(grid.mu),n)
  B[,2] <- grid.mu
  for (j in seq_len(n-2)) {
    B[,j+2] = exp(-0.5 * (grid.mu - myu[j])^2 / (h[1]^2))
  }

  # Calculate yhat = X*beta
  yhat <- apply(object@parameters$beta, 1, function(n) B%*%n)
  yhat.HPD <- coda::HPDinterval(coda::mcmc(t(yhat)), 0.95)

  df <- data.frame(mu = log(colMedians(object@parameters$mu)),
                   delta = log(colMedians(object@parameters$delta)),
                   included = !is.na(object@parameters$epsilon[1,]))
  rownames(df) <- colnames(object@parameters$mu)

  df2 <- data.frame(mu2 = grid.mu,
                    yhat = rowMedians(yhat),
                    yhat.upper = yhat.HPD[,2],
                    yhat.lower = yhat.HPD[,1])

  plot.out <- ggplot(df[df$included,],
                     ggplot2::aes_string(x = "mu", y = "delta")) +
    xlab(xlab) + ylab(ylab) +
    theme_minimal(base_size = 15)
  if(markExcludedGenes == TRUE) {
    plot.out <- plot.out + 
    ggplot2::geom_point(data = df[!df$included, ],
                        shape = pch,
                        colour = "purple",
                        alpha = 0.3)
  }

  if (smooth) {
    cols <- c("dark blue", "yellow", "dark red")
    plot.out <- plot.out +
      ggplot2::geom_hex(bins = 100) +
      ggplot2::scale_fill_gradientn(name = "",
                                    colours = grDevices::colorRampPalette(cols)(100),
                                    guide = FALSE)
  }
  else {
    plot.out <- plot.out +
      ggplot2::geom_point(shape = pch,
                          colour = colour)
  }
  if(!is.null(GenesSel)) {
    if(sum(GenesSel %in% rownames(df)) != length(GenesSel)) {
      stop("Some elements of `GenesSel` are not found in the data.")
    }
    plot.out <- plot.out +
      ggplot2::geom_point(data = df[rownames(df) %in% GenesSel,],
                          shape = pch,
                          colour = colourGenesSel)
  }
  
  plot.out <- plot.out + 
    ggplot2::geom_line(data = df2,
                       inherit.aes = FALSE,
                       mapping = ggplot2::aes_string(x = "mu2", y = "yhat"),
                       colour = "dark red") +
    ggplot2::geom_ribbon(data = df2,
                         inherit.aes = FALSE,
                         mapping = ggplot2::aes_string(x = "mu2",
                                             ymin = "yhat.lower",
                                             ymax = "yhat.upper"),
                         alpha = 0.5)
    
  return(plot.out)
}


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
#' 
#' @return A ggplot object.
#'
#' @examples
#'
#' # See
#' help(BASiCS_MCMC)
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
                            LogY = isTRUE(y %in% c("mu", "delta"))) {
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
  metric <- HiddenGetMeasure(object, Param, Measure)
  sX <- if (LogX) scale_x_log10() else scale_x_continuous()
  sY <- if (LogY) scale_y_log10() else scale_y_continuous()

  if (!is.null(x)) {
    xMat <- HiddenGetParam(object, x)
    yMat <- HiddenGetParam(object, y)
    df <- data.frame(
      x = colMedians(xMat),
      y = colMedians(yMat),
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
      ggplot2::geom_hex(ggplot2::aes_string(fill = "..density..")) +
      viridis::scale_fill_viridis(name = "Density") +
      sX + sY +
      ggplot2::labs(x = Param, y = HiddenScaleName(Measure))
  }
}

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
#' \code{'sigma2'} and \code{'epsilon'}.
#'
#' @return A ggplot object.
#'
#' @examples
#'
#' # See
#' help(BASiCS_MCMC)
#'
#' @seealso \code{\linkS4class{BASiCS_Chain}}
#'
#' @author Alan O'Callaghan \email{a.b.ocallaghan@sms.ed.ac.uk}
#'
#' @export
BASiCS_diagHist <- function(object, Param = NULL) {
  if (!inherits(object, "BASiCS_Chain")) {
    stop(paste0("Incorrect class for object:", class(object)))
  }
  Measure <- "effectiveSize"

  if (is.null(Param)) {
    metric <- lapply(names(object@parameters), function(param) {
      HiddenGetMeasure(object, param, Measure)
    })
    metric <- Reduce(c, metric)              
  } else {
    metric <- HiddenGetMeasure(object, Param, Measure)
  }
  if (length(metric) == 1) {
    stop(paste0("Cannot produce histogram of a single value (", metric, ")"))
  }
  ggplot2::ggplot(mapping = ggplot2::aes(x = metric)) + 
    ggplot2::geom_histogram(bins = grDevices::nclass.FD(metric)) +
    ggplot2::labs(x = HiddenScaleName(Measure, Param),
                  y = "Count")
}

HiddenScaleName <- function(Measure = c("effectiveSize",
                                  "geweke.diag"),
                      Param = NULL) {
  Measure <- match.arg(Measure)
  measure_name <- switch(Measure, 
    effectiveSize = "Effective sample size",
    geweke.diag = "Geweke diagnostic"
  )
  if (!is.null(Param)) {
    measure_name <- paste0(measure_name, ": ", Param)
  }
  measure_name
}
HiddenGetMeasure <- function(object, 
                       Param,
                       Measure = c("effectiveSize",
                                   "geweke.diag")) {
  Measure <- match.arg(Measure)
  MeasureFun <- match.fun(Measure)
  metric <- MeasureFun(coda::mcmc(HiddenGetParam(object, Param)))
  if (Measure == "geweke.diag") {
    metric <- metric$z
  }
  metric
}
