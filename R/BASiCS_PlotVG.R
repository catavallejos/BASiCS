#' Plots of HVG/LVG search.
#' @param object \linkS4class{BASiCS_ResultVG} object.
#' @param Plot Character scalar specifying the type of plot to be made.
#' Options are "Grid" and "VG".
#' @param ... Optional graphical parameters passed to \code{.VGPlot} 
#' (internal function). 
#' 
#' @return A plot.
#' @examples
#' data(ChainSC)
#'
#' # Highly and lowly variable genes detection (within a single group of cells)
#' DetectHVG <- BASiCS_DetectHVG(ChainSC, VarThreshold = 0.60,
#'                               EFDR = 0.10, Plot = TRUE)
#' BASiCS_PlotVG(DetectHVG)
#' @export
BASiCS_PlotVG <- function(object, Plot = c("Grid", "VG"), ...) {
  Plot <- match.arg(Plot)
  if (Plot == "Grid") {
    .GridPlot(
      EFDR = object@EFDR,
      ProbThresholds = object@ProbThresholds,
      EFDRgrid = object@EFDRgrid,
      EFNRgrid = object@EFNRgrid,
      ProbThreshold = object@ProbThreshold
    )
  } else if (Plot == "VG") {
    .VGPlot(
      Task = object@Name,
      Mu = object@Table$Mu,
      Prob = object@Table$Prob,
      OptThreshold = object@ProbThreshold,
      Hits = object@Table[[object@Name]],
      ...
    )
  }
}

.VGPlot <- function(
  Task,
  Mu,
  Prob,
  OptThreshold,
  Hits,
  ylim = c(0, 1),
  xlim = c(min(Mu), max(Mu)),
  cex = 1.5,
  pch = 16,
  col = 8,
  bty = "n",
  xlab = "Mean expression",
  ylab = paste(Task, "probability"),
  title = ""
) {
  
  df <- data.frame(Mu, Prob)
  ggplot(df, aes(x = .data$Mu, y = .data$Prob)) +
    geom_point(
      aes(colour = ifelse(Hits, Task, "Other")),
      pch = pch, cex = cex) +
    scale_colour_brewer(palette = "Set1", name = "") +
    geom_hline(
      yintercept = OptThreshold[[1]], lty = 2, col = "black"
    ) +
    scale_x_log10() +
    labs(
      x = xlab,
      y = ylab,
      title = title
    ) +
    theme_classic()
}
