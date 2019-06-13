#' Produce plots assessing differential expression results
#' 
#' @param ResultsDE A BASiCS_ResultsDE object
#' @param Which 
#' @param ... Passed to methods
setMethod("BASiCS_DEPlots", signature(object = "BASiCS_ResultsDE"),
  function(object, Which = c("MAPlot", "VolcanoPlot", "GridPlot"), ...) {
    l <- lapply(ResultsDE@Results, BASiCS_DEPlots, Which = Which, ...)
    if (length(Which) > 1) {
      nrow <- length(ResultsDE@Results)
      labels <- lapply(ResultsDE@Results, function(x) {
        cap(MeasureName(x@Name))
      })
      labels <- Reduce(c, labels)
      labels <- c(labels, rep("", length(ResultsDE@Results) * length(Which)))
    } else {
      nrow <- length(Which)
      labels <- sapply(ResultsDE@Results, function(x) cap(MeasureName(x@Name)))
    }
    cowplot::plot_grid(plotlist = l, nrow = nrow, labels = labels)
  }
)

setMethod("BASiCS_DEPlots", signature(object = "BASiCS_ResultDE"),
  function(
    object,
    GroupLabel1 = ResultDE@GroupLabel1,
    GroupLabel2 = ResultDE@GroupLabel2,
    ProbThresholds = seq(0.5, 0.9995, by = 0.00025),
    Epsilon = ResultDE@Epsilon,
    EFDR = ResultDE@EFDR,
    Table = ResultDE@Table,
    Measure = ResultDE@Name,
    EFDRgrid = ResultDE@EFDRgrid,
    EFNRgrid = ResultDE@EFNRgrid,
    ProbThreshold = ResultDE@ProbThreshold,
    Which = c("MAPlot", "VolcanoPlot", "GridPlot")
  ){

    BASiCS_DEPlots(
      GroupLabel1 = GroupLabel1,
      GroupLabel2 = GroupLabel2,
      ProbThresholds = ProbThresholds,
      Epsilon = Epsilon,
      EFDR = EFDR,
      Table = Table,
      Measure = Measure,
      EFDRgrid = EFDRgrid,
      EFNRgrid = EFNRgrid,
      ProbThreshold = ProbThreshold,
      Which = Which
    )
  }
)
setMethod("BASiCS_DEPlots", signature(object = "missing"),
  function(
    GroupLabel1,
    GroupLabel2,
    ProbThresholds = seq(0.5, 0.9995, by = 0.00025),
    Epsilon,
    EFDR,
    Table,
    Measure,
    EFDRgrid,
    EFNRgrid,
    ProbThreshold,
    Which = c("MAPlot", "VolcanoPlot", "GridPlot")
  ) {

    Which <- match.arg(Which, several.ok = TRUE)

    Plots <- list()

    if ("MAPlot" %in% Which) {
      Plots <- c(Plots, 
        list(MAPlot(Measure, Table, GroupLabel1, GroupLabel2, Epsilon))
      )
    }
    if ("VolcanoPlot" %in% Which) {
      Plots <- c(Plots, 
        list(VolcanoPlot(Measure, Table, GroupLabel1, GroupLabel2, Epsilon))
      )
    }
    if ("GridPlot" %in% Which) {
      Plots <- c(
        list(GridPlot(Measure, EFDR, ProbThresholds, EFDRgrid, EFNRgrid, ProbThreshold)),
        Plots
      )
    }
    if (length(Plots) == 1) {
      Plots <- Plots[[1]]
    } else {
      Plots <- cowplot::plot_grid(plotlist = Plots, nrow = 1)
    }
    Plots
  }
)



GridPlot <- function(Measure, 
                       EFDR,
                       ProbThresholds = seq(0.5, 0.9995, by = 0.00025),
                       EFDRgrid,
                       EFNRgrid,
                       ProbThreshold
                       ) {
  df <- data.frame(
    ProbThresholds, 
    EFDR = EFDRgrid, 
    EFNR = EFNRgrid)

  mdf <- reshape2::melt(df,measure.vars = c("EFDR", "EFNR"))
  ggplot2::ggplot(mdf, aes(x = ProbThresholds, y = value, color = variable)) +
    ggplot2::geom_line() +
    ggplot2::labs(
      y = "Error rate", 
      x = "Probability threshold"
      # ,title = paste("Differential", MeasureName(Measure))
    ) +
    ggplot2::ylim(c(0, 1)) +
    ggplot2::scale_color_discrete(name = ""
      # , palette = "Dark2"
      ) +
    ggplot2::geom_hline(yintercept = EFDR, col = "steelblue", lty = 2) +
    ggplot2::geom_vline(xintercept = ProbThreshold, col = "indianred", lty = 1)
}



VolcanoPlot <- function(Measure, Table, GroupLabel1, GroupLabel2, Epsilon) {
  IndDiff <- !Table[[paste0("ResultDiff", Measure)]] %in% c("ExcludedByUser", "NoDiff")
  # bins <- NClassFD2D(
  #   Table[[paste0(Measure, DistanceVar(Measure))]],
  #   Table[[paste0("ProbDiff", Measure)]]
  # )
  bins <- 100
  ggplot2::ggplot(
      Table, 
      ggplot2::aes_string(
        x = paste0(Measure, DistanceVar(Measure)), 
        y = paste0("ProbDiff", Measure))
    ) +
    # ggplot2::geom_point() +
    ggplot2::geom_hex(bins = bins, aes(fill = ..density..), na.rm = TRUE) +
    ggplot2::geom_point(
      data = Table[IndDiff, ], 
      shape = 16, 
      col = "violetred", 
      na.rm = TRUE) +
    ggplot2::geom_hline(
      yintercept = c(-Epsilon, Epsilon), 
      lty = "dashed", 
      color = "grey40", 
      na.rm = TRUE) +
    ggplot2::ylim(c(0, 1)) +
    viridis::scale_fill_viridis(name = "Density") +
    ggplot2::labs(
      x = paste(cap(LogDistanceName(Measure)), GroupLabel1, "vs", GroupLabel2),
      y = "Posterior probability"
      # ,title = paste("Differential", MeasureName(Measure), "test")
    )
}

