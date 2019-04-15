#' @export
BASiCS_PlotOffset <- function(OffsetCorrected,
                              Type = c(
                                "offset estimate", 
                                "before-after",
                                "MA plot"
                              ),
                              Print = FALSE
                              ) {
  Type <- match.arg(Type)

  plots <- list()
  if ("offset estimate" %in% Type) {
    # Offset uncertainty
    g <- ggplot2::ggplot(mapping = ggplot2::aes(y = OffsetCorrected@OffsetChain)) +
      ggplot2::geom_boxplot() +
      ggplot2::labs(title = "Offset MCMC chain", y = "Offset estimate", x = NULL)
    if (Print) {
      print(g)
    }
    plots[[Type]] <- g
  }

  # Mean expression parameters before/after offset correction
  if ("before-after" %in% Type) {
    df <- cbind(OffsetCorrected@Mu1_old, OffsetCorrected@Mu2)
    colnames(df) <- c(OffsetCorrected@GroupLabel1, OffsetCorrected@GroupLabel2)
    mdf <- reshape2::melt(df)
    ggplot2::ggplot(mdf, ggplot2::aes(x = Var2, y = value)) +
      ggplot2::labs(title = "Before correction", y = "Mean expression", x = NULL) +
      ggplot2::geom_boxplot() +
      ggplot2::scale_y_log10()

    df <- cbind(OffsetCorrected@Mu1, OffsetCorrected@Mu2)
    colnames(df) <- c(OffsetCorrected@GroupLabel1, OffsetCorrected@GroupLabel2)
    mdf <- reshape2::melt(df)
    g <- ggplot2::ggplot(mdf, ggplot2::aes(x = Var2, y = value)) +
      ggplot2::labs(title = "After correction", y = "Mean expression", x = NULL) +
      ggplot2::geom_boxplot() +
      ggplot2::scale_y_log10()

    if (Print) {
      print(g)
    }
    plots[[Type]] <- g
  }

  if ("MA plot" %in% Type) {
    # MA plot pre/after offset

    g1 <- ggplot2::ggplot(mapping = ggplot2::aes(x = log2(OffsetCorrected@MuBase_old), 
                                                 y = OffsetCorrected@MedianTau_old)) +
      # geom_hex() + 
      ggplot2::geom_point() +
      ggplot2::labs(
        x = "Mean expresssion (log2)",
        y = paste(
          "Log2 fold change", OffsetCorrected@GroupLabel1,
          "vs", OffsetCorrected@GroupLabel2),
        title = "Before correction") +
      viridis::scale_fill_viridis() +
      ggplot2::geom_hline(yintercept = 0, lty = 2) +
      ggplot2::geom_hline(yintercept = log2(OffsetCorrected@OffsetEst), lty = 1, col = "red")

    g2 <- ggplot2::ggplot(mapping = ggplot2::aes(x = log2(OffsetCorrected@MuBase), 
                                                 y = OffsetCorrected@MedianTau)) +
      # ggplot2::geom_hex() + 
      ggplot2::geom_point() +
      ggplot2::labs(
        x = "Mean expresssion (log2)",
        y = paste(
          "Log2 fold change", OffsetCorrected@GroupLabel1,
          "vs", OffsetCorrected@GroupLabel2),
        title = "After correction") +
      ggplot2::geom_hline(yintercept = 0, lty = 2) +
      viridis::scale_fill_viridis()

    g <- cowplot::plot_grid(g1, g2)
    if (Print) {
      print(g)
    }
    plots[[Type]] <- g
  }
  plots
}




DEPlots <- function(
    GroupLabel1,
    GroupLabel2,
    ProbThresholds = seq(0.5, 0.9995, by = 0.00025),
    Aux,
    OptThreshold,
    epsilon,
    EFDR,
    Table,
    Measure,
    Search = TRUE
    ) {

  plots <- list(
    MAPlot(Measure, Table, GroupLabel1, GroupLabel2, Epsilon),
    VolcanoPlot(Measure, Table, GroupLabel1, GroupLabel2, Epsilon)
  )

  if (Search) {
    plots <- c(
      list(PlotSearch(Measure, ProbThresholds, Aux, OptThreshold[[1]], EFDR)),
      plots
    )
  }

  cowplot::plot_grid(plotlist = plots, nrow = 1)
}



PlotSearch <- function(Measure, ProbThresholds, Aux, ProbThreshold, EFDR) {
  df <- data.frame(
    ProbThresholds, 
    EFDR = Aux$EFDRgrid, 
    EFNR = Aux$EFNRgrid)
  mdf <- reshape2::melt(df,measure.vars = c("EFDR", "EFNR"))
  ggplot2::ggplot(mdf, aes(x = ProbThresholds, y = value, color = variable)) +
    ggplot2::geom_line() +
    ggplot2::labs(
      y = "Error rate", 
      x = "Probability threshold",
      title = paste("Differential", MeasureName(Measure))
    ) +
    ggplot2::ylim(c(0, 1)) +
    ggplot2::scale_color_discrete(name = ""
      # , palette = "Dark2"
      ) +
    ggplot2::geom_hline(yintercept = EFDR, col = "steelblue", lty = 2) +
    ggplot2::geom_vline(xintercept = ProbThreshold, col = "indianred", lty = 1)
}



PlotGrids <- function(
    GroupLabel1,
    GroupLabel2,
    ProbThresholds = seq(0.5, 0.9995, by = 0.00025), 
    AuxMean, 
    AuxDisp, 
    AuxResDisp = NULL,
    OptThresholdM,
    OptThresholdD,
    OptThresholdE = NULL,
    EpsilonM,
    EpsilonD,
    EpsilonR,
    EFDR_M,
    EFDR_D,
    EFDR_R = NULL,
    TableMean,
    TableDisp,
    TableResDisp = NULL,
    IncludeEpsilon = !is.null(AuxResDisp)
    ) {

  plots <- list(
    PlotSearch("Mean", ProbThresholds, AuxMean, OptThresholdM[[1]], EFDR_M),
    PlotSearch("Disp", ProbThresholds, AuxDisp, OptThresholdD[[1]], EFDR_D)
  )

  if (IncludeEpsilon) {
    plots <- c(
      plots, 
      list(PlotSearch("ResDisp", ProbThresholds, AuxResDisp, OptThresholdE[[1]], EFDR_R))
    )
  }
  cowplot::plot_grid(plotlist = plots, nrow = 1)
}



PlotMAs <- function(
    GroupLabel1,
    GroupLabel2,
    ProbThresholds = seq(0.5, 0.9995, by = 0.00025), 
    AuxMean, 
    AuxDisp, 
    AuxResDisp = NULL,
    OptThresholdM,
    OptThresholdD,
    OptThresholdE = NULL,
    EpsilonM,
    EpsilonD,
    EpsilonR,
    EFDR_M,
    EFDR_D,
    EFDR_R = NULL,
    TableMean,
    TableDisp,
    TableResDisp = NULL,
    IncludeEpsilon = !is.null(AuxResDisp)
    ) {

  # MA plots
  plots <- list(
    MAPlot("Mean", TableMean, GroupLabel1, GroupLabel2, EpsilonM),
    MAPlot("Disp", TableDisp, GroupLabel1, GroupLabel2, EpsilonD)
  )

  if (IncludeEpsilon) {
    plots <- c(
      plots, 
      list(MAPlot("ResDisp", TableResDisp, GroupLabel1, GroupLabel2, EpsilonR))
    )
  }
  cowplot::plot_grid(plotlist = plots, nrow = 1)
}



PlotVolcanos <- function(
    GroupLabel1,
    GroupLabel2,
    ProbThresholds = seq(0.5, 0.9995, by = 0.00025), 
    AuxMean, 
    AuxDisp, 
    AuxResDisp = NULL,
    OptThresholdM,
    OptThresholdD,
    OptThresholdE = NULL,
    EpsilonM,
    EpsilonD,
    EpsilonR,
    EFDR_M,
    EFDR_D,
    EFDR_R = NULL,
    TableMean,
    TableDisp,
    TableResDisp = NULL,
    IncludeEpsilon = !is.null(AuxResDisp)
    ) {

  # MA plots
  plots <- list(
    VolcanoPlot("Mean", TableMean, GroupLabel1, GroupLabel2, EpsilonM),
    VolcanoPlot("Disp", TableDisp, GroupLabel1, GroupLabel2, EpsilonD)
  )

  if (IncludeEpsilon) {
    plots <- c(
      plots, 
      list(VolcanoPlot("ResDisp", TableResDisp, GroupLabel1, GroupLabel2, EpsilonR))
    )
  }
  cowplot::plot_grid(plotlist = plots, nrow = 1)

}

MAPlot <- function(Measure, Table, GroupLabel1, GroupLabel2, Epsilon) {
  IndDiff <- !Table[[paste0("ResultDiff", Measure)]] %in% c("ExcludedByUser", "NoDiff")
  ggplot2::ggplot(
      Table, 
      ggplot2::aes_string(
        x = paste0(Measure, "Overall"), 
        y = paste0(Measure, DistanceVar(Measure)))
    ) + 
    ggplot2::geom_hex(bins = 75, aes(fill = ..density..)) +
    ggplot2::geom_point() +
    ggplot2::scale_x_continuous(trans = "log2") +
    viridis::scale_fill_viridis(name = "Density") +
    ggplot2::labs(
      x = paste(cap(MeasureName(Measure)), " (log2)"),
      ylab = paste("Log2 fold change",
        GroupLabel1, "vs",
        GroupLabel2),
      title = paste("Differential", MeasureName(Measure))
    ) +
    ggplot2::geom_point(data = Table[IndDiff, ], shape = 16, colour = "violetred") +
    ggplot2::geom_hline(yintercept = -Epsilon, lty = "dashed", color = "grey40") +
    ggplot2::geom_hline(yintercept = Epsilon, lty = "dashed", color = "grey40")
}

VolcanoPlot <- function(Measure, Table, GroupLabel1, GroupLabel2, Epsilon) {
  IndDiff <- !Table[[paste0("ResultDiff", Measure)]] %in% c("ExcludedByUser", "NoDiff")
  ggplot2::ggplot(
      Table, 
      ggplot2::aes_string(
        x = paste0(Measure, DistanceVar(Measure)), 
        y = paste0("ProbDiff", Measure))
    ) +
    # ggplot2::geom_point() +
    ggplot2::geom_hex(bins = 75, aes(fill = ..density..)) +
    viridis::scale_fill_viridis(name = "Density") +
    ggplot2::ylim(c(0, 1)) +
    ggplot2::labs(y = "Posterior probability",
         x = paste("Log2 fold change", GroupLabel1, "vs", GroupLabel2),
         title = paste("Differential", MeasureName(Measure), "test")
    ) +
    ggplot2::geom_point(data = Table[IndDiff, ], shape = 16, col = "violetred") +
    ggplot2::geom_hline(yintercept = -Epsilon, lty = "dashed", color = "grey40") +
    ggplot2::geom_hline(yintercept = Epsilon, lty = "dashed", color = "grey40")
}

BASiCS_plotDE <- function(
    GroupLabel1,
    GroupLabel2,
    ProbThresholds = seq(0.5, 0.9995, by = 0.00025), 
    AuxMean, 
    AuxDisp, 
    AuxResDisp = NULL,
    OptThresholdM,
    OptThresholdD,
    OptThresholdE = NULL,
    EpsilonM,
    EpsilonD,
    EpsilonR,
    EFDR_M,
    EFDR_D,
    EFDR_R = NULL,
    TableMean,
    TableDisp,
    TableResDisp = NULL,
    IncludeEpsilon = !is.null(AuxResDisp),
    Search = TRUE) {

  if (Search) {

    ProbThresholds <- seq(0.5, 0.9995, by = 0.00025)

    plots <- list(
      PlotSearch("Mean", ProbThresholds, AuxMean, OptThresholdM[[1]], EFDR_M),
      PlotSearch("Disp", ProbThresholds, AuxDisp, OptThresholdD[[1]], EFDR_D)
    )

    if (IncludeEpsilon) {
      plots <- c(plots, list(PlotSearch("ResDisp", ProbThresholds, AuxResDisp, OptThresholdE[[1]], EFDR_R)))
    }
    cowplot::plot_grid(plotlist = plots, nrow = 1)
  }



  # MA plots
  plots <- list(
    MAPlot("Mean", TableMean, GroupLabel1, GroupLabel2, EpsilonM),
    MAPlot("Disp", TableDisp, GroupLabel1, GroupLabel2, EpsilonD)
  )

  if (IncludeEpsilon) {
    plots <- c(plots, list(MAPlot("ResDisp", TableResDisp, GroupLabel1, GroupLabel2, EpsilonR)))
  }
  cowplot::plot_grid(plotlist = plots, nrow = 1)





  plots <- list(
    VolcanoPlot("Mean", TableMean, GroupLabel1, GroupLabel2, EpsilonM),
    VolcanoPlot("Disp", TableDisp, GroupLabel1, GroupLabel2, EpsilonD)
  )

  if (IncludeEpsilon) {
    plots <- c(
      plots, 
      list(VolcanoPlot("ResDisp", TableResDisp, GroupLabel1, GroupLabel2, EpsilonR))
    )
  }
  cowplot::plot_grid(plotlist = plots, nrow = 1)
}
