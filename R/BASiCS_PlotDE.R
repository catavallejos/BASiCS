#' Produce plots assessing differential expression results
#' 
#' @param object A BASiCS_ResultsDE or BASiCS_ResultDE object.
#' @param Which Which plot to produce? Options: "MAPlot", "VolcanoPlot", 
#'  "GridPlot"
#' @param ... Passed to methods
#' @examples
#' data(ChainSC)
#' data(ChainRNA)
#'
#' Test <- BASiCS_TestDE(Chain1 = ChainSC, Chain2 = ChainRNA,
#'                       GroupLabel1 = 'SC', GroupLabel2 = 'P&S',
#'                       EpsilonM = log2(1.5), EpsilonD = log2(1.5),
#'                       OffSet = TRUE)
#' BASiCS_PlotDE(Test)
#' @return 
#'  Plot objects.
#' @author Catalina A. Vallejos \email{cnvallej@@uc.cl}
#' @author Nils Eling \email{eling@@ebi.ac.uk}
#' @author Alan O'Callaghan \email{a.b.o'callaghan@sms.ed.ac.uk}
#' @export
setMethod("BASiCS_PlotDE", signature(object = "BASiCS_ResultsDE"),
  function(object, Which = c("MAPlot", "VolcanoPlot", "GridPlot"), ...) {
    l <- lapply(object@Results, BASiCS_PlotDE, Which = Which, ...)
    if (length(Which) > 1) {
      nrow <- length(object@Results)
      labels <- sapply(object@Results, function(x) {
        cap(MeasureName(x@Name))
      })
      # labels <- Reduce(c, labels)
      labels <- c(labels, rep("", length(object@Results) * length(Which)))
      l <- lapply(l, 
        function(g) {
          g +
            ggplot2::theme(
              plot.margin = unit(c(0.100, 0, 0.05, 0), units = "npc")
            )
        }
      )
    } else {
      nrow <- length(Which)
      labels <- vapply(
        object@Results, 
        function(x) cap(MeasureName(x@Name)), 
        character(1)
      )
    }
    cowplot::plot_grid(plotlist = l, nrow = nrow, labels = labels, hjust = 0)
  }
)

#' @export
setMethod("BASiCS_PlotDE", signature(object = "BASiCS_ResultDE"),
  function(
    object,
    GroupLabel1 = object@GroupLabel1,
    GroupLabel2 = object@GroupLabel2,
    ProbThresholds = seq(0.5, 0.9995, by = 0.00025),
    Epsilon = object@Epsilon,
    EFDR = object@EFDR,
    Table = object@Table,
    Measure = object@Name,
    EFDRgrid = object@EFDRgrid,
    EFNRgrid = object@EFNRgrid,
    ProbThreshold = object@ProbThreshold,
    Which = c("MAPlot", "VolcanoPlot", "GridPlot")
  ){

    BASiCS_PlotDE(
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
setMethod("BASiCS_PlotDE", signature(object = "missing"),
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
    EFNR = EFNRgrid
  )
  mdf <- reshape2::melt(df,measure.vars = c("EFDR", "EFNR"))
  ggplot2::ggplot(mdf, 
      aes_string(x = "ProbThresholds", y = "value", color = "variable")
    ) +
    ggplot2::geom_line(na.rm = TRUE) +
    ggplot2::labs(
      y = "Error rate", 
      x = "Probability threshold"
      # ,title = paste("Differential", MeasureName(Measure))
    ) +
    ggplot2::ylim(c(0, 1)) +
    ggplot2::scale_color_brewer(name = "", palette = "Set2") +
    ggplot2::geom_hline(
      aes(yintercept = EFDR, colour = "Selected\nEFDR"),
      lty = 2,
      na.rm = TRUE
    ) +
    ggplot2::geom_vline(
      aes(colour = "Probability\nthreshold", xintercept = ProbThreshold),
      lty = 1,
      na.rm = TRUE
    )
}



VolcanoPlot <- function(Measure, Table, GroupLabel1, GroupLabel2, Epsilon) {
  IndDiff <- DiffExp(Table[[paste0("ResultDiff", Measure)]])
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
    ggplot2::geom_hex(bins = bins, aes_string(fill = "..density.."), na.rm = TRUE) +
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
    viridis::scale_fill_viridis(name = "Density", guide = FALSE) +
    ggplot2::labs(
      x = paste(cap(LogDistanceName(Measure)), GroupLabel1, "vs", GroupLabel2),
      y = "Posterior probability"
      # ,title = paste("Differential", MeasureName(Measure), "test")
    )
}

MAPlot <- function(Measure, Table, GroupLabel1, GroupLabel2, Epsilon) {

  IndDiff <- DiffExp(Table[[paste0("ResultDiff", Measure)]])
  # bins <- NClassFD2D(
  #   Table[[paste0(Measure, "Overall")]],
  #   Table[[paste0(Measure, DistanceVar(Measure))]]
  # )
  bins <- 100
  xscale <- ggplot2::scale_x_continuous(
    trans = if (Measure == "ResDisp") "identity" else "log2"
  )
  ggplot2::ggplot(
      Table, 
      ggplot2::aes_string(
        x = paste0(Measure, "Overall"), 
        y = paste0(Measure, DistanceVar(Measure)))
    ) + 
    ggplot2::geom_hex(bins = bins, aes_string(fill = "..density.."), na.rm = TRUE) +
    # ggplot2::geom_point() +
    ggplot2::geom_point(
      data = Table[IndDiff, ], 
      shape = 16, 
      colour = "violetred", 
      na.rm = TRUE) +
    ggplot2::geom_hline(
      yintercept = c(-Epsilon, Epsilon), 
      lty = "dashed", 
      color = "grey40") +
    xscale +
    viridis::scale_fill_viridis(name = "Density", guide = FALSE) +
    ggplot2::labs(
      x = paste(cap(MeasureName(Measure))),
      y = paste(cap(LogDistanceName(Measure)),
        GroupLabel1, "vs",
        GroupLabel2)
      # ,
      # title = paste("Differential", MeasureName(Measure))
    )
}

DiffExp <- function(res) {
  !res %in% c("ExcludedByUser", "ExcludedFromTesting", "NoDiff")
}
