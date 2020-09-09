#' Produce plots assessing differential expression results
#' 
#' @param object A \linkS4class{BASiCS_ResultsDE} or 
#' \linkS4class{BASiCS_ResultDE} object.
#' @param Plots Plots plot to produce? Options: "MA", "Volcano", 
#'  "Grid".
#' @param Parameters Character vector specifying the parameter(s) to produce 
#'  plots for,
#'  Available options are "Mean", (mu, mean expression),
#'  "Disp" (delta, overdispersion) and "ResDisp" 
#'  (epsilon, residual overdispersion).
#' @param MuX Use Mu (mean expression across both chains) as the X-axis for all
#'  MA plots? Default: TRUE.
#' @param Mu,GroupLabel1,GroupLabel2,ProbThresholds,Epsilon,EFDR,Table,Measure,EFDRgrid,EFNRgrid,ProbThreshold Internal arguments.
#' @param ... Passed to methods.
#' @examples
#' data(ChainSC)
#' data(ChainRNA)
#'
#' Test <- BASiCS_TestDE(Chain1 = ChainSC, Chain2 = ChainRNA,
#'                       GroupLabel1 = 'SC', GroupLabel2 = 'P&S',
#'                       EpsilonM = log2(1.5), EpsilonD = log2(1.5),
#'                       OffSet = TRUE)
#' BASiCS_PlotDE(Test)
#' @return A plot (possibly several combined using 
#' \code{\link[cowplot]{plot_grid}}).
#' 
#' @author Catalina A. Vallejos \email{cnvallej@@uc.cl}
#' @author Nils Eling \email{eling@@ebi.ac.uk}
#' @author Alan O'Callaghan \email{a.b.o'callaghan@sms.ed.ac.uk}
#' @rdname BASiCS_PlotDE
#' @export
setMethod("BASiCS_PlotDE", signature(object = "BASiCS_ResultsDE"),
  function(
      object,
      Plots = c("MA", "Volcano", "Grid"),
      Parameters = intersect(
        c("Mean", "Disp", "ResDisp"),
        names(object@Results)
      ),
      MuX = TRUE,
      ...) {

    Diff <- setdiff(Parameters, names(object@Results))
    if (length(Diff)) {
      stop(paste("Invalid Parameters selected:", Diff))
    }
    Mu <- if (MuX) object@Results$Mean@Table$MeanOverall else NULL
    l <- lapply(
      object@Results[Parameters],
      BASiCS_PlotDE,
      Plots = Plots,
      Mu = Mu,
      ...
    )
    if (length(Plots) > 1) {
      nrow <- length(object@Results[Parameters])
      labels <- vapply(object@Results[Parameters],
        function(x) .cap(.MeasureName(x@Name)),
        character(1)
      )
      # labels <- Reduce(c, labels)
      labels <- c(
        labels,
        rep("", length(object@Results[Parameters]) * length(Plots))
      )
    } else {
      nrow <- 1
      labels <- vapply(
        object@Results[Parameters],
        function(x) .cap(.MeasureName(x@Name)), 
        character(1)
      )
    }
    l <- lapply(l,
      function(g) {
        g + ggplot2::theme(
          plot.margin = grid::unit(c(0.1, 0, 0.05, 0), units = "npc")
        )
      }
    )
    if (length(l) > 1) {
      cowplot::plot_grid(plotlist = l, nrow = nrow, labels = labels, hjust = 0)
    } else {
      l[[1]]
    }
  }
)

#' @rdname BASiCS_PlotDE
#' @export
setMethod("BASiCS_PlotDE", signature(object = "BASiCS_ResultDE"),
  function(
    object,
    Plots = c("Grid", "MA", "Volcano"),
    Mu = NULL
  ){
    BASiCS_PlotDE(
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
      Mu = Mu,
      Plots = Plots
    )
  }
)
#' @rdname BASiCS_PlotDE
#' @export
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
    Mu,
    Plots = c("Grid", "MA", "Volcano")
  ) {

    Plots <- match.arg(Plots, several.ok = TRUE)
    PlotObjects <- list()

    if ("MA" %in% Plots) {
      PlotObjects <- c(PlotObjects, 
        list(
          .MAPlot(
            Measure,
            Table,
            GroupLabel1,
            GroupLabel2,
            Epsilon,
            Mu
          )
        )
      )
    }

    if ("Volcano" %in% Plots) {
      PlotObjects <- c(PlotObjects, 
        list(
          .VolcanoPlot(
            Measure,
            Table,
            GroupLabel1,
            GroupLabel2,
            Epsilon,
            ProbThreshold
          )
        )
      )
    }

    if ("Grid" %in% Plots) {
      PlotObjects <- c(
        list(
          .GridPlot(
            Measure,
            EFDR,
            ProbThresholds,
            EFDRgrid,
            EFNRgrid,
            ProbThreshold
          )
        ),
        PlotObjects
      )
    }
    if (length(PlotObjects) == 1) {
      PlotObjects <- PlotObjects[[1]]
    } else {
      PlotObjects <- cowplot::plot_grid(plotlist = PlotObjects, nrow = 1)
    }
    PlotObjects
  }
)



.GridPlot <- function(Measure, 
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
      # ,title = paste("Differential", .MeasureName(Measure))
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



.VolcanoPlot <- function(
    Measure,
    Table,
    GroupLabel1,
    GroupLabel2,
    Epsilon,
    ProbThreshold
  ) {

  dVar <- paste0(Measure, .DistanceVar(Measure))
  rVar <- paste0("ResultDiff", Measure)
  pVar <- paste0("ProbDiff", Measure)
  Table$IndDiff <- DiffExp(Table[[rVar]])
  # bins <- NClassFD2D(
  #   Table[[paste0(Measure, .DistanceVar(Measure))]],
  #   Table[[paste0("ProbDiff", Measure)]]
  # )
  Table <- Table[order(Table$IndDiff), ]
  bins <- 50
  ggplot2::ggplot(
      Table,
      ggplot2::aes_string(
        x = dVar,
        y = pVar)
    ) +
    ggplot2::geom_point(
      ggplot2::aes_string(color = rVar),
      shape = 16,
      alpha = 0.7
    ) +
    # ggplot2::geom_hex(
    #   bins = bins,
    #   aes_string(fill = "..density.."),
    #   na.rm = TRUE
    # ) +
    # ggplot2::geom_point(
    #   data = Table,
    #   shape = 16,
    #   col = "violetred",
    #   na.rm = TRUE,
    #   alpha = 0.8
    # ) +
    # viridis::scale_fill_viridis(name = "Density", guide = FALSE) +
    ggplot2::geom_vline(
      xintercept = c(-Epsilon, Epsilon),
      lty = "dashed",
      color = "grey40",
      na.rm = TRUE
    ) +
    ggplot2::geom_hline(
      yintercept = ProbThreshold, 
      lty = "dashed", 
      color = "grey40", 
      na.rm = TRUE
    ) +
    ggplot2::ylim(c(0, 1)) +
    # ggplot2::scale_color_brewer(palette = "Set1", name = NULL) +
    ggplot2::scale_color_manual(
      values = ColourMap(Table, dVar, rVar),
      drop = TRUE,
      name = NULL
    ) +
    ggplot2::labs(
      x = paste(
        .cap(.LogDistanceName(Measure)),
        GroupLabel1, "vs", GroupLabel2
      ),
      y = "Posterior probability"
    )
}

.MAPlot <- function(Measure, Table, GroupLabel1, GroupLabel2, Epsilon, Mu) {


  dVar <- paste0(Measure, .DistanceVar(Measure))
  rVar <- paste0("ResultDiff", Measure)
  Table$IndDiff <- DiffExp(Table[[rVar]])
  # bins <- NClassFD2D(
  #   Table[[paste0(Measure, "Overall")]],
  #   Table[[paste0(Measure, .DistanceVar(Measure))]]
  # )
  # bins <- 50
  xscale <- ggplot2::scale_x_continuous(
    trans = if (Measure == "ResDisp" & is.null(Mu)) "identity" else "log2"
  )
  Table$`_Mu` <- Mu
  Table <- Table[order(Table$IndDiff), ]
  ggplot2::ggplot(
      Table,
      # Table[!IndDiff, ],
      ggplot2::aes_string(
        # x = paste0(Measure, "Overall"), 
        x = if (!is.null(Mu)) "`_Mu`" else paste0(Measure, "Overall"),
        y = dVar
      )
    ) + 
    # ggplot2::geom_hex(
    #   bins = bins,
    #   aes_string(fill = "..density.."),
    #   na.rm = TRUE
    # ) +
    # ggplot2::geom_point(
    #   data = Table[IndDiff, ],
    #   shape = 16, 
    #   colour = "violetred", 
    #   na.rm = TRUE,
    #   alpha = 0.8
    # ) +
    # viridis::scale_fill_viridis(name = "Density", guide = FALSE) +
    ggplot2::geom_point(
      ggplot2::aes_string(color = rVar),
      shape = 16,
      alpha = 0.7
    ) +
    ggplot2::geom_hline(
      yintercept = c(-Epsilon, Epsilon), 
      lty = "dashed", 
      color = "grey40"
    ) +
    xscale +
    ggplot2::scale_color_manual(
      values = ColourMap(Table, dVar, rVar),
      drop = TRUE,
      name = NULL
    ) +
    ggplot2::labs(
      x = if (is.null(Mu)) paste(.cap(.MeasureName(Measure)))
        else "Mean expression",
      y = paste(.cap(.LogDistanceName(Measure)),
        GroupLabel1, "vs",
        GroupLabel2
      )
    )
}

DiffExp <- function(res) {
  !res %in% c(
    "ExcludedByUser", "ExcludedFromTesting", "ExcludedLowESS", "NoDiff"
  )
}

ColourMap <- function(Table, dVar, rVar) {
  colour_map <- c(
    "NoDiff" = "black",
    "ExcludedByUser" = "grey",
    "ExcludedLowESS" = "grey80",
    "ExcludedFromTesting" = "grey50"
  )
  ind_sig_up <- which(Table$IndDiff & Table[[dVar]] > 0)
  ind_sig_down <- which(Table$IndDiff & Table[[dVar]] < 0)
  if (length(ind_sig_up)) {
    up <- Table[[ind_sig_up[[1]], rVar]]
    colour_map[[up]] <- "firebrick"
  }
  if (length(ind_sig_down)) {
    down <- Table[[ind_sig_down[[1]], rVar]]
    colour_map[[down]] <- "dodgerblue"
  }
  colour_map
}
