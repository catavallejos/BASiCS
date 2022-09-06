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
#' @param TransLogit Logical scalar controlling whether a logit transform
#' is applied to the posterior probability in the y-axis of volcano plots.
#' As logit(0) and logit(1) are undefined, we clip these values near the range
#' of the data excluding 0 and 1.
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
#' @author Alan O'Callaghan
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
    Mu = NULL,
    TransLogit = FALSE
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
      Plots = Plots,
      TransLogit = TransLogit
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
    TransLogit = FALSE,
    Plots = c("Grid", "MA", "Volcano")
  ) {

    Plots <- match.arg(Plots, several.ok = TRUE)
    PlotObjects <- list()

    if ("MA" %in% Plots) {
      PlotObjects <- c(PlotObjects, 
        list(
          .MAPlot(
            Measure = Measure,
            Table = Table,
            GroupLabel1 = GroupLabel1,
            GroupLabel2 = GroupLabel2,
            Epsilon = Epsilon,
            Mu = Mu
          )
        )
      )
    }

    if ("Volcano" %in% Plots) {
      PlotObjects <- c(PlotObjects, 
        list(
          .VolcanoPlot(
            Measure = Measure,
            Table = Table,
            GroupLabel1 = GroupLabel1,
            GroupLabel2 = GroupLabel2,
            Epsilon = Epsilon,
            ProbThreshold = ProbThreshold,
            TransLogit = TransLogit
          )
        )
      )
    }

    if ("Grid" %in% Plots) {
      PlotObjects <- c(
        list(
          .GridPlot(
            EFDR = EFDR,
            ProbThresholds = ProbThresholds,
            EFDRgrid = EFDRgrid,
            EFNRgrid = EFNRgrid,
            ProbThreshold = ProbThreshold
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

.GridPlot <- function(
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
  mdf <- reshape2::melt(df, measure.vars = c("EFDR", "EFNR"))
  ggplot2::ggplot(mdf, 
      aes_string(x = "ProbThresholds", y = "value", colour = "variable")
    ) +
    ggplot2::geom_line(na.rm = TRUE) +
    ggplot2::labs(
      y = "Error rate", 
      x = "Probability threshold"
    ) +
    ggplot2::ylim(c(0, 1)) +
    ggplot2::scale_colour_brewer(name = "", palette = "Set2") +
    ggplot2::geom_hline(
      aes(yintercept = EFDR, colour = "Selected\nEFDR"),
      linetype = "dashed",
      na.rm = TRUE
    ) +
    ggplot2::geom_vline(
      aes(colour = "Probability\nthreshold", xintercept = ProbThreshold),
      linetype = "dashed",
      na.rm = TRUE
    ) +
    ggplot2::theme_bw()
}



.VolcanoPlot <- function(
    Measure,
    Table,
    GroupLabel1,
    GroupLabel2,
    Epsilon,
    ProbThreshold,
    TransLogit = FALSE
  ) {

  dVar <- paste0(Measure, .DistanceVar(Measure))
  rVar <- paste0("ResultDiff", Measure)
  pVar <- paste0("ProbDiff", Measure)
  Table$IndDiff <- .DiffExp(Table[[rVar]])
  Table <- Table[order(Table$IndDiff), ]
  bins <- 50
  if (TransLogit) {
    Table[[pVar]] <- .logit_nudge(Table[[pVar]])
  }
  g <- ggplot2::ggplot(
      Table,
      ggplot2::aes_string(
        x = dVar,
        y = pVar
      )
    ) +
    ggplot2::geom_point(
      ggplot2::aes_string(colour = rVar),
      shape = 16,
      alpha = 0.6
    ) +
    ggplot2::geom_vline(
      xintercept = c(-Epsilon, Epsilon),
      lty = "dashed",
      colour = "grey40",
      na.rm = TRUE
    ) +
    ggplot2::geom_hline(
      yintercept = ProbThreshold,
      lty = "dashed",
      colour = "grey40",
      na.rm = TRUE
    ) +
    ggplot2::scale_colour_manual(
      values = .ColourMap(Table, dVar, rVar),
      drop = TRUE,
      name = NULL
    ) +
    ggplot2::labs(
      x = paste(
        .cap(.LogDistanceName(Measure)),
        GroupLabel1, "vs", GroupLabel2
      ),
      y = "Posterior probability"
    ) +
    ggplot2::theme_bw() +
    if (TransLogit) {
      ggplot2::scale_y_continuous(trans = "logit",
        # limits = c(.Machine$double.xmin, 1 - .Machine$double.xmin)
        # limits = c(0.01, 0.99)
        breaks = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99)
      )
    } else {
      ggplot2::ylim(c(0, 1))
    }
  g
}

## Can't have logit(0) or logit(1) as these are -Inf, Inf
.logit_nudge <- function(x) {
  x2 <- x
  x[x==0] <- min(x[x != 0], na.rm = TRUE) - .Machine$double.xmin
  x[x==1] <- max(x[x != 1], na.rm = TRUE) + .Machine$double.xmin
  x
}

.MAPlot <- function(Measure, Table, GroupLabel1, GroupLabel2, Epsilon, Mu) {

  dVar <- paste0(Measure, .DistanceVar(Measure))
  rVar <- paste0("ResultDiff", Measure)
  Table$IndDiff <- .DiffExp(Table[[rVar]])
  xscale <- ggplot2::scale_x_continuous(
    trans = if (Measure == "ResDisp" & is.null(Mu)) "identity" else "log2"
  )
  Table$`_Mu` <- Mu
  Table <- Table[order(Table$IndDiff), ]
  ggplot2::ggplot(
      Table,
      ggplot2::aes_string(
        x = if (!is.null(Mu)) "`_Mu`" else paste0(Measure, "Overall"),
        y = dVar
      )
    ) + 
    ggplot2::geom_point(
      ggplot2::aes_string(colour = rVar),
      shape = 16,
      alpha = 0.6
    ) +
    ggplot2::geom_hline(
      yintercept = c(-Epsilon, Epsilon),
      lty = "dashed", 
      colour = "grey40"
    ) +
    xscale +
    ggplot2::scale_colour_manual(
      values = .ColourMap(Table, dVar, rVar),
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
    ) +
    ggplot2::theme_bw()
}

.DiffExp <- function(res) {
  !res %in% c(
    "ExcludedByUser", "ExcludedFromTesting", "ExcludedLowESS", "NoDiff",
    "ExcludedSmallDiff"
  )
}

.ColourMap <- function(Table, dVar, rVar) {
  colour_map <- c(
    "NoDiff" = "black",
    "ExcludedByUser" = "grey",
    "ExcludedLowESS" = "grey80",
    "ExcludedFromTesting" = "grey50",
    "ExcludedSmallDiff" = "grey30"
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
