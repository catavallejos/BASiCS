#' Plot variance decomposition results.
#' @param Decomp The output of \code{\link{BASiCS_VarianceDecomp}}.
#' @param beside If \code{TRUE}, bars are placed beside each other.
#' If \code{FALSE}, bars are stacked.
#' @param nBatch Number of batches.
#' @param main Plot title.
#' @param xlabs x-axis labels. Defaults to "Batch 1", "Batch 2", etc.
#' @param ylab y axis label.
#' @return A ggplot object.
BASiCS_PlotVarianceDecomp <- function(
    Decomp,
    beside = FALSE,
    nBatch = ((ncol(Decomp) - 2) / 3) - 1,
    main = "Overall variance decomposition",
    xlabs = if (nBatch == 1) "Overall"
      else c(
        "Overall",
        paste("Batch", seq_len(nBatch))
      ),
    ylab = "% of variance"
  ) {
  outmat <- 100 * matrix(
    apply(Decomp[, -c(1,2)], 2, mean),
    nrow = 3, byrow = FALSE
  )
  rownames(outmat) <- c("Shot noise", "Technical", "Biological")
  colnames(outmat) <- xlabs
  mdf <- reshape2::melt(outmat)
  ggplot(mdf,
      aes(x = Var2, y = value, fill = Var1)
    ) +
    geom_col(position = if (!beside) "fill" else "dodge") +
    scale_fill_brewer(palette = "Set1", name = NULL) +
    labs(main = main, x = NULL, y = ylab)
}
