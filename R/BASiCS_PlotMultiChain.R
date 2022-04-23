#' todo document
#' todo rethink x/y axes plus colour
BASiCS_PlotMultiChain <- function(
        object,
        Param = "mu"
    ) {
    xlab <- if (Param %in% .GeneParams()) "Cene"
        else if (Param %in% .CellParams()) "Cell"
    mat <- sapply(object@chains,
        function(chain) {
            colMedians(chain@parameters[[Param]])
        }
    )
    means <- rowMeans(mat)
    range <- apply(mat, 1, range)
    rhat <- .Rhat(object, Param)
    ggplot() +
        aes(means, rhat) +
        ggplot2::geom_pointrange(xmin = range[1, ], xmax = range[2, ]) +
        geom_hline(yintercept = 1.05) +
        # viridis::scale_colour_viridis() +
        # ggplot2::labs(x = xlab, y = Param)
        ggplot2::labs(x = Param, y = "Rhat")
}
