#' @rdname BASiCS_ShowFit
#'
#' @aliases BASiCS_showFit
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
#' @param Uncertainty logical indicator. If true, statistical uncertainty
#' around the regression fit is shown in the plot. 
#'
#' @examples
#' data(ChainRNAReg)
#' BASiCS_ShowFit(ChainRNAReg)
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
BASiCS_ShowFit <- function(object,
                           xlab = "log(mu)",
                           ylab = "log(delta)",
                           pch = 16,
                           smooth = TRUE,
                           variance = 1.2,
                           colour = "dark blue",
                           markExcludedGenes = TRUE,
                           GenesSel = NULL,
                           colourGenesSel = "dark red",
                           Uncertainty = TRUE) {

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
  if (is.null(myu <- object@parameters$RBFLocations)) {
    myu <- seq(min(grid.mu), by = range/(n-3), length.out = n-2)
  }
  h <- diff(myu)*variance

  B <- matrix(1, length(grid.mu), n)
  B[, 2] <- grid.mu
  for (j in seq_len(n - 2)) {
    B[, j+2] = exp(-0.5 * (grid.mu - myu[j])^2 / (h[1]^2))
  }

  # Calculate yhat = X*beta
  yhat <- apply(object@parameters$beta, 1, function(n) B%*%n)
  yhat.HPD <- coda::HPDinterval(coda::mcmc(t(yhat)), 0.95)

  df <- data.frame(mu = log(colMedians(object@parameters$mu)),
                   delta = log(colMedians(object@parameters$delta)),
                   included = !is.na(object@parameters$epsilon[1, ]))
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
      ggplot2::scale_fill_gradientn(
        name = "",
        colours = grDevices::colorRampPalette(cols)(100),
        guide = "none"
      )
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
                       colour = "dark red") #
  if(Uncertainty == TRUE) {
    plot.out <- plot.out + geom_ribbon(
      data = df2,
      inherit.aes = FALSE,
      mapping = aes_string(x = "mu2",
      ymin = "yhat.lower",
      ymax = "yhat.upper"),
      alpha = 0.5
    )
  }
    
  return(plot.out)
}

#' @export
BASiCS_showFit <- function(...) {
  .Deprecated("BASiCS_ShowFit")
  BASiCS_ShowFit(...)
}

#' todo document
BASiCS_ShowMultiFit <- function(
    object,
    variance = 1.2,
    Uncertainty = TRUE, ## HPD of curve
    Variability = TRUE ## chain-to-chain variability
  ) {

  if (!inherits(object, "BASiCS_MultiChain")) {
    stop(paste0("Incorrect class for object:", class(object)))
  }

  if (!("beta" %in% names(object@chains[[1]]@parameters))) {
    stop("'beta' is missing. Regression was not performed.")
  }
  mu_mat <- sapply(object@chains,
    function(chain) {
      colMedians(chain@parameters$mu)
    }
  )
  delta_mat <- sapply(object@chains,
    function(chain) {
      colMedians(chain@parameters$delta)
    }
  )
  mu <- rowMeans(mu_mat)
  delta <- rowMeans(delta_mat)
  mu_r <- log(apply(mu_mat, 1, range))
  delta_r <- log(apply(delta_mat, 1, range))

  geoms <- lapply(seq_along(object@chains),
    function(i) {
      chain <- object@chains[[i]]
      m <- log(chain@parameters$mu[1, ])
      grid.mu <- seq(round(min(m), digits = 2),
                     round(max(m), digits = 2), length.out = 1000)

      # Create design matrix across the grid
      n <- ncol(chain@parameters$beta)
      range <- diff(range(grid.mu))
      if (is.null(myu <- chain@parameters$RBFLocations)) {
        myu <- seq(min(grid.mu), by = range/(n-3), length.out = n-2)
      }
      h <- diff(myu) * variance

      B <- matrix(1, length(grid.mu), n)
      B[, 2] <- grid.mu
      for (j in seq_len(n - 2)) {
        B[, j+2] = exp(-0.5 * (grid.mu - myu[j])^2 / (h[1]^2))
      }

      # Calculate yhat = X*beta
      yhat <- apply(chain@parameters$beta, 1, function(n) B%*%n)
      yhat.HPD <- coda::HPDinterval(coda::mcmc(t(yhat)), 0.95)

      df2 <- data.frame(
        mu2 = grid.mu,
        yhat = rowMedians(yhat),
        yhat.upper = yhat.HPD[, 2],
        yhat.lower = yhat.HPD[, 1],
        Chain = as.character(i)
      )
      geom <- list(
        ggplot2::geom_line(
          data = df2,
          inherit.aes = FALSE,
          mapping = ggplot2::aes_string(
            x = "mu2", y = "yhat",
            colour = "Chain"
          )
        )
      )
      if (Uncertainty) {
        geom <- c(geom,
          list(
            geom_ribbon(
              data = df2,
              inherit.aes = FALSE,
              mapping = aes_string(
                x = "mu2",
                ymin = "yhat.lower",
                ymax = "yhat.upper",
                fill = "Chain"
              ),
              alpha = 0.25
            )
          )
        )
      }
      geom
    }
  )
  g <- ggplot() +
    aes(log(mu), log(delta))
  if (Variability) {
    g <- g + 
      ggplot2::geom_pointrange(aes(xmin = mu_r[1, ], xmax = mu_r[2, ])) +
      ggplot2::geom_pointrange(aes(ymin = delta_r[1, ], ymax = delta_r[2, ]))
  } else {
    g <- g + geom_point()
  }
  g +
    geoms +
    ggplot2::scale_colour_brewer(palette = "Set1", aesthetics = c("fill", "colour"))
}
