#' Run divide and conquer MCMC with BASiCS
#' 
#' @param Data SingleCellExperiemnt object
#' @param NSubsets Number of partitions to make
#' @param SubsetBy Passed to generateSubsets
#' @param Alpha Passed to generateSubsets
#' @param WithSpikes Similar to the argument in BASiCS_MCMC
#' @param mc.cores Number of cores for parallel processing
#' @param ... Passed to BASiCS_MCMC
#' 
#' 
#' @export
multi_MCMC <- function(
    Data,
    NSubsets = 5,
    SubsetBy = c("cell", "gene"),
    Alpha = 0.05,
    WithSpikes = FALSE,
    Regression = TRUE,
    mc.cores = min(detectCores(all.tests = FALSE, logical = TRUE), NSubsets),
    Seed = 42,
    Prior,
    MinGenesPerRBF = 100,
    ...
    ) {

  if ("package:RevoUtilsMath" %in% search() && getMKLthreads() > 1) {
    setMKLthreads(1)
    message("Setting matrix parallel threads to 1 with `setMKLthreads(1)`")
  }

  set.seed(Seed)
  ## These values should be passed in.
  L <- 12
  eta <- 5
  astwo <- 2
  bstwo <- 2
  PP <- BASiCS_PriorParam(Data, FixLocations = TRUE)
  Start <- BASiCS:::HiddenBASiCS_MCMC_Start(
    Data,
    Regression = TRUE,
    PriorParam = PP,
    WithSpikes = WithSpikes
  )

  lm <- log(Start$mu0)
  RBFLocations <- BASiCS:::.estimateRBFLocations(lm, L, FALSE)

  d <- (RBFLocations[[2]] - RBFLocations[[1]]) / 2
  retain <- vapply(
    RBFLocations,
    function(RBFLocation) {
      sum(lm > RBFLocation - d & lm < RBFLocation + d) > MinGenesPerRBF
    },
    logical(1)
  )

  # png("rbflocs.png")
  # plot(log(Start$mu0), log(Start$delta0), pch = 16, cex = 0.75, col = "grey80")
  # plot_distn <- function(mean, sd) {
  #   x <- seq(0, 10, length=200)
  #   y <- dnorm(x, mean=mean, sd=sd)
  #   lines(x, y, type="l", lwd=2)
  # }
  # abline(v = RBFLocations)
  # abline(v = RBFLocations + d, col = "grey60", lty = "dashed")
  # abline(v = RBFLocations - d, col = "grey60", lty = "dashed")
  # tmp <- sapply(RBFLocations, function(RBFLocation) plot_distn(mean = RBFLocation, sd = 0.5))
  # dev.off()

  PP$RBFLocations <- PP$RBFLocations[retain]
  retain_prior <- c(TRUE, TRUE, retain)
  PP$m <- PP$m[retain_prior]
  PP$V <- PP$V[retain_prior, retain_prior]
  PP$k <- length(PP$m)

  ## To do:
  ## Automatically determine NSubsets?
  ## Balance subsets for variables/batches?
  SubsetBy <- match.arg(SubsetBy)
  cat("Generating partitions...\n")
  Subsets <- generateSubsets(
    Data = Data,
    NSubsets = NSubsets,
    SubsetBy = SubsetBy,
    Alpha = Alpha,
    WithSpikes = WithSpikes
  )
  cat("Starting BASiCS...\n")
  output <- future.apply::future_lapply(
    Subsets,
    function(Subset) {
      PP$p.phi <- rep(1, ncol(Subset))
      PP$CellExponent <- ncol(Subset) / ncol(Data)
      PP$GeneExponent <- nrow(Subset) / nrow(Data)
      BASiCS_MCMC(
        Data = Subset,
        WithSpikes = WithSpikes,
        Regression = Regression,
        PriorParam = PP,
        ...
      )
    }
  )
  if (any(sapply(output, class) == "try-error")) {
    warning("Some MCMC runs failed!")
  }
  output
}


## Balance subsets for total NCounts (sample-wise) or quantile of mean count (gene-wise)
##   1. Divide genes/samples into quantile bins
##   2. Sample from each quantile
##   3. Check that anova(lm(count ~ batch)) is non-significant
##   4. If it is, back to step 2, else go ahead
#' Generate balanced subsets for divide and conquer BASiCS
#' 
#' Partitions data based on either cells or genes. Attempts to find a 
#' partitioning scheme which is "balanced" for either total reads per cell
#' across all genes (partitioning by gene)
#' or total expression per gene across all cells (partitioning by gene).
#' When partitioning by cell, at least 20 cells must be in each partition
#' or BASiCS_MCMC will fail.
#' If this partitioning fails, it will continue recursively up to a maximum 
#' number of iterations (20 by default).
#' 
#' @param Data a SingleCellExperiment object
#' @param NSubsets Number of partitions into which to divide Data 
#' @param SubsetBy Partition by "cell" or by "gene"
#' @param Alpha p-value threshold for ANOVA testing of "balance"
#' @param WithSpikes Similar to argument for BASiCS_MCMC - do the Data contain 
#'  spikes?
#' @param MaxDepth Maximum number of recursive
#' @param .depth Internal parameter. Do not set.
#' 
#' @return A list of SingleCellExperiment objects
generateSubsets <- function(
    Data,
    NSubsets,
    SubsetBy = c("cell", "gene"),
    Alpha = 0.05,
    WithSpikes = "spike-ins" %in% altExpNames(Data),
    MaxDepth = 20,
    .depth = 1
  ) {

  if (SubsetBy == "cell" && ((ncol(Data) / NSubsets) < 20)) {
    stop("Cannot generate subsets with at least 20 cells per subset!")
  }

  if (.depth > MaxDepth) {
    stop(paste("Unable to find balanced subset after", MaxDepth, "iterations!"))
  }

  SubsetBy <- match.arg(SubsetBy)

  if (SubsetBy == "cell") {
    balance_by <- colSums(counts(Data))
  } else {
    if (WithSpikes) {
      dimnames <- dimnames(Data)
      Spikes <- altExp(Data)
      BioData <- Data
    } else {
      Spikes <- NULL
      BioData <- Data
    }
    balance_by <- rowSums(counts(BioData))
  }

  ## How many quantiles?
  deciles <- quantile(
    balance_by,
    probs = seq(0, 1, length.out = 10)
  )
  bins <- cut(
    balance_by,
    breaks = deciles,
    include.lowest = TRUE
  )

  subsets <- as.character(bins)

  for (level in levels(bins)) {
    inds <- which(bins == level)
    subsets[inds] <- sample(seq_len(NSubsets), length(inds), replace = TRUE)
  }

  anova <- anova(lm(balance_by ~ subsets))
  balanced <- anova[["Pr(>F)"]][[1]] > Alpha &&
    ## scran::computeSumFactors fails if the condition below not met
    (SubsetBy == "gene" || all(table(subsets) > 20))

  if (balanced) {
    lapply(unique(subsets),
      function(subset) {
        if (SubsetBy == "cell") {
          ind <- subsets == subset
          Data <- Data[, ind]
          ## Sometimes spikes will be zero when subsampling.
          Data <- Data[rowSums(counts(Data)) != 0, ]
          Counts <- counts(Data)
          if (WithSpikes) {
            removeSpikes <- rowSums(assay(altExp(Data)) == 0)
            altExp(Data) <- altExp(Data)[!removeSpikes, ]
            SpikeInput <- Data@metadata$SpikeInput[removeSpikes, ]
            Tech <- isSpike(Data)
          } else {
            SpikeInput <- NULL
            Tech <- rep(FALSE, nrow(Counts))
          }
          BatchInfo <- Data@colData$BatchInfo
          suppressMessages(
            newBASiCS_Data(
              Counts = Counts, 
              Tech = Tech, 
              SpikeInfo = SpikeInput,
              BatchInfo = BatchInfo
            )
          )
        } else {
          DataOut <- BioData[subsets == subset, ]
          ind_keep <- colSums(counts(DataOut)) != 0
          DataOut <- DataOut[, ind_keep]
          if (WithSpikes) {
            Spikes <- Spikes[, ind_keep]
            Counts <- rbind(counts(DataOut), assay(Spikes))
            Tech <- c(rep(FALSE, nrow(DataOut)), rep(TRUE, nrow(Spikes)))
            SpikeInput <- Data@metadata$SpikeInput
          } else {
            Counts <- counts(DataOut)
            SpikeInput <- NULL
            Tech <- rep(FALSE, nrow(Counts))
          }
          BatchInfo <- DataOut@colData$BatchInfo
          # Filtered <- BASiCS_Filter(Counts = Counts, 
          #                           Tech = Tech, 
          #                           SpikeInput = SpikeInput, 
          #                           BatchInfo = BatchInfo)
          # if (!is.null(Filtered[["SpikeInput"]])) {
          #   Filtered[["SpikeInput"]] <- data.frame(
          #     rownames(Filtered[["Counts"]])[Filtered[["SpikeInput"]]],
          #     Filtered[["SpikeInput"]]
          #   )
          # }
          suppressMessages(
            newBASiCS_Data(
              Counts = Counts, 
              Tech = Tech, 
              SpikeInfo = SpikeInput,
              BatchInfo = BatchInfo
            )
          )
        }
      }
    )
  } else {
    generateSubsets(
      Data,
      NSubsets = NSubsets,
      SubsetBy = SubsetBy,
      Alpha = Alpha,
      WithSpikes = WithSpikes,
      .depth = .depth + 1
    )
  }
}
