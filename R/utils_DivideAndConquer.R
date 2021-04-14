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
#' @param NSubsets Integer specifying the number of batches into which to 
#' divide Data for divide and conquer inference.
#' @param SubsetBy Partition by "cell" or by "gene".
#' @param Alpha p-value threshold for ANOVA testing of "balance"
#' @param WithSpikes Similar to argument for BASiCS_MCMC - do the Data contain 
#'  spikes?
#' @param MaxDepth Maximum number of recursive
#' @param .Depth Internal parameter. Do not set.
#' 
#' @return A list of SingleCellExperiment objects
.generateSubsets <- function(
    Data,
    NSubsets,
    SubsetBy = c("cell", "gene"),
    Alpha = 0.05,
    WithSpikes = FALSE,
    MaxDepth = 20,
    .Depth = 1
  ) {

  if (SubsetBy == "cell" && ((ncol(Data) / NSubsets) < 20)) {
    stop("Cannot generate subsets with at least 20 cells per subset!")
  }

  if (.Depth > MaxDepth) {
    stop(paste("Unable to find balanced subset after", MaxDepth, "iterations!"))
  }

  SubsetBy <- match.arg(SubsetBy)

  if (SubsetBy == "cell") {
    balance_by <- colSums(counts(Data))
  } else {
    balance_by <- rowSums(counts(Data))
  }
  ## How many quantiles?
  nquantiles <- 10
  quantiles <- quantile(
    balance_by,
    probs = seq(0, 1, length.out = 10)
  )
  while (any(duplicated(quantiles))) {
    message(
      "Cannot find a balanced split with ",
      nquantiles,
      " quantiles, trying again with ",
      nquantiles - 1
    )
    nquantiles <- nquantiles - 1
    quantiles <- quantile(
      balance_by,
      probs = seq(0, 1, length.out = nquantiles)
    )
  }
  bins <- cut(
    balance_by,
    breaks = quantiles,
    include.lowest = TRUE
  )

  Subsets <- as.character(bins)

  for (level in levels(bins)) {
    inds <- which(bins == level)
    Subsets[inds] <- sample(seq_len(NSubsets), length(inds), replace = TRUE)
  }

  anova <- anova(lm(balance_by ~ Subsets))
  balanced <- anova[["Pr(>F)"]][[1]] > Alpha &&
    ## scran::computeSumFactors fails if the condition below not met
    (SubsetBy == "gene" || all(table(Subsets) > 20))

  if (balanced) {
    lapply(
      unique(Subsets),
      .subset,
      Data = Data,
      SubsetBy = SubsetBy,
      Subsets = Subsets,
      WithSpikes = WithSpikes
    )
  } else {
    .generateSubsets(
      Data,
      NSubsets = NSubsets,
      SubsetBy = SubsetBy,
      Alpha = Alpha,
      WithSpikes = WithSpikes,
      .Depth = .Depth + 1
    )
  }
}

.subset <- function(subset, Data, SubsetBy, Subsets, WithSpikes) {
  if (SubsetBy == "cell") {
    ind <- Subsets == subset
    Data <- Data[, ind]

    if (WithSpikes) {
      ## Sometimes spikes will be zero when subsampling.
      removeSpikes <- rowSums(assay(altExp(Data))) == 0
      altExp(Data) <- altExp(Data)[!removeSpikes, ]
    }
  }

  if (SubsetBy == "gene") {
    Data <- Data[Subsets == subset, ]
    ind_keep <- colSums(counts(Data)) != 0
    if (WithSpikes) {
      ind_keep <- ind_keep & (colSums(assay(altExp(Data))) != 0)
    }
    Data <- Data[, ind_keep]
  }
  Data
}

.consensus_average <- function(...) {
  .combine_subposteriors(..., CombineMethod = "consensus")
}

.combine_subposteriors <- function(
    Chains,
    GeneOrder = NULL,
    CellOrder = NULL,
    CombineMethod = c("consensus", "pie"),
    SubsetBy = c("gene", "cell"),
    BPPARAM = BiocParallel::bpparam(),
    ...
  ) {

  if (missing(SubsetBy)) {
    stop("Must supply SubsetBy")
  }
  SubsetBy <- match.arg(SubsetBy)

  if (SubsetBy == "cell") {
    ReferenceChain <- Chains[[1]]
    Chains[] <- bplapply(
      Chains,
      function(Chain) {
        .offset_correct(Chain = Chain, ReferenceChain = ReferenceChain)
      },
      BPPARAM = BPPARAM
    )    
  }
  CombineMethod <- match.arg(CombineMethod)
  Fun <- switch(CombineMethod,
    "consensus" = .weighted_posterior_average,
    "pie" = .quantile_average
  )

  Params <- names(Chains[[1]]@parameters)
  Params <- setdiff(Params, "RefFreq")

  mean_params <- bplapply(
    Params,
    .iterate_chains,
    Chains = Chains,
    Fun = Fun,
    SubsetBy = SubsetBy,
    BPPARAM = BPPARAM,
    ...
  )

  names(mean_params) <- Params
  new("BASiCS_Chain",
    parameters = .reorder_params(
      Params = mean_params,
      GeneOrder = GeneOrder,
      CellOrder = CellOrder
    )
  )
}


.weighted_posterior_average <- function(
    Colname,
    Chains,
    Param,
    Weighting = c("naive", "n_weight", "inverse_variance"),
    SubsetBy,
    Sort = FALSE
  ) {

  Weighting <- match.arg(Weighting)
  NSamples <- nrow(Chains[[1]]@parameters[[Param]])

  if (Weighting == "n_weight") {
    SubsetBy <- match.arg(SubsetBy, choices = c("cell", "gene"))
    if (SubsetBy == "cell") {
      n_param <- "s"
      all_colnames <- lapply(Chains, function(x) colnames(x@parameters[["s"]]))
    } else {
      n_param <- "mu"
      all_colnames <- lapply(Chains, function(x) colnames(x@parameters[["mu"]]))      
    }
    all_colnames <- Reduce(union, all_colnames)
    num <- length(all_colnames)      
  }

  weights <- numeric(length(Chains))
  subposterior_matrix <- matrix(
    NA,
    ncol = length(Chains),
    nrow = NSamples
  )
  for (i in seq_along(Chains)) {
    Chain <- Chains[[i]]
    mat <- Chain@parameters[[Param]]
    if (Colname %in% colnames(mat) || is.numeric(Colname)) {
      out <- mat[, Colname]
    } else {
      next
    }
    weights[[i]] <- switch(Weighting,
      "n_weight" = length(colnames(Chain@parameters[[n_param]])) / num, 
      "inverse_variance" = {
        ## what if infinite variance????
        ## Why infinite variance???? - Because invariant.
        ## Or NA (for epsilon)
        v <- 1 / var(out)
        if (!is.finite(v) & Param != "epsilon") {
          stop(paste0("Undefined variance for ", Param, " (", Colname, ")"))
        }
        v
      },
      "naive" = 1 / length(Chains)
    )
    if (!(
          (SubsetBy == "gene" && Param %in% c("mu", "delta", "epsilon")) ||
          (SubsetBy == "cell" && Param %in% c("phi", "nu", "s"))
        )) {
      out <- out * weights[[i]]
    }
    if (Param != "RefFreq" && !all(is.na(out))) {
      subposterior_matrix[, i] <- out
    }
  }
  ind_all_na <- apply(
    subposterior_matrix,
    2,
    function(col) all(is.na(col))
  )
  ## if all NA there's nothing fancy to do here
  ## this is the case where the only chain that has values are all missing
  if (all(ind_all_na)) {
    return(rep(NA, NSamples))
  }
  weights <- weights[!ind_all_na]
  subposterior_matrix <- subposterior_matrix[, !ind_all_na, drop = FALSE]

  if ((SubsetBy == "gene" && Param %in% c("mu", "delta", "epsilon")) ||
      (SubsetBy == "cell" && Param %in% c("phi", "nu", "s"))
      ) {
    sums <- subposterior_matrix 
    if (ncol(sums) > 1) {
      stop(
        "Too many draws for parameter ", Param, ", ", SubsetBy, ": ", Colname
      )
    }
  } else {
    if (Sort) {
      subposterior_matrix[] <- apply(
        subposterior_matrix[],
        2,
        sort
      )
    }
    sums <- rowSums(subposterior_matrix, na.rm = TRUE)
    sums <- sums / sum(weights, na.rm = TRUE)
  }
  sums
}

.quantile_average <- function(...) {
  .weighted_posterior_average(..., Sort = TRUE) 
}

.iterate_chains <- function(Param, Chains, Fun, ...) {
  message("Combining batch posteriors for ", Param, "\n")
  param_vals <- Chains[[1]]@parameters[[Param]]

  all_colnames <- lapply(
    Chains,
    function(Chain) {
      colnames(Chain@parameters[[Param]])
    }
  )

  if (all(vapply(all_colnames, is.null, logical(1)))) {
    all_colnames <- lapply(
      Chains, 
      function(Chain) {
        seq_len(ncol(Chain@parameters[[Param]]))
      }
    )
  }

  all_colnames <- Reduce(union, all_colnames)
  output <- matrix(
    NA,
    nrow = nrow(param_vals),
    ncol = length(all_colnames),
    dimnames = list(NULL, all_colnames)
  )
  output[, ] <- vapply(
    all_colnames,
    Fun,
    Chains = Chains,
    Param = Param,
    ...,
    numeric(nrow(output))
  )
  output
}

#' @export
.reorder_params <- function(Params, GeneOrder = NULL, CellOrder = NULL) {

  if (!is.null(GeneOrder)) {
    row_params <- intersect(names(Params), c("mu", "delta", "epsilon"))
    Params[row_params] <- lapply(
      Params[row_params],
      function(x) {
        x[, GeneOrder, drop = FALSE]
      }
    )
  }

  if (!is.null(CellOrder)) {
    col_params <- intersect(names(Params), c("s", "nu", "phi"))
    Params[col_params] <- lapply(
      Params[col_params],
      function(x) {
        x[, CellOrder, drop = FALSE]
      }
    )
  }

  Params
}

.offset_correct <- function(Chain, ReferenceChain) {
  offset <- matrixStats::rowSums2(Chain@parameters$mu) /
    matrixStats::rowSums2(ReferenceChain@parameters$mu)
  offset <- median(offset)

  Chain@parameters$mu <- Chain@parameters$mu / offset
  Chain
}


bplapply <- function(..., BPPARAM) {
  if (inherits(BPPARAM, "MulticoreParam")) {
    stop(
      "Cannot use BiocParallel::MultiCoreParam due to openMP code in BASiCS! ",
      "Try SnowParam instead."
    )
  }
  BiocParallel::bplapply(..., BPPARAM = BPPARAM)
}
