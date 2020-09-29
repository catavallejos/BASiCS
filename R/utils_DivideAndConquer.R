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
BASiCS_DivideAndConquer <- function(
    Data,
    NSubsets = 5,
    SubsetBy = c("cell", "gene"),
    Alpha = 0.05,
    WithSpikes = FALSE,
    Regression = TRUE,
    mc.cores = min(detectCores(all.tests = FALSE, logical = TRUE), NSubsets),
    ...
    ) {

  ## These values should be passed in.
  L <- 12
  eta <- 5
  astwo <- 2
  bstwo <- 2
  PP <- BASiCS_PriorParam(
    Data,
    PriorMu = "EmpiricalBayes"
  )
  start <- BASiCS:::.BASiCS_MCMC_Start(
    Data,
    PriorParam = PP,
    Regression = TRUE,
    WithSpikes = WithSpikes
  )

  Locations <- BASiCS:::.estimateRBFLocations(
    start$mu0,
    L,
    RBFMinMax = FALSE
  )
  PP$RBFLocations <- Locations

  ## To do:
  ## Automatically determine NSubsets?
  ## Balance subsets for variables/batches?
  SubsetBy <- match.arg(SubsetBy)
  cat("Generating partitions...\n")
  Subsets <- .generateSubsets(
    Data = Data,
    NSubsets = NSubsets,
    SubsetBy = SubsetBy,
    Alpha = Alpha,
    WithSpikes = WithSpikes
  )
  cat("Starting BASiCS...\n")
  output <- lapply(
    Subsets,
    function(Subset) {
      if (SubsetBy == "gene") {
        ind <- match(rownames(Subset), rownames(Data))
        PP$mu.mu <- PP$mu.mu[ind]
      } else if (SubsetBy == "cell") {
        ind <- match(colnames(Subset), colnames(Data))
        PP$p.phi <- PP$p.phi[ind]
      }
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
.generateSubsets <- function(
    Data,
    NSubsets,
    SubsetBy = c("cell", "gene"),
    Alpha = 0.05,
    WithSpikes = FALSE,
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
    lapply(unique(Subsets),
      .subset,
      Data = Data,
      SubsetBy = SubsetBy,
      Subsets = Subsets
    )
  } else {
    .generateSubsets(
      Data,
      NSubsets = NSubsets,
      SubsetBy = SubsetBy,
      Alpha = Alpha,
      WithSpikes = WithSpikes,
      .depth = .depth + 1
    )
  }
}

.subset <- function(subset, Data, SubsetBy, Subsets) {
  if (SubsetBy == "cell") {
    ind <- Subsets == subset
    Data <- Data[, ind]

    ## Sometimes spikes will be zero when subsampling.
    removeSpikes <- rowSums(assay(altExp(Data)) == 0)
    altExp(Data) <- altExp(Data)[!removeSpikes, ]
    Data <- Data[rowSums(counts(Data)) != 0, ]
  }

  if (SubsetBy == "gene") {
    Data <- Data[Subsets == subset, ]
    ind_keep <- colSums(counts(Data)) != 0
    Data <- Data[, ind_keep]
  }
  Data
}

.consensus_average <- function(...) {
  .combine_subposteriors(..., method = "consensus")
}

.combine_subposteriors <- function(
    chains,
    gene_order = NULL,
    cell_order = NULL,
    method = c("consensus", "pie"),
    subset_by = c("gene", "cell"),
    mc.cores = min(detectCores(all.tests = FALSE, logical = TRUE), 16),
    ...
  ) {

  if (missing(subset_by)) {
    stop("Must supply subset_by")
  }
  subset_by <- match.arg(subset_by)

  if (subset_by == "cell") {
    reference_chain <- chains[[1]]
    chains[] <- lapply(
      chains,
      function(chain) {
        .offset_correct(chain = chain, reference_chain = reference_chain)
      }
      # , mc.cores = mc.cores
    )    
  }
  method <- match.arg(method)
  fun <- switch(method,
    "consensus" = .weighted_posterior_average,
    "pie" = .quantile_average
  )

  params <- names(chains[[1]]@parameters)
  params <- setdiff(params, "RefFreq")

  mean_params <- lapply(
    params,
    .iterate_chains,
    chains = chains,
    fun = fun,
#    mc.cores = mc.cores,
    subset_by = subset_by,
    ...
  )

  names(mean_params) <- params
  new("BASiCS_Chain",
    parameters = .reorder_params(
      params = mean_params,
      gene_order = gene_order,
      cell_order = cell_order
    )
  )
}


.weighted_posterior_average <- function(
    colname,
    chains,
    param,
    weight_method = c("naive", "n_weight", "inverse_variance"),
    subset_by,
    sort = FALSE
  ) {

  weight_method <- match.arg(weight_method)

  if (weight_method == "n_weight") {
    subset_by <- match.arg(subset_by, choices = c("cell", "gene"))
    if (subset_by == "cell") {
      n_param <- "s"
      all_colnames <- lapply(chains, function(x) colnames(x@parameters[["s"]]))
    } else {
      n_param <- "mu"
      all_colnames <- lapply(chains, function(x) colnames(x@parameters[["mu"]]))      
    }
    all_colnames <- Reduce(union, all_colnames)
    num <- length(all_colnames)      
  }

  weights <- numeric(length(chains))
  subposterior_matrix <- matrix(
    NA,
    ncol = length(chains),
    nrow = nrow(chains[[1]]@parameters[[param]])
  )
  for (i in seq_along(chains)) {
    chain <- chains[[i]]
    mat <- chain@parameters[[param]]
    if (colname %in% colnames(mat) || is.numeric(colname)) {
      out <- mat[, colname]
    } else {
      next
    }
    weights[[i]] <- switch(weight_method,
      "n_weight" = length(colnames(chain@parameters[[n_param]])) / num, 
      "inverse_variance" = {
        ## what if infinite variance????
        ## Why infinite variance???? - Because invariant.
        ## Or NA (for epsilon)
        v <- 1 / var(out)
        if (!is.finite(v) & param != "epsilon") {
          stop(paste0("Undefined variance for ", param, " (", colname, ")"))
        }
        v
      },
      "naive" = 1 / length(chains)
    )
    if (!(
          (subset_by == "gene" && param %in% c("mu", "delta", "epsilon")) ||
          (subset_by == "cell" && param %in% c("phi", "nu", "s"))
        )) {
      out <- out * weights[[i]]
    }
    if (param != "RefFreq" && !all(is.na(out))) {
      subposterior_matrix[, i] <- out
    }
  }
  ind_all_na <- apply(
    subposterior_matrix,
    2,
    function(col) all(is.na(col))
  )

  weights <- weights[!ind_all_na]
  subposterior_matrix <- subposterior_matrix[, !ind_all_na, drop = FALSE]

  if ((subset_by == "gene" && param %in% c("mu", "delta", "epsilon")) ||
      (subset_by == "cell" && param %in% c("phi", "nu", "s"))
      ) {
    sums <- subposterior_matrix 
  } else {
    if (sort) {
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
  .weighted_posterior_average(..., sort = TRUE) 
}

.iterate_chains <- function(param, chains, fun, ...) {
  message("Combining batch posteriors for ", param, " parameter\n")
  param_vals <- chains[[1]]@parameters[[param]]

  all_colnames <- lapply(
    chains,
    function(chain) {
      colnames(chain@parameters[[param]])
    }
  )

  if (all(sapply(all_colnames, is.null))) {
    all_colnames <- lapply(
      chains, 
      function(chain) {
        seq_len(ncol(chain@parameters[[param]]))
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

  output[, ] <- sapply(
    all_colnames,
    fun,
    chains = chains,
    param = param,
    ...
  )
  output
}

#' @export
.reorder_params <- function(params, gene_order = NULL, cell_order = NULL) {

  if (!is.null(gene_order)) {
    row_params <- intersect(names(params), c("mu", "delta", "epsilon"))
    params[row_params] <- lapply(
      params[row_params],
      function(x) {
        x[, gene_order, drop = FALSE]
      }
    )
  }

  if (!is.null(cell_order)) {
    col_params <- intersect(names(params), c("s", "nu", "phi"))
    params[col_params] <- lapply(
      params[col_params],
      function(x) {
        x[, cell_order, drop = FALSE]
      }
    )
  }

  params
}

.offset_correct <- function(chain, reference_chain) {
  offset <- matrixStats::rowSums2(chain@parameters$mu) /
    matrixStats::rowSums2(reference_chain@parameters$mu)
  offset <- median(offset)

  chain@parameters$mu <- chain@parameters$mu / offset
  chain
}
