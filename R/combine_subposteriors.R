#' @export
consensus_average <- function(...) {
  combine_subposteriors(..., method = "consensus")
}
#' @export
combine_subposteriors <- function(
    chains,
    gene_order = NULL,
    cell_order = NULL,
    method = c("consensus", "gaussian", "pie", "spde", "pconsensus", "consensusIndep"),
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
    chains[] <- mclapply(
      chains,
      function(chain) {
        offset_correct(reference_chain = reference_chain, chain = chain)
      },
      mc.cores = mc.cores
    )    
  }
  method <- match.arg(method)
  fun <- switch(method,
    "consensus" = weighted_posterior_average,
    "gaussian" = gaussian_approximation,
    "pie" = quantile_average,
    "spde" = spde,
    "pconsensus" = consensus, 
    "consensusIndep" = consensusIndep
  )

  params <- names(chains[[1]]@parameters)
  params <- setdiff(params, "RefFreq")

  mean_params <- lapply(
    params,
    iterate_chains,
    chains = chains,
    fun = fun,
#    mc.cores = mc.cores,
    subset_by = subset_by,
    ...
  )

  names(mean_params) <- params
  new("BASiCS_Chain",
    parameters = reorder_params(
      params = mean_params,
      gene_order = gene_order,
      cell_order = cell_order
    )
  )
}


weighted_posterior_average <- function(
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

quantile_average <- function(...) {
  weighted_posterior_average(..., sort = TRUE) 
}

iterate_chains <- function(param, chains, fun, ...) {
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
reorder_params <- function(params, gene_order = NULL, cell_order = NULL) {

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
