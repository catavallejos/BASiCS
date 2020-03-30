HiddenHeaderBASiCS_Sim <- function(Mu,
                                   Mu_spikes,
                                   Delta,
                                   Phi, 
                                   S,
                                   Theta,
                                   BatchInfo) {
  
  # Data dimensions - number of cells
  n <- length(S)
  # Data dimensions - number of biological genes
  q.bio <- length(Mu)
  # Data dimensions - total number of genes
  if(!is.null(Mu_spikes)) {
    q <- q.bio + length(Mu_spikes)
  } else {
    q <- q.bio
  }
  # Data dimensions - number of batches
  if(!is.null(BatchInfo)) {
    nBatch <- length(unique(BatchInfo))  
  } else {
    nBatch <- 1
  }

  if (!all(c(n, nBatch, q, q.bio) > 0)) {
    stop("Arguments' dimensions are not valid")
  }
  
  # Basic checks for parameters that are always present
  if (!(is.vector(Mu) & is.numeric(Mu) & all(Mu > 0))) {
    stop("Invalid argument value for 'Mu'.")
  }
  if (!(is.vector(Delta) & is.numeric(Delta) & all(Delta >= 0) & 
        (length(Delta) == q.bio) )) {
    stop("Invalid argument value for 'Delta'.")
  }  
  if (!(is.vector(S) & is.numeric(S) & all(S > 0))) {
    stop("Invalid argument value for 'S'.")
  }
  if (!(is.numeric(Theta) & (length(Theta) == nBatch) & all(Theta >= 0))) {
    stop("Invalid argument value for 'Theta'.")
  }
  
  # Basic checks for parameters that are only present when spike-ins are there
  if (!is.null(Mu_spikes) &&
      (!(is.vector(Mu_spikes) & 
        is.numeric(Mu_spikes) & 
        all(Mu_spikes > 0)))) { 
    stop("Invalid argument value for 'Mu_spikes'.")
  }
  if (!is.null(Phi) &&
      (!(is.vector(Phi) & 
        is.numeric(Phi) & 
        all(Phi > 0) & 
        (length(Phi) == n) &
        isTRUE(all.equal(sum(Phi), n, tolerance = 0.01))))) {
    stop("Invalid argument value for 'Phi'.")
  }
  if (!is.null(BatchInfo) && !(length(BatchInfo) == n)) {
    stop("Invalid argument value for 'BatchInfo'. Must be same length as 'S'")
  }
  
  # Warnings when only some of the optional arguments are provided
  if (is.null(Mu_spikes) + is.null(Phi) == 1) {
    if (is.null(Mu_spikes)) {
      Phi <- NULL
      warning("Because 'Phi = NULL', 'Mu_spikes' will be ignored")
    }
    if (is.null(Phi)) {
      Mu_spikes = NULL
      warning("Because 'Mu_spikes = NULL', 'Phi' will be ignored")
    }
  }
  
  # Check that batches are present if spike-ins not included
  # Potentially remove this in next releases
  if (is.null(Mu_spikes)) {
    if (is.null(BatchInfo)) {
      stop("When spike-ins are not included, 'BatchInfo' is required.")
    } else {
      if (length(unique(BatchInfo)) <= 1) {
        stop(
          "When spike-ins are not included, 'BatchInfo' must contain",
          " multiple batches (i.e. length(unique(BatchInfo)) > 1)."
        )
      }
    }
  }
}
