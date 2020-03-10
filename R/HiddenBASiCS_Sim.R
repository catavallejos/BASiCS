HiddenBASiCS_Sim <- function(Mu, Mu_spikes = NULL, 
                             Delta, Phi = NULL, S, Theta) {
  
  # Number of cells
  n <- length(S)
  # Total number of genes, including biological and technical ones
  if (!is.null(Mu_spikes)) {
    q <- length(Mu) + length(Mu_spikes) 
    q.bio <- length(Mu)
    # Merge biological and technical genes
    Mu <- c(Mu, Mu_spikes)
  } else {
    q <- length(Mu)
    q.bio <- length(Mu)
  }
  if (is.null(Phi)) {
    Phi <- rep(1, times = n)
  }
  
  
  # Matrix where simulated counts will be stored
  Counts.sim <- matrix(0, ncol = n, nrow = q)
  # Matrix where gene-cell specific simulated random effects will be stored
  rho <- matrix(1, ncol = n, nrow = q.bio)
  # Simulated cell-specific random effects
  if (all(Theta > 0)) {
    Nu <- rgamma(n, shape = 1 / Theta, rate = 1 / (S * Theta))
  } else {
    Nu <- S
  }
  # Simulated counts data - biological genes
  for (i in seq_len(q.bio)) {
    if (Delta[i] > 0) {
      rho[i, ] <- rgamma(n, shape = 1 / Delta[i], rate = 1 / Delta[i])
      Counts.sim[i, ] <- rpois(n, lambda = Phi * Nu * rho[i, ] * Mu[i])
    }
  }
  ## BASiCS_Data will fail if a cell has no intrinsic counts
  while (any(ind_zero <- colSums(Counts.sim) == 0)) {
    for (i in seq_len(q.bio)) {
      if (Delta[i] > 0) {
        Counts.sim[, ind_zero] <- rpois(n,
          lambda = Phi[ind_zero] * Nu[ind_zero] * rho[i, ind_zero] * Mu[i]
        )
      }
    }
  }

  # Simulated counts data - spike-in genes
  if (!is.null(Mu_spikes)) {
    for (i in (seq_len(q - q.bio) + q.bio)) {
      Counts.sim[i, ] <- rpois(n, lambda = Nu * Mu[i])  
    }
  }
  
  rownames(Counts.sim) <- paste0("Gene", seq_len(q))
  
  if (!is.null(Mu_spikes)) {
    SpikeInfo <- data.frame(
      paste0("Gene", seq(q.bio + 1, q)),
      Mu[seq(q.bio + 1, q)]
    )
    Tech <- seq_len(q) > q.bio
  } else {
    SpikeInfo <- NULL
    Tech <- rep(FALSE, q)
  }
  list("Counts" = Counts.sim, "SpikeInfo" = SpikeInfo, "Tech" = Tech)
}
