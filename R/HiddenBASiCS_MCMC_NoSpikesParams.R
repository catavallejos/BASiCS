HiddenBASiCS_MCMC_NoSpikesParam <- function(Counts, ConstrainType,
                                            StochasticRef, q.bio, mu0,
                                            PriorDelta, ConstrainProp)
{
  if (PriorDelta == "gamma")
    stop("PriorDelta = 'gamma' is not supported for the no-spikes case")

  # 1: Full constrain; 2: Genes with average count >= 1
  if (ConstrainType == 1) {
    ConstrainGene <- seq_len(q.bio) - 1
    NotConstrainGene <- 0
    NonZero <- which(matrixStats::rowSums2(Counts) > 0)
  }
  if (ConstrainType == 2) {
    ConstrainGene <- which(matrixStats::rowSums2(Counts) >= 1) - 1
    NotConstrainGene <- which(matrixStats::rowSums2(Counts) < 1) - 1
    NonZero <- ConstrainGene + 1
  }
  Constrain <- mean(log(mu0[NonZero]))

  # Whether or not a stochatic reference is used
  # If stochastic, range of possible reference values only includes
  # the nearest 10% genes located around the constrain

  if (StochasticRef) {
    aux.ref <- cbind(ConstrainGene, abs(log(mu0[ConstrainGene+1]) - Constrain))
    aux.ref <- aux.ref[order(aux.ref[, 2]), ]
    # In total ConstrainProp*100% of genes to be used as reference candidates
    CandidateRef <- round(ConstrainProp * q.bio)
    # Fix for the code to run on the synthetic small dataset
    # generaed by makeExample_BASiCS function (less than 200 genes)
    if (length(ConstrainGene) > CandidateRef) {
      RefGenes <- aux.ref[seq_len(CandidateRef), 1]
    }
    else {
      RefGenes <- aux.ref[, 1]
    }
    RefGene <- RefGenes[1]
  }
  else {
    aux.ref <- which(abs(log(mu0[ConstrainGene+1]) - Constrain) ==
                       min(abs(log(mu0[ConstrainGene+1]) - Constrain)))[1]
    RefGene <- ConstrainGene[aux.ref]
    RefGenes <- RefGene
  }

  list(ConstrainGene = ConstrainGene, NotConstrainGene = NotConstrainGene,
       Constrain = Constrain, RefGenes = RefGenes, RefGene = RefGene)
}
