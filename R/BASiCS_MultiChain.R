BASiCS_MultiMCMC <- function(
        ...,
        NChains = 4,
        Seeds = seq_len(NChains),
        BPPARAM = BiocParallel::bpparam()
    ) {
    Chains <- bplapply(seq_len(NChains),
        function(i) {
            set.seed(Seeds[[i]])
            BASiCS_MCMC(...)
        },
        BPPARAM = BPPARAM
    )
    new(
        "BASiCS_MultiChain",
        chains = Chains 
    )
}
