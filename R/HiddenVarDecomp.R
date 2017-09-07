# Used in BASiCS_VarianceDecomp
HiddenVarDecomp <- function(Chain) {
    if (!is(Chain, "BASiCS_Chain")) 
        stop("'Chain' is not a BASiCS_Chain class object.")
    
    N = nrow(Chain@delta)
    q.bio = ncol(Chain@delta)
    UniqueBatch = colnames(Chain@theta)
    nBatch = length(UniqueBatch)
    
    if (nBatch > 1) {
        Theta = apply(Chain@theta, 1, median)
    } else {
        Theta = as.vector(Chain@theta)
    }
    
    # To store global values (uses median values across all cells)
    PhiS = apply(Chain@phi * Chain@s, 1, median)
    Aux = (1/(PhiS * Chain@mu[, 1:q.bio])) + Chain@delta * (Theta + 1)
    TechVarGlobal = Theta/(Aux + Theta)
    BioVarGlobal = (Chain@delta * (Theta + 1))/(Aux + Theta)
    
    # To store batch specific values (in arrays)
    TechVarBatch = array(0, dim = c(N, q.bio, nBatch))  # Technical
    BioVarBatch = array(0, dim = c(N, q.bio, nBatch))  # Biological
    
    if (nBatch > 1) {
        for (Batch in 1:nBatch) {
            PhiSBatch = apply(Chain@phi[, grep(UniqueBatch[Batch], colnames(Chain@phi))] * Chain@s[, grep(UniqueBatch[Batch], 
                colnames(Chain@phi))], 1, median)
            Aux = (1/(PhiSBatch * Chain@mu[, 1:q.bio])) + Chain@delta * (Chain@theta[, Batch] + 1)
            TechVarBatch[, , Batch] = Chain@theta[, Batch]/(Aux + Chain@theta[, Batch])
            BioVarBatch[, , Batch] = (Chain@delta * (Chain@theta[, Batch] + 1))/(Aux + Chain@theta[, Batch])
        }
    }
    
    if (nBatch > 1) {
        list(TechVarGlobal = TechVarGlobal, BioVarGlobal = BioVarGlobal, TechVarBatch = TechVarBatch, BioVarBatch = BioVarBatch)
    } else {
        list(TechVarGlobal = TechVarGlobal, BioVarGlobal = BioVarGlobal)
    }
}
