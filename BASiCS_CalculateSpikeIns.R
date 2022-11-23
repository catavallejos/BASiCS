#' Convert concentration in moles per microlitre to molecule counts
#' @param ERCCConcentration A vector of concentrations, in moles per litre,
#' of the spike-in genes in the mixture used.
#' @param DilutionFactor The dilution factor applied to the spike-in mixture.
#' e.g., 1μl per 50ml would be a 1/50000 \code{DilutionFactor}.
#' @param VolumePerCell The volume of spike-in mixture added to each well, or to
#' each cell.
#' @return The molecule counts per well, or per cell, based on the input
#' parameters.
BASiCS_CalculateSpikeIns <- function(
        ERCCConcentration,
        DilutionFactor,
        VolumePerCell
    ) {
    ## Molecule count per L
    ## (1 mole comprises 6.02214076 x 10^{23} molecules)
    ERCC_countmul <- ERCC_mmul * (6.02214076 * (1e23))
    ## Application of the dilution factor (1:50,000)
    ERCC_count <- ERCC_countmul * DilutionFactor
    ## Multiplying by the volume added into each well
    ERCC_count_final <- ERCC_count * VolumePerCell
    ERCC_count_final
}
