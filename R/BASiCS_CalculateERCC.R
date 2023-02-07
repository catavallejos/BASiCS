#' Convert concentration in moles per microlitre to molecule counts
#' @param Mix The name of the spike-in mix to use.
#' @param DilutionFactor The dilution factor applied to the spike-in mixture.
#' e.g., 1 microlitre per 50ml would be a 1/50000 \code{DilutionFactor}.
#' @param VolumePerCell The volume of spike-in mixture added to each well, or to
#' each cell.
#' @return The molecule counts per well, or per cell, based on the input
#' parameters.
#' @export
BASiCS_CalculateERCC <- function(
        Mix,
        DilutionFactor,
        VolumePerCell
    ) {
    if (!length(Mix) == 1 && Mix %in% c(1, 2)) {
        stop("Invalid value for Mix.")
    }
    if (!requireNamespace("scRNAseq", quietly = TRUE)) {
        stop("scRNAseq package is required to retrieve ERCC mix concentrations")
    }
    Col <- sprintf("concentration in Mix %s (attomoles/ul)", Mix)

    # Moles per micro litre
    ERCC_mmul <- scRNAseq::ERCCSpikeInConcentrations()[[Col]] * (10^(-18))
    ## Molecule count per L
    ## (1 mole comprises 6.02214076 x 10^{23} molecules)
    ERCC_countmul <- ERCC_mmul * (6.02214076 * (1e23))
    ## Application of the dilution factor (1:50,000)
    ERCC_count <- ERCC_countmul * DilutionFactor
    ## Multiplying by the volume added into each well
    ERCC_count_final <- ERCC_count * VolumePerCell
    ERCC_count_final
}
