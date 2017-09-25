HiddenEFDR <- function(EviThreshold, VarThreshold, Prob) {
    return(sum((1 - Prob) * I(Prob > EviThreshold))/sum(I(Prob > EviThreshold)))
}

HiddenEFNR <- function(EviThreshold, VarThreshold, Prob) {
    return(sum(Prob * I(EviThreshold >= Prob))/sum(I(EviThreshold >= Prob)))
}





