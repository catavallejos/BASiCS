HiddenBASiCS_MCMC_RefFreqStore <- function(Data, Chain, 
                                           RefGene, RefGenes, ConstrainType, 
                                           StoreDir, RunName)
{
  # Directory
  OldDir <- getwd()
  setwd(StoreDir)
  
  message("-------------------------------------------------------------\n", 
          paste("BASiCS version", packageVersion("BASiCS"), 
                ": horizontal integration (no-spikes case)"), "\n", 
          "-------------------------------------------------------------\n", 
          "ConstrainType: ", ConstrainType, "\n")
  
  if (length(RefGenes) == 1) {
    
    message("Reference gene:", RefGene + 1, "\n", 
            "Information stored as a .txt file in \n", 
            "'", StoreDir, "' directory ... \n", 
            "-------------------------------------------------------------\n")
    
    TableRef <- cbind.data.frame(GeneNames = rownames(counts(Data))[RefGene + 1], 
                                 GeneIndex = RefGene + 1, 
                                 stringsAsFactors = FALSE)
    write.table(TableRef, paste0("TableRef_", RunName, ".txt"), 
                col.names = TRUE, row.names = FALSE)
  } 
  else {
    
    message("Randomly, 1 out of ", length(RefGenes), "\n",
            "genes was left as reference at each iteration \n", 
            "List of reference genes and frequencies stored as a .txt file in\n", 
            "'", StoreDir, "' directory ... \n", 
            "-------------------------------------------------------------\n")

    TableRef <- cbind.data.frame(GeneNames = rownames(counts(Data))[RefGenes + 1], 
                                GeneIndex = RefGenes + 1, 
                                ReferenceFreq = Chain$RefFreq[RefGenes + 1], 
                                stringsAsFactors = FALSE)
    write.table(TableRef, paste0("TableRef_", RunName, ".txt"), 
                col.names = TRUE, row.names = FALSE)
  }
  
  # Restore working directory
  setwd(OldDir)
}