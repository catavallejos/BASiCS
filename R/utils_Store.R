.BASiCS_MCMC_RefFreqStore <- function(Data,
                                      Chain, 
                                      RefGene,
                                      RefGenes,
                                      ConstrainType, 
                                      StoreDir,
                                      RunName) {
  
  # Directory
  OldDir <- getwd()
  setwd(StoreDir)
  ## in case function exits early/fails
  on.exit(setwd(StoreDir))
  
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
  
}

.BASiCS_MCMC_OutputStore <- function(ChainClass,
                                     Chain,
                                     StoreChains,
                                     StoreAdapt,
                                     StoreDir,
                                     RunName) {
  # Directory
  OldDir <- getwd()
  setwd(StoreDir)
  ## in case function exits early/fails
  on.exit(setwd(OldDir))

  # Store MCMC chain
  if (StoreChains) {
    message("-------------------------------------------------------------\n",
            "BASiCS_Chain object stored as ",
            paste0("chain_", RunName, ".Rds"), " file in", "\n",
            paste0("'", StoreDir, "' directory ... "), "\n",
            "-------------------------------------------------------------\n")
    saveRDS(ChainClass, file = paste0("chain_", RunName, ".Rds"))
  }

  # Store adaptive variances
  if (StoreAdapt) {
    message("-------------------------------------------------------------\n",
            "Storing trajectories of adaptive proposal variances (log-scale) as",
            "chain_ls_", RunName, ".Rds file in \n",
            "'", StoreDir, "' directory ... \n",
            "-------------------------------------------------------------\n")
    ChainLS <- list(ls.mu = Chain$ls.mu, ls.delta = Chain$ls.delta,
                    ls.nu = Chain$ls.nu, ls.theta = Chain$ls.theta)
    if ("ls.phi" %in% names(Chain)) {
      ChainLS$ls.phi <- Chain$ls.phi
    }

    saveRDS(ChainLS, file = paste0("chain_ls_", RunName, ".Rds"))
  }
}
