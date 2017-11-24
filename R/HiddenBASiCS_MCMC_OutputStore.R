HiddenBASiCS_MCMC_OutputStore <- function(ChainClass, Chain, 
                                          StoreChains, StoreAdapt, 
                                          StoreDir, RunName)
{
  # Directory
  OldDir <- getwd()
  setwd(StoreDir)
  
  # Store MCMC chain
  if(StoreChains) {
    message("-------------------------------------------------------------\n", 
            "BASiCS_Chain object stored as ", 
            paste0("chain_", RunName, ".Rds"), " file in", "\n", 
            paste0("'", StoreDir, "' directory ... "), "\n", 
            "-------------------------------------------------------------\n")
    saveRDS(ChainClass, file = paste0("chain_", RunName, ".Rds"))
  }
  
  # Store adaptive variances
  if(StoreAdapt) {
    message("-------------------------------------------------------------\n", 
            "Storing trajectories of adaptive proposal variances (log-scale) as", 
            "chain_ls_", RunName, ".Rds file in \n", 
            "'", StoreDir, "' directory ... \n", 
            "-------------------------------------------------------------\n")
    ChainLS <- list(ls.mu = Chain$ls.mu, ls.delta = Chain$ls.delta, 
                    ls.phi = Chain$ls.phi, ls.nu = Chain$ls.nu, 
                    ls.theta = Chain$ls.theta)
    saveRDS(ChainLS, file = paste0("chain_ls_", RunName, ".Rds"))
  }
  
  # Restore working directory
  setwd(OldDir)
}