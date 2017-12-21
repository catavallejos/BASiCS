HiddenBASiCS_MCMC_InputCheck <- function(Data, N, Thin, Burn)
{
  if (!is(Data, "SingleCellExperiment")) 
    stop("'Data' is not a SingleCellExperiment class object.")
  # The following checks are only relevant when the input data was
  # not created using the `newBASiCS_Data` function
  if(!("SpikeInput" %in% names(metadata(Data))))
    stop("'Data' does not contained all the required information \n",
         "See: https://github.com/catavallejos/BASiCS/wiki/2.-Input-preparation")
  if(!("BatchInfo" %in% names(metadata(Data))))
    stop("'Data' does not contained all the required information \n",
         "See: https://github.com/catavallejos/BASiCS/wiki/2.-Input-preparation")
  errors <- HiddenChecksBASiCS_Data(assay(Data), isSpike(Data), 
                                    metadata(Data)$SpikeInput, 
                                    rownames(assay(Data)), 
                                    metadata(Data)$BatchInfo)
  if (length(errors) > 0) stop(errors) 
  
  if (!(length(N) == 1 | length(Thin) == 1 | length(Burn) == 1)) 
    stop("Invalid parameter values.")
  if (!(N%%Thin == 0 & N >= max(4, Thin))) 
    stop("Please use an integer value for N (N>=4); multiple of thin.")
  if (!(Thin%%1 == 0 & Thin >= 2)) 
    stop("Please use an integer value for Thin (Thin>=2).")
  if (!(Burn%%Thin == 0 & Burn < N & Burn >= 1)) 
    stop("Please use an integer value for Burn (1<=Burn<N); multiple of thin.")
}