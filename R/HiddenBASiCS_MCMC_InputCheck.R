HiddenBASiCS_MCMC_InputCheck <- function(Data, N, Thin, Burn, Regression)
{
  if (!is(Data, "SingleCellExperiment")) 
    stop("'Data' is not a SingleCellExperiment class object.")
  # The following checks are only relevant when the input data was
  # not created using the `newBASiCS_Data` function
  if(!("SpikeInput" %in% names(metadata(Data))))
    stop("'Data' does not contained all the required information \n",
         "See: https://github.com/catavallejos/BASiCS/wiki/2.-Input-preparation")
  if(!("BatchInfo" %in% names(colData(Data))))
    stop("'Data' does not contained all the required information \n",
         "See: https://github.com/catavallejos/BASiCS/wiki/2.-Input-preparation")
  if(!("counts" %in% names(assays(Data))))
    stop("'Data' does not contain a 'counts' slot. \n",
         "Please make sure to include the raw data in the SingleCellExperiment object under the name 'countss' \n",
         "See: https://github.com/catavallejos/BASiCS/wiki/2.-Input-preparation")
  errors <- HiddenChecksBASiCS_Data(counts(Data), isSpike(Data), 
                                    metadata(Data)$SpikeInput, 
                                    rownames(counts(Data)), 
                                    colData(Data)$BatchInfo)
  if (length(errors) > 0) stop(errors) 
  
  if (!(length(N) == 1 | length(Thin) == 1 | length(Burn) == 1)) 
    stop("Invalid parameter values.")
  if (!(N%%Thin == 0 & N >= max(4, Thin))) 
    stop("Please use an integer value for N (N>=4); multiple of thin.")
  if (!(Thin%%1 == 0 & Thin >= 2)) 
    stop("Please use an integer value for Thin (Thin>=2).")
  if (!(Burn%%Thin == 0 & Burn < N & Burn >= 1)) 
    stop("Please use an integer value for Burn (1<=Burn<N); multiple of thin.")
  if (!is.logical(Regression)) 
    stop("Please use a logical value for the Regression parameter.")  
}