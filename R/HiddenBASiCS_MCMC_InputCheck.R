HiddenBASiCS_MCMC_InputCheck <- function(Data, N, Thin,
                                         Burn, Regression, WithSpikes)
{

  if (!is(Data, "SingleCellExperiment")) 
    stop("'Data' is not a SingleCellExperiment class object.\n")
  
  # The following checks are only relevant when the input data was
  # not created using the `newBASiCS_Data` function

  # Check if `Data` contains more than one type of types
  if( (length(SingleCellExperiment::altExpNames(Data)) > 1))
    stop("'Data' contains multiple 'altExp' assays but only one is allowed. \n",
         "The latter should include information for spike-in genes. \n",
         "Please remove unwanted the remaining elements. See help(altExp). \n")
  
  if(WithSpikes) {
    message("altExp '", altExpNames(Data),"' is assumed to contain spike-in genes.\n",
            "see help(altExp) for details.")  
  }
  
  # If SpikeInput slot is missing and WithSpikes == TRUE
  if((length(altExpNames(Data)) > 0) & WithSpikes & 
     is.null(metadata(Data)$SpikeInput))
    stop("'Data' does not contain the 'SpikeInput' slot. \n",
         "See: https://github.com/catavallejos/BASiCS/wiki/2.-Input-preparation\n")
  
  # If isSpike slot is missing and WithSpikes == TRUE
  if((length(altExpNames(Data)) == 0)  & WithSpikes)
    stop("'Data' does not contain information about spike-in fenes \n", 
         "Please indicate include this information using 'altExp' \n",
         "or set 'WithSpikes = FALSE' \n.",
         "See: https://github.com/catavallejos/BASiCS/wiki/2.-Input-preparation\n")
  
 # # If isSpike slot does not contain spikes and WithSpikes == TRUE
 # if(sum(Spikes) == 0 & WithSpikes)
 #   stop("'isSpike(Data)' does not contain TRUE values, meaning the sce object \n",
 #        "does not contain spikes. Please indicate in 'isSpike(Data)' which \n",
 #        "genes are spike-ins or set 'WithSpikes = FALSE' \n.",
 #        "See: https://github.com/catavallejos/BASiCS/wiki/2.-Input-preparation\n")
  
  # If BatchInfo slot is missing and WithSpikes == FALSE
  if(!WithSpikes & is.null(colData(Data)$BatchInfo))
    stop("'Data' does not contain a BatchInfo vector needed when 'WithSpikes = FALSE'. \n", 
         "Please assign the batch information to: 'colData(Data)$BatchInfo = BatchInfo'. \n",
         "See: https://github.com/catavallejos/BASiCS/wiki/2.-Input-preparation\n")

  # Checking how counts are stored
  if(!("counts" %in% assayNames(Data)))
    stop("'Data' does not contain a 'counts' slot. \n",
         "Please make sure to include the raw data in the \n", 
         "SingleCellExperiment object under the name 'counts'. \n",
         "See: https://github.com/catavallejos/BASiCS/wiki/2.-Input-preparation\n")
  
  # Further checks on input
  errors <- HiddenChecksBASiCS_Data(Data, WithSpikes)
  
  if (length(errors) > 0) stop(errors) 
  
  if (!(length(N) == 1 | length(Thin) == 1 | length(Burn) == 1)) 
    stop("Invalid parameter values.\n")
  if (!(N%%Thin == 0 & N >= max(4, Thin))) 
    stop("Please use an integer value for N (N>=4); multiple of thin.\n")
  if (!(Thin%%1 == 0 & Thin >= 2)) 
    stop("Please use an integer value for Thin (Thin>=2).\n")
  if (!(Burn%%Thin == 0 & Burn < N & Burn >= 1)) 
    stop("Please use an integer value for Burn (1<=Burn<N); multiple of thin.\n")
  if (!is.logical(Regression)) 
    stop("Please use a logical value for the Regression parameter.\n")  
  
}
