HiddenBASiCS_MCMC_ExtraArgs <- function(Args, Data, Burn, n)
{
  # Whether the regression case is used
  Regression <- ifelse("Regression" %in% names(Args), Args$Regression, FALSE)
  if (!("PriorDelta" %in% names(Args)) & Regression == FALSE) {
    message("-------------------------------------------------------------\n", 
            "NOTE: default choice PriorDelta = 'log-normal'  (recommended value). \n",
            "Vallejos et al (2015) used a 'gamma' prior instead.\n", 
            "-------------------------------------------------------------\n")
  }
  if(Regression == TRUE) {
    k <- ifelse("k" %in% names(Args), Args$k, 12)
    if(k <= 3) stop("The number of basis functions needs to be >= 4.")
    variance <- ifelse("Var" %in% names(Args), Args$Var, 1.2)
  } 
  else {
    k <- NULL; variance <- NULL; eta <- NULL
  }
  
  # Whether spike-ins are in use
  WithSpikes <- ifelse("WithSpikes" %in% names(Args), Args$WithSpikes, TRUE)
  # Whether a sthocastic reference is used (no spikes case only)
  if(WithSpikes == FALSE)
  {
    StochasticRef <- ifelse("StochasticRef" %in% names(Args), 
                            Args$StochasticRef, TRUE)  
    ConstrainType <- ifelse("ConstrainType" %in% names(Args), 
                            Args$ConstrainType, 1)
  } else { StochasticRef <- FALSE; ConstrainType <- NULL }
  
  # MCMC and storage parameters
  AR <- ifelse("AR" %in% names(Args), Args$AR, 0.44)
  StopAdapt <- ifelse("StopAdapt" %in% names(Args), Args$StopAdapt, Burn)
  StoreChains <- ifelse("StoreChains" %in% names(Args), Args$StoreChains, FALSE)
  StoreAdapt <- ifelse("StoreAdapt" %in% names(Args), Args$StoreAdapt, FALSE)
  StoreDir <- ifelse("StoreDir" %in% names(Args), Args$StoreDir, getwd())
  RunName <- ifelse("RunName" %in% names(Args), Args$RunName, "")
  PrintProgress <- ifelse("PrintProgress" %in% names(Args), Args$PrintProgress, TRUE)
  
  # Prior parameters
  if ("PriorParam" %in% names(Args)) { PriorParam <- Args$PriorParam } 
  else 
  {
    PriorParam <- list(s2.mu = 0.5, s2.delta = 0.5, a.delta = 1, 
                       b.delta = 1, p.phi = rep(1, times = n), 
                       a.s = 1, b.s = 1, a.theta = 1, b.theta = 1)
    if(Regression == TRUE) {
      PriorParam$m <- rep(0, k); PriorParam$V <- diag(k) 
      PriorParam$a.sigma2 <- 2; PriorParam$b.sigma2 <- 2  
      PriorParam$eta <- ifelse("eta" %in% names(Args), Args$eta, 5)
    }
  }
  PriorDelta <- ifelse("PriorDelta" %in% names(Args), Args$PriorDelta, "log-normal")
  PriorDeltaNum <- ifelse(PriorDelta == "gamma", 1, 2)
  # Validity checks
  if (!(PriorParam$s2.mu > 0 & length(PriorParam$s2.mu) == 1 & 
        PriorParam$s2.delta > 0 & length(PriorParam$s2.delta) == 1 & 
        PriorParam$a.delta > 0 & length(PriorParam$a.delta) == 1 & 
        PriorParam$b.delta > 0 & length(PriorParam$b.delta) == 1 & 
        all(PriorParam$p.phi > 0) & length(PriorParam$p.phi) == n & 
        PriorParam$a.s > 0 & length(PriorParam$a.s) == 1 & 
        PriorParam$b.s > 0 & length(PriorParam$b.s) == 1 & 
        PriorParam$a.theta > 0 & length(PriorParam$a.theta) == 1 & 
        PriorParam$b.theta > 0 & length(PriorParam$b.theta) == 1 )) 
    stop("Invalid prior hyper-parameter values.")
  if(Regression == TRUE) {
    if(!(length(PriorParam$m) == k & 
         ncol(PriorParam$V) == k & nrow(PriorParam$V) == k &
         PriorParam$a.sigma2 > 0 & length(PriorParam$a.sigma2) == 1 & 
         PriorParam$b.sigma2 > 0 & length(PriorParam$b.sigma2) == 1 ))
      stop("Invalid prior hyper-parameter values.")
  }
    
  # Starting values for MCMC
  if ("Start" %in% names(Args)) { Start <- Args$Start } 
  else { Start <- HiddenBASiCS_MCMC_Start(Data, PriorParam, WithSpikes) }
  
  if (!(AR > 0 & AR < 1 & length(AR) == 1)) 
    stop("Invalid AR value. Recommended value: AR = 0.44.")
  if (!(StopAdapt > 0)) 
    stop("Invalid StopAdapt value.")
  if (!(is.logical(StoreChains) & length(StoreChains) == 1)) 
    stop("Invalid StoreChains value.")
  if (!(is.logical(StoreAdapt) & length(StoreAdapt) == 1)) 
    stop("Invalid StoreAdapt value.")
  if (!(file.info(StoreDir)["isdir"])) 
    stop("Invalid StoreDir value.")
  if (!(PriorDelta %in% c("gamma", "log-normal"))) 
    stop("Invalid PriorDelta value.")

 list(AR = AR, StopAdapt = StopAdapt, StoreChains = StoreChains,
      StoreAdapt = StoreAdapt, StoreDir = StoreDir, 
      RunName = RunName, PrintProgress = PrintProgress, 
      PriorParam = PriorParam, PriorDeltaNum = PriorDeltaNum, 
      PriorDelta = PriorDelta,
      WithSpikes = WithSpikes, StochasticRef = StochasticRef,
      ConstrainType = ConstrainType, 
      Regression = Regression, k = k, variance = variance, 
      Start = Start)
}