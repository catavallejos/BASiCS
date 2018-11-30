## ---- echo=FALSE, results="hide", message=FALSE----------------------------
require(knitr)
opts_chunk$set(error=FALSE, message=FALSE, warning=FALSE, dpi = 30)
knitr::opts_chunk$set(dev="png")

## ----library, echo=FALSE---------------------------------------------------
library(BASiCS)
library(BiocStyle)

## ----quick-start-MCMC------------------------------------------------------
Data <- makeExampleBASiCS_Data()
Chain <- BASiCS_MCMC(Data = Data, N = 1000, Thin = 10, Burn = 500, 
                     PrintProgress = FALSE, Regression = TRUE)

## ----ExampleDataTest-------------------------------------------------------
set.seed(1)
Counts <- matrix(rpois(50*40, 2), ncol = 40)
rownames(Counts) <- c(paste0("Gene", 1:40), paste0("Spike", 1:10))
Tech <- c(rep(FALSE,40),rep(TRUE,10))
set.seed(2)
SpikeInput <- rgamma(10,1,1)
SpikeInfo <- data.frame("SpikeID" = paste0("Spike", 1:10), 
                        "SpikeInput" = SpikeInput)

# No batch structure
DataExample <- newBASiCS_Data(Counts, Tech, SpikeInfo)

# With batch structure
DataExample <- newBASiCS_Data(Counts, Tech, SpikeInfo, 
                              BatchInfo = rep(c(1,2), each = 20)) 

## ----ExampleDataNoSpikes---------------------------------------------------
set.seed(1)
CountsNoSpikes <- matrix(rpois(50*40, 2), ncol = 40)
rownames(CountsNoSpikes) <- paste0("Gene", 1:50)

# With batch structure
DataExampleNoSpikes <- SingleCellExperiment(assays = list(counts = CountsNoSpikes), 
                      colData = data.frame(BatchInfo = rep(c(1,2), each = 20))) 

## ----LoadSingleData--------------------------------------------------------
data(ChainSC)

## ----quick-start-HVGdetection, fig.height = 20, fig.width = 20-------------
par(mfrow = c(2,2))
HVG <- BASiCS_DetectHVG(ChainSC, VarThreshold = 0.6, Plot = TRUE)
LVG <- BASiCS_DetectLVG(ChainSC, VarThreshold = 0.2, Plot = TRUE)

## ----quick-start-HVGdetectionTable-----------------------------------------
head(HVG$Table)
head(LVG$Table)

## ----quick-start-HVGdetectionPlot, fig.width=8, fig.height=8---------------
SummarySC <- Summary(ChainSC)
plot(SummarySC, Param = "mu", Param2 = "delta", log = "xy")
with(HVG$Table[HVG$Table$HVG,], points(Mu, Delta))
with(LVG$Table[LVG$Table$LVG,], points(Mu, Delta))

## ----quick-start-LoadBothData----------------------------------------------
data(ChainSC)
data(ChainRNA)

## ----quick-start-TestDE, fig.width=16, fig.height=8------------------------
Test <- BASiCS_TestDE(Chain1 = ChainSC, Chain2 = ChainRNA,
                      GroupLabel1 = "SC", GroupLabel2 = "PaS",
                      EpsilonM = log2(1.5), EpsilonD = log2(1.5),
                      EFDR_M = 0.10, EFDR_D = 0.10,
                      Offset = TRUE, PlotOffset = FALSE, Plot = TRUE)

## ----quick-start-DE--------------------------------------------------------
head(Test$TableMean)
head(Test$TableDisp)

## ----quick-start-TestDE-2, fig.width=16, fig.height=8----------------------
Test <- BASiCS_TestDE(Chain1 = ChainSC, Chain2 = ChainRNA,
                      GroupLabel1 = "SC", GroupLabel2 = "PaS",
                      EpsilonM = 0, EpsilonD = log2(1.5),
                      EFDR_M = 0.10, EFDR_D = 0.10,
                      Offset = TRUE, PlotOffset = FALSE, Plot = FALSE)

## ---- fig.height = 8, fig.width = 8----------------------------------------
DataNoSpikes <- newBASiCS_Data(Counts, Tech, SpikeInfo = NULL, 
                              BatchInfo = rep(c(1,2), each = 20))

# Alternatively
DataNoSpikes <- SingleCellExperiment(assays = list(counts = Counts), 
                      colData = data.frame(BatchInfo = rep(c(1,2), each = 20))) 

ChainNoSpikes <- BASiCS_MCMC(Data = DataNoSpikes, N = 1000, 
                             Thin = 10, Burn = 500, 
                             WithSpikes = FALSE,  Regression = TRUE,
                             PrintProgress = FALSE)

## ---- fig.height = 8, fig.width = 8----------------------------------------
DataRegression <- newBASiCS_Data(Counts, Tech, SpikeInfo, 
                              BatchInfo = rep(c(1,2), each = 20))
ChainRegression <- BASiCS_MCMC(Data = DataRegression, N = 1000, 
                               Thin = 10, Burn = 500, 
                               Regression = TRUE, 
                               PrintProgress = FALSE)

## ---- fig.height = 8, fig.width = 8----------------------------------------
data("ChainRNAReg")
BASiCS_showFit(ChainRNAReg)

## ---- fig.height = 8, fig.width = 16---------------------------------------
data("ChainSCReg")
Test <- BASiCS_TestDE(Chain1 = ChainSCReg, Chain2 = ChainRNAReg,
                      GroupLabel1 = "SC", GroupLabel2 = "PaS",
                      EpsilonM = log2(1.5), EpsilonD = log2(1.5),
                      EpsilonR = log2(1.5)/log2(exp(1)),
                      EFDR_M = 0.10, EFDR_D = 0.10,
                      Offset = TRUE, PlotOffset = FALSE, Plot = FALSE)

## --------------------------------------------------------------------------
head(Test$TableResDisp, n = 2)

## ----MCMCNotRun------------------------------------------------------------
Data <- makeExampleBASiCS_Data()
Chain <- BASiCS_MCMC(Data, N = 1000, Thin = 10, Burn = 500, Regression = TRUE,
                     PrintProgress = FALSE, StoreChains = TRUE, 
                     StoreDir = tempdir(), RunName = "Example")

## ----LoadChainNotRun-------------------------------------------------------
Chain <- BASiCS_LoadChain("Example", StoreDir = tempdir()) 

## ----Traceplots, fig.height = 8, fig.width = 16----------------------------
plot(Chain, Param = "mu", Gene = 1, log = "y")

## ----AccessChains----------------------------------------------------------
displayChainBASiCS(Chain, Param = "mu")[1:2,1:2]

## ----Summary---------------------------------------------------------------
ChainSummary <- Summary(Chain)

## ----SummaryExample--------------------------------------------------------
head(displaySummaryBASiCS(ChainSummary, Param = "mu"), n = 2)

## ----OtherHPD, fig.width = 16, fig.height = 8------------------------------
par(mfrow = c(1,2))
plot(ChainSummary, Param = "mu", main = "All genes", log = "y")
plot(ChainSummary, Param = "delta", 
     Genes = c(2,5,10,50), main = "5 customized genes")

## ----Normalisation, fig.width = 16, fig.height = 8-------------------------
par(mfrow = c(1,2))
plot(ChainSummary, Param = "phi")
plot(ChainSummary, Param = "s", Cells = 1:5)

## ----DispVsExp, fig.width = 16, fig.height = 8-----------------------------
par(mfrow = c(1,2))
plot(ChainSummary, Param = "mu", Param2 = "delta", log = "x", SmoothPlot = FALSE)
plot(ChainSummary, Param = "mu", Param2 = "delta", log = "x", SmoothPlot = TRUE)

## ----DenoisedCounts--------------------------------------------------------
DenoisedCounts <- BASiCS_DenoisedCounts(Data = Data, Chain = Chain)
DenoisedCounts[1:2, 1:2]

## ----DenoisedProp----------------------------------------------------------
DenoisedRates <- BASiCS_DenoisedRates(Data = Data, Chain = Chain, 
                                      Propensities = FALSE)
DenoisedRates[1:2, 1:2]

## ----DenoisedRates---------------------------------------------------------
DenoisedProp <- BASiCS_DenoisedRates(Data = Data, Chain = Chain, 
                                    Propensities = TRUE)
DenoisedProp[1:2, 1:2]

## ----SessionInfo-----------------------------------------------------------
sessionInfo()

