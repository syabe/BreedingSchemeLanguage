#'Define and create species
#'
#'@param loadData if null create a new species (default), else the file name of previously created species like "fileName.RData".
#'@param saveDataFileName the name to save the species data with double-quotation, like "species1_1". (default: "previousData")
#'@param nSim the number of simulation trials
#'@param nCore simulation processed in parallel over this number of CPUs (Check computer capacity before setting above 1.)
#'@param nChr the number of chromosomes
#'@param lengthChr the length of each chromosome (cM; all chromosomes are the same length)
#'@param effPopSize the effective population size in the base population
#'@param nMarkers the number of markers, which is used especially for genomic selection
#'@param nQTL the number of QTLs controlling the target trait
#'@param propDomi the probability of dominant QTL among the all QTL
#'@param nEpiLoci the expected number of epistatic loci for each effect
#'
#'@return Species information and input values for the simulation (list)
#'
#'@export
defineSpecies <- function(loadData = NULL, saveDataFileName = "previousData", nSim = 1, nCore = 1, nChr = 7, lengthChr = 150, effPopSize = 100, nMarkers = 1000, nQTL = 50, propDomi = 0, nEpiLoci = 0){
  defineSpecies.func <- function(data, nChr, lengthChr, effPopSize, nMarkers, nQTL, propDomi, nEpiLoci){
    seed <- round(runif(1, 0, 1e9))
    nLoci <- nMarkers + nQTL * (nEpiLoci + 1) * 2
    minMAF <- 0.01
    piecesPerM <- 10000
    nPiecesPerChr <- lengthChr / 100 * piecesPerM
    recBTpieces <- 1 / piecesPerM
    coalSim <- getCoalescentSim(effPopSize = 2 * effPopSize, nMrkOrMut = nLoci, nChr = nChr, nPiecesPerChr = nPiecesPerChr, recBTpieces = recBTpieces, minMAF = minMAF, seed = seed)
    markers <- coalSim$markers
    map <- coalSim$map
    ancestralState <- rbinom(nLoci, 1, 0.5)
    markers[,ancestralState == 1] <- 1 - markers[,ancestralState == 1]
    mapData <- makeMap(map = map, nLoci = nLoci, nMarkers = nMarkers, nQTL = nQTL, propDomi = propDomi, interactionMean = nEpiLoci)
    return(list(mapData = mapData, founderHaps = markers))
  }
  if(is.null(loadData)){
    lists0 <- list()
    for(sim in 1:nSim){
      lists0[[sim]] <- 0
    }
    lists <<- lapply(lists0, defineSpecies.func, nChr = nChr, lengthChr = lengthChr, effPopSize = effPopSize, nMarkers = nMarkers, nQTL = nQTL, propDomi = propDomi, nEpiLoci = nEpiLoci)
    nSim <<- nSim
    nCore <<- nCore
    save(lists, nSim, nCore, file = paste(saveDataFileName, ".RData", sep = ""))
  }else{
    load(paste(loadData, ".RData", sep = ""))
    lists <<- lists
    nSim <<- nSim
    nCore <<- nCore
  }
}
