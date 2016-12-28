#'Define and create species
#'
#'@param loadData if null create a new species (default), else the file name of previously created species like "fileName.RData".
#'@param importFounderHap if null create new founder haplotypes (default), else the file name of externally generated founder haplotypes in hapmap format.
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
defineSpecies <- function(loadData=NULL, importFounderHap=NULL, saveDataFileName="previousData", nSim=1, nCore=1, nChr=7, lengthChr=150, effPopSize=100, nMarkers=1000, nQTL=50, propDomi=0, nEpiLoci=0){
  defineSpecies.func <- function(simNum, nChr, lengthChr, effPopSize, nMarkers, nQTL, propDomi, nEpiLoci, founderHaps=NULL){
    seed <- round(runif(1, 0, 1e9))
    nLoci <- nMarkers + nQTL * (nEpiLoci + 1) * 2
    if (is.null(founderHaps)){
      minMAF <- 0.01
      piecesPerM <- 10000
      nPiecesPerChr <- lengthChr / 100 * piecesPerM
      recBTpieces <- 1 / piecesPerM
      coalSim <- getCoalescentSim(effPopSize=2 * effPopSize, nMrkOrMut=nLoci, nChr=nChr, nPiecesPerChr=nPiecesPerChr, recBTpieces=recBTpieces, minMAF=minMAF, seed=seed)
      markers <- coalSim$markers
      map <- coalSim$map
    } else{
      markers <- founderHaps$markers
      map <- founderHaps$map
      if (nrow(map) < nLoci) print("Not enough loci in imported founder haplotypes for both markers and QTL")
      # This function will not allow missing data so anything missing will be
      # set to the major allele
      markers[is.na(markers)] <- 1
    }
    ancestralState <- rbinom(nLoci, 1, 0.5)
    markers[,ancestralState == 1] <- 1 - markers[,ancestralState == 1]
    mapData <- makeMap(map=map, nLoci=nLoci, nMarkers=nMarkers, nQTL=nQTL, propDomi=propDomi, interactionMean=nEpiLoci)
    return(list(mapData=mapData, founderHaps=markers))
  }#END defineSpecies.func
  
  if(is.null(loadData)){
    if (is.null(importFounderHap)){
    sims <- lapply(1:nSim, defineSpecies.func, nChr=nChr, lengthChr=lengthChr, effPopSize=effPopSize, nMarkers=nMarkers, nQTL=nQTL, propDomi=propDomi, nEpiLoci=nEpiLoci)
    }else{
      foundHap <- read.table(file=paste(importFounderHap, ".hmp", sep=""))
      foundHap <- phasedHapMap2mat(foundHap)
      sims <- lapply(1:nSim, defineSpecies.func, nChr=nChr, lengthChr=lengthChr, effPopSize=effPopSize, nMarkers=nMarkers, nQTL=nQTL, propDomi=propDomi, nEpiLoci=nEpiLoci, founderHaps=foundHap)
    }
    save(sims, nSim, nCore, file=paste(saveDataFileName, ".RData", sep=""))
  }else{
    load(paste(loadData, ".RData", sep=""))
  }
  # list of objects to remove before returning the environment
  toRemove <- c(setdiff(ls(), c("sims", "nSim", "nCore")), "toRemove")
  rm(list=toRemove)
  return(environment())
}
