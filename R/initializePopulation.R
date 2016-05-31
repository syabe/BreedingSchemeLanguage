#'Create a founder population
#'
#'@param nPop population size
#'@param gVariance genetic variance in the initial population
#'
#'@return initial population informationand the all information created before (list)
#'
#'@export
initializePopulation <- function(nPop = 100, gVariance = 1){
  initializePopulation.func <- function(data, nPop, gVariance){
    mapData <- data$mapData
    founderHaps <- data$founderHaps
    seed <- round(runif(1, 0, 1e9))
    doubleGametes <- function(gametes){
      genotypes <- matrix(NA, 2 * nrow(gametes), ncol(gametes))
      genotypes[1:nrow(gametes) * 2, ] <- gametes
      genotypes[1:nrow(gametes) * 2 - 1, ] <- gametes
    }
    geno <- doubleGametes(founderHaps)
    geno <- randomMate(popSize = nrow(geno) / 2, geno = geno, pos = mapData$map$Pos)$progenies
    geno <- randomMate(popSize = nPop, geno = geno, pos = mapData$map$Pos)$progenies
    geno <- geno * 2 - 1
    g <- calcGenotypicValue(geno = geno, mapData = mapData)
    vg <- var(g)
    vg.ideal <- gVariance
    coef <- sqrt(vg.ideal / vg)
    mapData$effects <- as.matrix(mapData$effects * coef, ncol = 1)
    GID <- 1:nPop
    popID <- rep(0, nPop)
    gValue <- calcGenotypicValue(geno = geno, mapData = mapData)
    breedingData <- list(geno = geno, GID = GID, popID = popID, popIDsel = popID, gValue = gValue)
    return(list(mapData = mapData, breedingData = breedingData, score = NULL, selCriterion = NULL))
  }
  if(nCore > 1){
    sfInit(parallel=T, cpus=nCore)
    lists <<- sfLapply(lists, initializePopulation.func, nPop = nPop, gVariance = gVariance)
    sfStop()
  }else{
    lists <<- lapply(lists, initializePopulation.func, nPop = nPop, gVariance = gVariance)
  }
}
