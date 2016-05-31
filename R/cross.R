#'Random mate
#'
#'@param nProgeny the number of progenies
#'@param equalContribution if T all individuals used the same number of times as parents, if F individuals chosen at random to be parents
#'@param popID population ID to be crossed (default: the latest population)
#'
#'@return sequence information of progenies and the all information created before (list)
#'
#'@export
cross <- function(nProgeny = 100, equalContribution = F, popID = NULL){
  cross.func <- function(data, nProgeny, equalContribution, popID){
    mapData <- data$mapData
    breedingData <- data$breedingData
    score <- data$score
    selCriterion <- data$selCriterion
    if(is.null(popID)){
      popID <- max(breedingData$popID)
    }
    tf <- rep(F, length(breedingData$GID))
    for(i in popID){
      tf[breedingData$popID == i] <- T
    }
    GID.now <- breedingData$GID[tf]
    geno.now <- breedingData$geno[sort(c(GID.now * 2 - 1, GID.now * 2)), ]
    if(equalContribution){
      geno.progeny <- randomMateAll(popSize = nProgeny, geno = geno.now, pos = mapData$map$Pos)$progenies
    }else{
      geno.progeny <- randomMate(popSize = nProgeny, geno = geno.now, pos = mapData$map$Pos)$progenies
    }
    gValue <- calcGenotypicValue(geno = geno.progeny, mapData = mapData)
    GID.progeny <- max(breedingData$GID) + 1:nProgeny
    popID.progeny <- rep(max(breedingData$popID) + 1, nProgeny)
    breedingData$geno <- rbind(breedingData$geno, geno.progeny)
    breedingData$GID <- c(breedingData$GID, GID.progeny)
    breedingData$popID <- c(breedingData$popID, popID.progeny)
    breedingData$popIDsel <- c(breedingData$popIDsel, popID.progeny)
    breedingData$gValue <- c(breedingData$gValue, gValue)
    return(list(mapData = mapData, breedingData = breedingData, score = score, selCriterion = selCriterion))
  }
  if(nCore > 1){
    sfInit(parallel=T, cpus=nCore)
    lists <<- sfLapply(lists, cross.func, nProgeny = nProgeny, equalContribution = equalContribution, popID = popID)
    sfStop()
  }else{
    lists <<- lapply(lists, cross.func, nProgeny = nProgeny, equalContribution = equalContribution, popID = popID)
  }
}
