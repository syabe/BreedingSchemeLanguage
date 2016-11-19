#'Self-fertilize
#'
#'@param nProgeny the number of progeny
#'@param popID population ID to be self-fertilized (default: the latest population)
#'
#'@return sequence information of progenies and the all information created before (list)
#'
#'@export
selfFertilize <- function(nProgenyPerInd = 1, popID = NULL){
  selfFertilize.func <- function(data, nProgenyPerInd, popID){
    mapData <- data$mapData
    breedingData <- data$breedingData
    score <- data$score
    selCriterion <- data$selCriterion
    if(is.null(popID)){
      popID <- max(breedingData$popID)
    }
    tf <- breedingData$popID %n% popID
    GID.now <- breedingData$GID[tf]
    geno.now <- breedingData$geno[sort(c(GID.now * 2 - 1, GID.now * 2)), ]
    geno.progeny <- makeSelfs(popSize = nProgeny, geno = geno.now, pos = mapData$map$Pos)$progenies
    gValue <- calcGenotypicValue(geno = geno.progeny, mapData = mapData)
    GID.progeny <- max(breedingData$GID) + 1:(nrow(geno.progeny) / 2)
    popID.progeny <- rep(max(breedingData$popID) + 1, nrow(geno.progeny) / 2)
    breedingData$geno <- rbind(breedingData$geno, geno.progeny)
    breedingData$GID <- c(breedingData$GID, GID.progeny)
    breedingData$popID <- c(breedingData$popID, popID.progeny)
    breedingData$popIDsel <- c(breedingData$popIDsel, popID.progeny)
    breedingData$gValue <- c(breedingData$gValue, gValue)
    return(list(mapData = mapData, breedingData = breedingData, score = score, selCriterion = selCriterion))
  }
  if(nCore > 1){
    sfInit(parallel=T, cpus=nCore)
    lists <<- sfLapply(lists, selfFertilize.func, nProgenyPerInd = nProgenyPerInd, popID = popID)
    sfStop()
  }else{
    lists <<- lapply(lists, selfFertilize.func, nProgenyPerInd = nProgenyPerInd, popID = popID)
  }
}
